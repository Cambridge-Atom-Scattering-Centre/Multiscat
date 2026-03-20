from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path
import tempfile
from typing import Any

import numpy as np
from slate_core import array
from scipy.constants import (  # type: ignore[import-untyped]
    angstrom,
    atomic_mass,
    electron_volt,
    physical_constants,
)
from slate_quantum import operator

from multiscat.basis import (
    close_coupling_basis,
    scattering_metadata_from_stacked_delta_x,
    split_scattering_metadata,
)
from multiscat.config import OptimizationConfig, ScatteringCondition
from multiscat.interpolate import ScatteringOperator

ROOT = Path(__file__).resolve().parents[1]
TESTS_DIR = Path(__file__).resolve().parent


def _parse_raw_intensities(output_file: Path) -> dict[tuple[int, int], float]:
    # Regex for lines without the '#' prefix: two ints and one float
    pattern = re.compile(r"^\s*(-?\d+)\s+(-?\d+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*")
    intensities: dict[tuple[int, int], float] = {}
    with output_file.open("r") as f:
        for line in f:
            stripped = line.strip()
            # Skip empty lines or actual comments
            if not stripped or stripped.startswith("#"):
                continue
                
            match = pattern.match(line)
            if match:
                h = int(match.group(1))
                k = int(match.group(2))
                val = float(match.group(3))
                intensities[(h, k)] = val
    return intensities

def _parse_intensities(output_file: Path) -> dict[tuple[int, int], float]:
    # Parse diffraction rows whether they are emitted as "h k I" or "# h k I".
    pattern = re.compile(
        r"^\s*#?\s*(-?\d+)\s+(-?\d+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*"
    )
    intensities: dict[tuple[int, int], float] = {}
    with output_file.open("r") as f:
        for line in f:
            match = pattern.match(line)
            if match:
                h = int(match.group(1))
                k = int(match.group(2))
                val = float(match.group(3))
                intensities[(h, k)] = val
    return intensities


def _condition_to_input_file(condition: ScatteringCondition) -> str:
    scattering_vector = np.asarray(condition.incident_k)
    scattering_magnitude = float(np.linalg.norm(scattering_vector))
    assert scattering_magnitude > 0, "Incident wavevector magnitude must be non-zero"

    mass_amu = condition.mass / atomic_mass
    energy_meV = condition.incident_energy / (electron_volt * 10**-3)
    theta_degrees = np.degrees(condition.theta)
    phi_degrees = np.degrees(condition.phi)

    scat_cond_lines = [
        "Comment line: Energy, theta, phi   This file defines the combination of conditions to be used by multiscat",
        f"{mass_amu:.10g}       !helium mass",
        f"{energy_meV:.10g},{theta_degrees:.10g},{phi_degrees:.10g}",
    ]
    return "\n".join(scat_cond_lines) + "\n"


def _optimization_to_input_file(config: OptimizationConfig) -> str:
    max_negative_energy = config.max_negative_energy / (electron_volt * 10**3) if config.max_negative_energy is not None else 120
    max_channel_index = config.max_channel_index if config.max_channel_index is not None else 120
    lines = [
        "1       !itest=1 enables output of each diffraction intensity; itest=0 outputs specular only",
        "0       !gmres preconditioner flag (ipc)",
        f"{int(np.log10(1 / config.precision))}       !number of significant figures convergence (nsf)",
        f"{max_negative_energy:.10g}       !max -ve energy of closed channels (dmax)",
        f"{max_channel_index}       !max index of channels (imax)",
    ]
    return "\n".join(lines) + "\n"


def _load_potential_file_as_array(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    data_start = next(i for i, line in enumerate(lines) if line.strip().startswith("("))
    lines = lines[data_start:]
    return _load_potential_lines_as_array(lines)


def _load_potential_lines_as_array(lines: list[str]) -> np.ndarray:
    values = []
    for line in lines:
        real_str, imag_str = line.strip()[1:-1].split(",")
        values.append(complex(float(real_str), float(imag_str)))
    return np.asarray(values, dtype=np.complex128)



def _raw_potential_in_input_file_convention(
    potential: ScatteringOperator,
) -> np.ndarray[Any, np.dtype[np.complex128]]:

    potential_diagonal = array.extract_diagonal(potential)
    nx, ny, nz = potential_diagonal.basis.metadata().shape
    basis_weights = potential_diagonal.basis.metadata().children[2].basis_weights
    basis = close_coupling_basis(potential_diagonal.basis.metadata())

    data = potential_diagonal.with_basis(basis).raw_data
    data = data.reshape((nx, ny, nz)) * (basis_weights[np.newaxis, np.newaxis, :])

    # Multiscat uses a slightly different fourier convention
    return data.ravel() / (electron_volt * 10**-3 * np.sqrt(nx * ny))

def _potential_to_input_file(potential: ScatteringOperator) -> str:
    potential_lobatto = _raw_potential_in_input_file_convention(potential)
    metadata = potential.basis.metadata().children[0]
    nx, ny, nz = metadata.shape
    z_domain = metadata.children[2].domain
    z_start_angstrom = z_domain.start / angstrom
    z_end_angstrom = (z_domain.start + z_domain.delta) / angstrom
    nfc = nx * ny
    metadata_x01, _ = split_scattering_metadata(metadata)
    directions = metadata.extra.vectors
    x_vector = np.asarray(directions[0]) * metadata_x01.children[0].domain.delta
    y_vector = np.asarray(directions[1]) * metadata_x01.children[1].domain.delta
    ax_angstrom = x_vector[0] / angstrom
    ay_angstrom = x_vector[1] / angstrom
    bx_angstrom = y_vector[0] / angstrom
    by_angstrom = y_vector[1] / angstrom

    header_lines = [
        "Generated from ScatteringCondition.potential",
        "Metadata line follows: nfc nkx nky nzlobatto",
        f"{nfc} {nx} {ny} {nz}",
        "Unit cell vectors in Angstrom: ax ay bx by",
        f"{ax_angstrom:.10g} {ay_angstrom:.10g} {bx_angstrom:.10g} {by_angstrom:.10g}",
        "Integration range in Angstrom: zmin zmax",
        f"{z_start_angstrom:.10g} {z_end_angstrom:.10g}",
        "Format: (real, imag)",
        "Ordering: Fourier component then z-slice",
        "Do not edit by hand",
    ]
    data_lines = [
        f"({value.real:+.6e}, {value.imag:+.6e})"
        for value in potential_lobatto.reshape(-1)
    ]
    return "\n".join([*header_lines, *data_lines]) + "\n"


def _simple_example_condition() -> tuple[ScatteringCondition, OptimizationConfig]:

    HELIUM_MASS = physical_constants["alpha particle mass"][0]
    UNIT_CELL = 2.84 * angstrom
    Z_HEIGHT = 8 * angstrom

    MORSE_PARAMETERS = operator.build.CorrugatedMorseParameters(
        depth=7.63 * electron_volt * 10**-3,
        height=(1.0 / 1.1) * angstrom,
        offset=3.0 * angstrom,
        beta=0.10,
    )

    metadata = scattering_metadata_from_stacked_delta_x(
        (
            np.array([UNIT_CELL, 0, 0]),
            np.array([0, UNIT_CELL, 0]),
            np.array([0, 0, Z_HEIGHT]),
        ),
        (32, 32, 550),
    )
    # This is taken from https://doi.org/10.1039/FT9908601641
    # and is a reproduction of the Wolken 4He-LiF problem in table 1,
    # originally simulated in https://doi.org/10.1063/1.1679617.
    condition = ScatteringCondition.from_angles(
        mass=HELIUM_MASS,
        energy=20 * electron_volt * 10**-3,
        theta=np.deg2rad(30),
        phi=np.deg2rad(0),
        potential=operator.build.corrugated_morse_potential(
            metadata,
            MORSE_PARAMETERS,
        ),
    )
    config = OptimizationConfig(precision=1e-5, max_iterations=1000)
    return condition, config


def _run_multiscat_cli(
    condition: ScatteringCondition, config: OptimizationConfig
) -> dict[tuple[int, int], float]:

    binary = ROOT / "multiscat"
    if not binary.exists():
        subprocess.run(["make", "multiscat"], cwd=ROOT, check=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_path = Path(temp_dir)
        (tmp_path / "potential.in").write_text(
            _potential_to_input_file(condition.potential)
        )
        (tmp_path / "condition.in").write_text(_condition_to_input_file(condition))
        (tmp_path / "optimization.conf").write_text(_optimization_to_input_file(config))

        subprocess.run(
            [
                str(binary),
                "--potential",
                "potential.in",
                "--condition",
                "condition.in",
                "--optimization",
                "optimization.conf",
            ],
            cwd=tmp_path,
            check=True,
        )
        output_file = tmp_path / "diffrac.out"
        assert output_file.exists(), "Expected diffrac.out to be generated"
        return _parse_intensities(output_file)


def test_simple_system() -> None:

    condition, config = _simple_example_condition()
    intensities = _run_multiscat_cli(condition, config)

    assert math.isclose(sum(intensities.values()), 1.0, abs_tol=1e-6)

    expected_from_file = _parse_raw_intensities(
        TESTS_DIR / Path("expected_intensities.txt")
    )
    for spot, expected_value in expected_from_file.items():
        assert spot in intensities, f"Missing diffraction spot {spot}"
        assert math.isclose(intensities[spot], expected_value, abs_tol=1e-5)


def _rotated_example_condition() -> tuple[ScatteringCondition, OptimizationConfig]:

    HELIUM_MASS = physical_constants["alpha particle mass"][0]
    UNIT_CELL = 2.84 * angstrom
    Z_HEIGHT = 8 * angstrom

    MORSE_PARAMETERS = operator.build.CorrugatedMorseParameters(
        depth=7.63 * electron_volt * 10**-3,
        height=(1.0 / 1.1) * angstrom,
        offset=3.0 * angstrom,
        beta=0.10,
    )
    rotation = np.deg2rad(20.0)
    cos_t = np.cos(rotation)
    sin_t = np.sin(rotation)
    rotation_matrix = np.array(
        [
            [cos_t, -sin_t, 0.0],
            [sin_t, cos_t, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    x_vector = rotation_matrix @ np.array([UNIT_CELL, 0.0, 0.0])
    y_vector = rotation_matrix @ np.array([0.0, UNIT_CELL, 0.0])

    metadata = scattering_metadata_from_stacked_delta_x(
        (
            x_vector,
            y_vector,
            np.array([0.0, 0.0, Z_HEIGHT]),
        ),
        (10, 10, 550),
    )

    condition = ScatteringCondition.from_angles(
        mass=HELIUM_MASS,
        energy=20 * electron_volt * 10**-3,
        theta=np.deg2rad(30),
        phi=np.deg2rad(0),
        potential=operator.build.corrugated_morse_potential(
            metadata,
            MORSE_PARAMETERS,
        ),
    )
    config = OptimizationConfig(precision=1e-5, max_iterations=1000)
    return condition, config


def test_rotated_system() -> None:

    condition, config = _rotated_example_condition()
    intensities = _run_multiscat_cli(condition, config)

    assert intensities, "Expected at least one diffraction intensity"
    assert math.isclose(sum(intensities.values()), 1.0, abs_tol=1e-6)


def test_raw_potential_in_input_file_convention() -> None:
    condition, _ = _simple_example_condition()
    from_condition = _raw_potential_in_input_file_convention(condition.potential)

    reference_potential = TESTS_DIR / Path("pot10001.in")
    expected = _load_potential_file_as_array(reference_potential)

    assert expected.shape == from_condition.shape
    np.testing.assert_allclose((from_condition), (expected), rtol=1e-5, atol=1e-10)
