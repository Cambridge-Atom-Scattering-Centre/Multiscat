from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path
import tempfile

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
    scattering_metadata_from_stacked_delta_x,
    split_scattering_metadata,
)
from multiscat.config import OptimizationConfig, ScatteringCondition


ROOT = Path(__file__).resolve().parents[1]
TESTS_DIR = Path(__file__).resolve().parent


def _parse_raw_intensities(output_file: Path) -> dict[tuple[int, int], float]:
    # Regex for lines without the '#' prefix: two ints and one float
    pattern = re.compile(r"^\s*(-?\d+)\s+(-?\d+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*")
    intensities: dict[tuple[int, int], float] = {}
    print(output_file.read_text())
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
    print(output_file.read_text())
    with output_file.open("r") as f:
        for line in f:
            match = pattern.match(line)
            if match:
                h = int(match.group(1))
                k = int(match.group(2))
                val = float(match.group(3))
                intensities[(h, k)] = val
    return intensities


def _scat_cond_from_condition(condition: ScatteringCondition) -> str:
    scattering_vector = np.asarray(condition.incident_k)
    scattering_magnitude = float(np.linalg.norm(scattering_vector))
    assert scattering_magnitude > 0, "Incident wavevector magnitude must be non-zero"

    energy_meV = condition.incident_energy / (electron_volt * 10**-3)
    theta_degrees = np.degrees(condition.theta)
    phi_degrees = np.degrees(condition.phi)

    scat_cond_lines = [
        "Comment line: Energy, theta, phi   This file defines the combination of conditions to be used by multiscat",
        f"{energy_meV:.10g},{theta_degrees:.10g},{phi_degrees:.10g}",
    ]
    return "\n".join(scat_cond_lines) + "\n"


def _multiscat_conf_from_condition(
    condition: ScatteringCondition,
    config: OptimizationConfig,
) -> str:
    mass_amu = condition.mass / atomic_mass
    _a = condition.metadata.children[2].domain
    z_start_angstrom = _a.start / angstrom
    z_end_angstrom = (_a.start + _a.delta) / angstrom
    metadata_x01, _ = split_scattering_metadata(condition.metadata)
    directions = condition.metadata.extra.vectors
    x_vector = np.asarray(directions[0]) * metadata_x01.children[0].domain.delta
    y_vector = np.asarray(directions[1]) * metadata_x01.children[1].domain.delta
    a1_angstrom = x_vector[0] / angstrom
    a2_angstrom = y_vector[0] / angstrom
    b2_angstrom = y_vector[1] / angstrom

    lines = [
        "scatCond.in\t! The scattering conditions input file",
        "1       !itest=1 enables output of each diffraction intensity; itest=0 outputs specular only",
        "0       !gmres preconditioner flag (ipc)",
        f"{int(np.log10(1 / config.precision))}       !number of significant figures convergence (nsf)",
        f"{z_start_angstrom:.10g},{z_end_angstrom:.10g}       !integration range (zmin,zmax)",
        "120       !max -ve energy of closed channels (dmax)",
        "120       !max index of channels (imax)",
        f'{a1_angstrom:.10g}       !a1 (see subroutine basis in "scatsub.f")',
        f"{a2_angstrom:.10g}       !a2",
        f"{b2_angstrom:.10g}       !b2",
        "10001       !startindex",
        "10001       !endindex",
        f"{mass_amu:.10g}       !helium mass",
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


def _potential_from_condition(condition: ScatteringCondition) -> str:
    potential_lobatto = _raw_potential_in_input_file_convention(condition)
    nx, ny, nz = condition.metadata.shape
    nfc = nx * ny

    header_lines = [
        "Generated from ScatteringCondition.potential",
        "Metadata line follows: nfc nkx nky nzlobatto",
        f"{nfc} {nx} {ny} {nz}",
        "Format: (real, imag)",
        "Ordering: Fourier component then z-slice",
        "Do not edit by hand",
    ]
    data_lines = [
        f"({value.real:+.6e}, {value.imag:+.6e})"
        for value in potential_lobatto.reshape(-1)
    ]
    return "\n".join([*header_lines, *data_lines]) + "\n"


def _raw_potential_in_input_file_convention(
    condition: ScatteringCondition,
) -> np.ndarray:

    nx, ny, nz = condition.metadata.shape

    # Potential data from the Lobatto-basis operator is weighted in z.
    # Convert to unweighted potential values in meV while preserving
    # Lobatto z nodes.
    potential_lobatto = array.extract_diagonal(condition.potential).raw_data.reshape(
        (nx, ny, nz)
    )
    potential_lobatto = (
        potential_lobatto
        * (condition.metadata.children[2].basis_weights[np.newaxis, np.newaxis, :])
        / (electron_volt * 10**-3)
    )

    # Convert real-space potential samples to Fourier components at each
    # Lobatto z node.
    return np.fft.fft2(potential_lobatto, axes=(0, 1)).reshape(-1) / (nx * ny)


def _manual_example_condition() -> tuple[ScatteringCondition, OptimizationConfig]:

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
        phi=0,
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
        (tmp_path / "pot10001.in").write_text(_potential_from_condition(condition))
        (tmp_path / "scatCond.in").write_text(_scat_cond_from_condition(condition))
        (tmp_path / "Multiscat.conf").write_text(
            _multiscat_conf_from_condition(condition, config)
        )

        subprocess.run([str(binary), "Multiscat.conf"], cwd=tmp_path, check=True)
        output_file = tmp_path / "diffrac10001.out"
        assert output_file.exists(), "Expected diffrac10001.out to be generated"
        return _parse_intensities(output_file)


def test_manual_lif_exercise_intensities() -> None:

    condition, config = _manual_example_condition()
    intensities = _run_multiscat_cli(condition, config)

    assert math.isclose(sum(intensities.values()), 1.0, abs_tol=1e-6)

    expected_from_file = _parse_raw_intensities(
        TESTS_DIR / Path("expected_intensities.txt")
    )
    for spot, expected_value in expected_from_file.items():
        assert spot in intensities, f"Missing diffraction spot {spot}"
        assert math.isclose(intensities[spot], expected_value, abs_tol=1e-5)


def test_raw_potential_in_input_file_convention() -> None:
    condition, _ = _manual_example_condition()
    from_condition = _raw_potential_in_input_file_convention(condition)

    reference_potential = TESTS_DIR / Path("pot10001.in")
    reference_potential.write_text(_potential_from_condition(condition))
    expected = _load_potential_file_as_array(reference_potential)

    assert expected.shape == from_condition.shape
    np.testing.assert_allclose((from_condition), (expected), rtol=1e-5)
