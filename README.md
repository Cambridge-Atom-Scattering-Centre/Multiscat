# Multiscat
Please refer to the manual attached as a pdf file called multiscat_manual.

## Testing

Build the `multiscat` executable first:

```bash
make
```

Then run the test suite with pytest from the project virtual environment:

```bash
./.venv/bin/python -m pytest -q
```

Run only the manual example regression test:

```bash
./.venv/bin/python -m pytest -q tests/test_manual_example.py
```

Clean all gitignored files (build artifacts, caches, outputs, etc.):

```bash
git clean -fdX
```
