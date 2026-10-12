# Setup: Python, uv and git

Python runs the ICM scripts in `scripts/` (standard library only) and the N-body adapters. git and the GitHub CLI handle commits and pull requests.

## Install

1. Python 3.12 from https://www.python.org/downloads/ (installed on the main machine).
2. uv from https://docs.astral.sh/uv/ for virtual environments: on the main machine the system Python lacks the `venv` module, so use `uv venv` instead of `python -m venv`.
3. Git for Windows from https://git-scm.com/ and the GitHub CLI from https://cli.github.com/ ; run `gh auth login` once.

## Verify

- `python --version` prints 3.12.x; `uv --version` prints a version.
- `git --version` and `gh auth status` succeed.

## How this repo uses it

- `python scripts/<name>.py` for ICM mechanics (no third-party packages).
- `uv venv` plus `uv pip install` for the N-body environments (`rebound.md`).
- Every push, pull request or submission needs explicit human approval.
