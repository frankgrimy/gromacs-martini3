# Scripts to set up a Martini 3 environment for protein simulations

## Requirements

- A working installation of GROMACS on Linux or WSL
- Python 3.9 or later
- wget (get it from your package manager)
- unzip (get it from your package manager)

## Configuration

- Clone this repo and change directory to it.
- Run `./get-martini3.sh` to configure the Martini 3 environment. This will:
  - Download the [Martini 3 force field](https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/particle-definitions.html) and add a custom path for it.
  - Create a Python virtual environment in `.venv` and install [vermouth/martinize2](https://github.com/marrink-lab/vermouth-martinize?tab=readme-ov-file#installation) and [mdtraj](https://pypi.org/project/mdtraj/) to provide DSSP functionality.

## Recommended usage

Create a folder named _workdir_ and place all the needed files there. Change directory to _workdir_ and then you should be able to run `martinize2` and `gmx` commands from there.

## Notes

- If you have exited the shell after configuring the environment and want to resume the work, you can reactivate the Python virtual environment by running `source .venv/bin/activate` on the repository root, instead of running the configuration script again.