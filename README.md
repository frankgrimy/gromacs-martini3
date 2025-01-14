# Scripts to set up a Martini 3 environment for protein simulations

## Requirements

- A working installation of GROMACS on Linux or WSL
- Python 3.9 or later
- wget (get it from your package manager)
- unzip (get it from your package manager)

## Configuration

- Clone this repo and change directory to it.
- Run `./get-martini3.sh` to configure the Martini 3 environment. This will:
  - Download the [Martini 3 force field](https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/particle-definitions.html) and add a custom path for it, inside the repo folder.
  - Create a Python virtual environment in `.venv` and install [vermouth/martinize2](https://github.com/marrink-lab/vermouth-martinize?tab=readme-ov-file#installation) and [mdtraj](https://pypi.org/project/mdtraj/) to provide DSSP functionality.

## Recommended usage

Create a folder named _workdir_ and place all the needed files there. Change directory to _workdir_ and then you should be able to run `martinize2`, `insane` and `gmx` commands from there.

## Notes

- If you exited the shell after configuring the environment and want to resume your work, you can reactivate the Python virtual environment by running `source .venv/bin/activate` inside the repository root, instead of re-running the configuration script.
- The Martini 3 force field will be installed in the `share/gromacs/top` directory of the repository that is created by the script, and the `GMXLIB` environment variable will be set to this path. If you want to use other force fields remember to set the `GMXLIB` variable accordingly, or unset it to revert to the default GROMACS force fields. 
