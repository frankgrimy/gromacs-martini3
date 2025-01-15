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

## [Martini tutorial: Proteins - Part 1](https://cgmartini.nl/docs/tutorials/Legacy/martini3/ProteinsI/)

The script located in `martini-basics/martini3-tutorial1.sh` runs the simulation example from the [Martini 3 tutorial](https://cgmartini.nl/docs/tutorials/Legacy/martini3/ProteinsI/).  
You can use it as a base to verify the installation and even run your own simulations.

## Notes

- If you exited the shell after configuring the environment and want to resume your work, you can reactivate the Python virtual environment by running `source .venv/bin/activate` inside the repository root, instead of re-running the configuration script.