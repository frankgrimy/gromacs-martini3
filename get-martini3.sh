#!/bin/bash

MARTINI_VER="v300"
MARTINI_VERDOTS="v3.0.0"

# Set up the file structure
GMXLIB=$(pwd)
FFDIR=$GMXLIB/share/gromacs/top
/usr/bin/mkdir -p $FFDIR

# Tell GROMACS where the custom force fields are
export GMXLIB=$GMXLIB

# Download Martini 3 force field
## Particle definitions
/usr/bin/wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_${MARTINI_VER}.zip
/usr/bin/unzip martini_${MARTINI_VER}.zip -d martini_${MARTINI_VER}

/usr/bin/cp -r martini_${MARTINI_VER}/martini_${MARTINI_VER}/* $FFDIR
/usr/bin/rm -rf martini_${MARTINI_VER} martini_${MARTINI_VER}.zip # Clean up

# Configure Python venv
/usr/bin/python3 -m venv .venv # Create a virtual environment
source .venv/bin/activate # Activate the virtual environment
pip install -r requirements.txt # Install the required Python packages