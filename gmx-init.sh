#!/bin/bash

MARTINI_VER="v300"
MARTINI_VERDOTS="v3.0.0"

# Set up the file structure
GMXLIB=$(pwd)
FFDIR=$GMXLIB/share/gromacs/top
mkdir -p $FFDIR

# Tell GROMACS where the custom force fields are
export GMXLIB=$GMXLIB

# Download Martini 3 force field
## Particle definitions
wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_${MARTINI_VER}.zip
unzip martini_${MARTINI_VER}.zip -d martini_${MARTINI_VER}

cp -r martini_${MARTINI_VER}/martini_${MARTINI_VER}/* $FFDIR
rm -rf martini_${MARTINI_VER} martini_${MARTINI_VER}.zip # Clean up