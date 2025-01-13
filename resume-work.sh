#!/bin/bash

GMXLIB=$(pwd)/share/gromacs/top
export GMXLIB=$GMXLIB

source .venv/bin/activate
echo "Virtual environment activated"
$SHELL