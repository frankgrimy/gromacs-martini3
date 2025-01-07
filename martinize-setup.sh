#!/bin/bash

# Configure Python venv
python3 -m venv .venv # Create a virtual environment
source .venv/bin/activate # Activate the virtual environment
pip install -r requirements.txt # Install the required Python packages