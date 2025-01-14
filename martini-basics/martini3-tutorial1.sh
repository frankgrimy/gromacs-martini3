#!/bin/bash

# Output colors
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[0;33m'
CLEAR='\033[0m'

printf "This script is based on the official Martini 3 tutorial (available at ${CYAN}https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/${CLEAR})\n"
printf "${YELLOW}Remember to activate the virtual environment before running this script.${CLEAR}\n\n"

# Download tutorial files
printf "${GREEN}Downloading tutorial files...${CLEAR}\n"
#/usr/bin/wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/0_Tutorials/m3_tutorials/ProteinsI/M3_proteins_tutorial_part1.zip
cp /tmp/M3_proteins_tutorial_part1.zip ./
/usr/bin/unzip M3_proteins_tutorial_part1.zip "tutorial_2/*" -d ./
/usr/bin/rm -rf M3_proteins_tutorial_part1.zip # Clean up
printf "${GREEN}Downloading tutorial files...OK${CLEAR}\n\n"

# Change to the tutorial directory
/usr/bin/mkdir tutorial_2/workdir
cd tutorial_2/workdir
/usr/bin/cp ../template/181L.pdb . # Copy the PDB file to the workdir
/usr/bin/cp ../template/martini_v3.0.0.itp .
/usr/bin/cp ../template/*.mdp .

# Clean the structure
/usr/bin/grep "^ATOM" 181L.pdb > 181L_clean.pdb

# Convert to Coarse-Grain (martinize2)
printf "${GREEN}Converting to Coarse-Grain witn martinize2...${CLEAR}\n"
martinize2 -f 181L_clean.pdb \
  -o t4l_only.top \
  -x t4l_cg.pdb \
  -p backbone \
  -ff martini3001 \
  -dssp
printf "${GREEN}Converting to Coarse-Grain witn martinize2...OK${CLEAR}\n\n"

# Minimization in vacuum
printf "${GREEN}Preparing minimization in vacuum..."
printf "  Creating the box...${CLEAR}\n"
gmx editconf -f t4l_cg.pdb \
  -d 1.0 \
  -bt dodecahedron \
  -o t4l_cg.gro
printf "${GREEN}  Creating the box...OK${CLEAR}\n\n"

# Ask the user to edit the topology file
printf "${YELLOW}You have to edit the topology file to fix the .itp file path.${CLEAR}\nReplace the Martini include directive with the following line:\n${CYAN}#include \"../template/martini_v3.0.0.itp\"${CLEAR}\n"
read -p "Press Enter to edit the topology file with nano..."
nano t4l_only.top

printf "  Generating minimization .tpr file...${CLEAR}\n"
gmx grompp -p t4l_only.top \
  -f minimization.mdp \
  -c t4l_cg.gro \
  -o minimization-vac.tpr \
  -r t4l_cg.gro
printf "${GREEN}  Generating minimization .tpr file...OK${CLEAR}\n"
printf "${GREEN}Preparing minimization in vacuum...OK${CLEAR}\n\n"

printf "Running minimization in vacuum...${CLEAR}\n"
gmx mdrun -v -deffnm minimization-vac
printf "${GREEN}Running minimization in vacuum...OK${CLEAR}\n\n"

# # Solvating the system
# printf "${GREEN}Solvating the system with gmx solvate...${CLEAR}\n"
# /usr/bin/cp ../template/water.gro .
# gmx solvate -cp minimization-vac.gro \
#   -cs water.gro \
#   -radius 0.21 \
#   -o solvated.gro
# printf "${GREEN}Solvating the system with gmx solvate...OK${CLEAR}\n\n"

# # Adding ions with insane
# printf "${GREEN}Adding ions with insane...${CLEAR}\n"
# insane -f solvated.gro \
#   -o solvated_ions.gro \
#   -p t4l_only.top \
#   -pbc dodecahedron \

# Solvating and adding ions with insane
printf "${GREEN}Solvating and adding ions with insane...${CLEAR}\n"
insane -f t4l_cg.gro \
  -o t4l_cg_sol_ions.gro \
  -p t4l_only.top \
  -pbc keep \
  -d 1.0 \
  -sol W \
  -salt 0.15 \
  -box 10,10,10 \
  -center
printf "${GREEN}Solvating and adding ions with insane...OK${CLEAR}\n\n"

# Ask the user to edit the topology file
printf "${YELLOW}You have to edit the topology file to fix the .itp file path.${CLEAR}\nReplace the Martini include directive with the following lines:\n\
${CYAN}#include \"../template/martini_v3.0.0.itp\"\n\
#include \"../template/martini_v3.0.0_solvents_v1.itp\"\n\
#include \"../template/martini_v3.0.0_ions_v1.itp\"\n\
#include \"molecule_0.itp\"${CLEAR}\n"

printf "${YELLOW}You also have to remove the charge signs from the ions, and also change the name of the protein molecule to 'molecule_0' in the [molecules] section.${CLEAR}\n"
read -p "Press Enter to edit the topology file with nano..."
nano t4l_only.top

# Ask the user to remove the charge signs from the ions in the .gro file
printf "${YELLOW}You have to remove the charge signs from the ions in the .gro file.${CLEAR}\n"
read -p "Press Enter to edit the .gro file with nano..."
nano t4l_cg_sol_ions.gro

# Minimization in water
printf "${GREEN}Preparing minimization in water...${CLEAR}\n"
printf "  Generating minimization .tpr file...${CLEAR}\n"
gmx grompp -p t4l_only.top \
  -c t4l_cg_sol_ions.gro \
  -f minimization.mdp \
  -o minimization.tpr \
  -r t4l_cg_sol_ions.gro
printf "${GREEN}  Generating minimization .tpr file...OK${CLEAR}\n"

printf "Running minimization in water...${CLEAR}\n"
gmx mdrun -v -deffnm minimization
printf "${GREEN}Running minimization in water...OK${CLEAR}\n\n"

