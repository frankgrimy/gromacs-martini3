#!/bin/bash
# This script runs a simulation of a protein using the unmodified Martini 3 force field.

# Output colors
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[0;33m'
CLEAR='\033[0m'
STRNAME="fold_foxp1_solo_qrich_2025_01_07_11_21_model_0"

# Repo root
REPOROOT=$(git rev-parse --show-toplevel)

# Unmodified Martini 3 workdir
BASEDIR=$REPOROOT/comparison/workdir/base

# Check if the workdir exists, if not, create it
if [ ! -d $REPOROOT/comparison/workdir ]; then
  /usr/bin/mkdir $REPOROOT/comparison/workdir
fi

# Create the workdir for the Martini 3 unmodified force field run
/usr/bin/mkdir $REPOROOT/comparison/workdir/base

# Unzip the protein structure
printf "\n${GREEN}Unzipping the protein structure...${CLEAR}\n"
/usr/bin/unzip $REPOROOT/comparison/src/fold_foxp1_solo_qrich_2025_01_07_11_21.zip "${STRNAME}.cif" -d $REPOROOT/comparison/workdir/base
#/usr/bin/unzip src/fold_foxp1_solo_qrich_2025_01_07_11_21.zip "${STRNAME}" -d workdir/base
printf "${GREEN}Unzipping the protein structure...OK${CLEAR}\n\n"

# Convert the CIF file to PDB
printf "\n${GREEN}Converting the CIF file to PDB...${CLEAR}\n"
cd $BASEDIR
/usr/bin/obabel ./$STRNAME.cif -O $STRNAME.pdb
printf "${GREEN}Converting the CIF file to PDB...OK${CLEAR}\n\n"

# Clean the structure
printf "\n${GREEN}Cleaning the structure...${CLEAR}\n"
/usr/bin/grep -v HETATM $STRNAME.pdb > $STRNAME-clean.pdb
printf "${GREEN}Cleaning the structure...OK${CLEAR}\n\n"

# Convert to Coarse-Grain
printf "${GREEN}Converting to Coarse-Grain witn martinize2...${CLEAR}\n"
martinize2 -f $STRNAME-clean.pdb \
  -x $STRNAME-cg.pdb \
  -o topol.top \
  -ff martini3001 \
  -cys auto \
  -p backbone \
  -elastic \
  -ef 700.0 \
  -el 0.5 \
  -eu 0.9 \
  -dssp mkdssp
printf "${GREEN}Converting to Coarse-Grain witn martinize2...OK${CLEAR}\n\n"

# Prepare the water box
printf "\n${GREEN}Preparing the water box...${CLEAR}\n"
gmx editconf -f $STRNAME-cg.pdb \
  -o $STRNAME-cg-box.gro \
  -bt dodecahedron \
  -d 4.0
printf "${GREEN}Preparing the water box...OK${CLEAR}\n\n"

# Solvate with insane (Python library)
printf "${GREEN}Solvating the system with insane...${CLEAR}\n"
insane -f $STRNAME-cg-box.gro \
  -o system.gro \
  -p topol.top \
  -sol W \
  -salt 0.15 \
  -center \
  -pbc keep
printf "${GREEN}Solvating the system with insane...OK${CLEAR}\n\n"

# Ask if the user wants to check the topology file
printf "${YELLOW}You should edit the topology file to add the CG beads\ndefinitions from Martini 3 force field.\nReplace the Martini include with the following lines, in the topology file:\n\n"
printf "${CYAN}#include \"martini_v3.0.0.itp\"\n"
printf "#include \"martini_v3.0.0_solvents_v1.itp\"\n"
printf "#include \"martini_v3.0.0_ions_v1.itp\"\n"
printf "#include \"molecule_0.itp\"\n${CLEAR}"

printf "\n${YELLOW}Also change the protein name to ${CYAN}molecule_0${YELLOW}, and remove the charge signs from the ions.${CLEAR}\n"

read -n 1 -p "Do you want to edit the topology file? [y/N] " checktopology
echo

if [[ "$checktopology" =~ ^[Yy]$ ]]; then
  nano topol.top
fi

printf "\n${YELLOW}You have to remove the charge signs from the ions in the .gro file (second column).\nThe format of the ions lines should be:\n\n"
printf "${CYAN}  51xxxNA+    ${YELLOW}NA 51xxx${CYAN}  16.548   9.165  11.022\n${CLEAR}"
read -p "Press Enter to edit the .gro file with nano..."
nano system.gro
printf "\n"

# Energy minimization
## Ask the user if they have the minim.mdp file, if not, suggest to configure it and stop the run
if [ ! -f minim.mdp ]; then
  printf "${RED}The minim.mdp file is missing. Please create it and configure it for the energy minimization.${CLEAR}\n"
  exit 1
fi

## Ask the user if they want to run the energy minimization
read -n 1 -p "Do you want to run the energy minimization? [y/N] " runminim

if [[ "$runminim" =~ ^[Yy]$ ]]; then
  gmx grompp -f minim.mdp \
  -c system.gro \
  -r system.gro \
  -p topol.top \
  -o minim.tpr

  #gmx mdrun -v -s minim.tpr
  gmx mdrun -v -deffnm minim

  read -n 1 -p "Do you want to generate the energy minimization results? [y/N] " checkminim

  if [[ "$checkminim" =~ ^[Yy]$ ]]; then
    printf "\n${YELLOW}When prompted, type ${CYAN}8 0${YELLOW} to select the Potential Energy.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx energy -f minim.edr -o potential.xvg
    printf "${GREEN}The potential energy has been saved to ${CYAN}potential.xvg\n${GREEN}You can plot it with ${CYAN}xmgrace potential.xvg${CLEAR}\n\n"
  fi
fi

# Equilibration
## NPT equilibration
if [ ! -f npt.mdp ]; then
  printf "${RED}The npt.mdp file is missing. Please create it and configure it for the NPT equilibration.${CLEAR}\n"
  exit 1
fi

read -n 1 -p "Do you want to run the NPT equilibration? [y/N] " runnpt

if [[ "$runnpt" =~ ^[Yy]$ ]]; then
  gmx grompp -f npt.mdp \
  -c minim.gro \
  -r minim.gro \
  -p topol.top \
  -o npt.tpr \
  -maxwarn 2

  gmx mdrun -v -deffnm npt

  read -n 1 -p "Do you want to generate the NPT equilibration results? [y/N] " checknpt

  if [[ "$checknpt" =~ ^[Yy]$ ]]; then
    printf "\n${YELLOW}When prompted, type ${CYAN}13 0${YELLOW} to select the Pressure.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx energy -f npt.edr -o pressure.xvg
    printf "${GREEN}The pressure data has been saved to ${CYAN}pressure.xvg\n${GREEN}You can plot it with ${CYAN}xmgrace pressure.xvg${CLEAR}\n"

    printf "\n${YELLOW}When prompted, type ${CYAN}19 0${YELLOW} to select the Density.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx energy -f npt.edr -o density.xvg
    printf "${GREEN}The density data has been saved to ${CYAN}density.xvg\n${GREEN}You can plot it with ${CYAN}xmgrace density.xvg${CLEAR}\n"
  fi
fi

# Production run
if [ ! -f prod.mdp ]; then
  #echo "The md.mdp file is missing. Please create it and configure it for the production run."
  printf "${RED}The prod.mdp file is missing. Please create it and configure it for the production run.${CLEAR}\n"
  exit 1
fi

read -n 1 -p "Do you want to run the production run? This may take a LONG time: [y/N] " runmd

if [[ "$runmd" =~ ^[Yy]$ ]]; then
  gmx grompp -f prod.mdp \
  -c npt.gro \
  -t npt.cpt \
  -p topol.top \
  -o prod.tpr

  gmx mdrun -v -deffnm prod

  read -n 1 -p "Do you want to generate the production run results? [y/N] " checkmd

  if [[ "$checkmd" =~ ^[Yy]$ ]]; then
    printf "\n${GREEN}Correcting the trajectory...${CLEAR}\n"
    printf "${YELLOW}When prompted, type ${CYAN}1 (Protein), then 0 (System)${CLEAR}\n"
    read -p "Press Enter to continue..."

    gmx trjconv -pbc mol -center -ur compact -s prod.tpr -f prod.xtc -o temp.xtc

    printf "\n${YELLOW}When prompted, type ${CYAN}a BB${YELLOW}, then ${CYAN}q${YELLOW} to align the protein to the backbone.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx make_ndx -f system.gro -o index.ndx
    gmx trjconv -fit rot+trans -s prod.tpr -f temp.xtc -o trajfitted.xtc -n index.ndx

    printf "\n${GREEN}Calculating the RMSD...${CLEAR}\n"
    printf "${YELLOW}When prompted, select ${CYAN}BB${YELLOW} for the group, and ${CYAN}0${YELLOW} for System.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx rms -f trajfitted.xtc -s prod.tpr -n index.ndx -o rmsd.xvg
    printf "${GREEN}The RMSD data has been saved to ${CYAN}rmsd.xvg\n${GREEN}You can plot it with ${CYAN}xmgrace rmsd.xvg${CLEAR}\n"

    printf "\n${GREEN}Calculating the gyration radius...${CLEAR}\n"
    printf "${YELLOW}When prompted, select ${CYAN}BB${YELLOW} for the group, and ${CYAN}0${YELLOW} for System.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx gyrate \
    -s prod.tpr \
    -f trajfitted.xtc \
    -o gyrate.xvg \
    -n index.ndx
    printf "${GREEN}The gyration radius data has been saved to ${CYAN}gyrate.xvg\n${GREEN}You can plot it with ${CYAN}xmgrace gyrate.xvg${CLEAR}\n"

    printf "\n${GREEN}  Calculating the average gyration radius...${CLEAR}\n"
    awk '{print $1, $2}' gyrate.xvg > rg_clean.xvg
    gmx analyze -f rg_clean.xvg

    printf "\n${GREEN}Calculating RMSF...${CLEAR}\n"
    printf "${YELLOW}When prompted, select ${CYAN}BB${YELLOW} for the group, and ${CYAN}0${YELLOW} for System.${CLEAR}\n"
    read -p "Press Enter to continue..."
    gmx rmsf \
    -s prod.tpr \
    -f trajfitted.xtc \
    -o rmsf.xvg \
    -n index.ndx
    printf "${GREEN}The RMSF data has been saved to ${CYAN}rmsf.xvg\n${GREEN}"
    printf "Prior to plotting, you need to modify the second column of the file, since the numbering must be sequential.${CLEAR}\n"
    read -p "Do you want to edit the file now? [y/N] " editrmsf

    if [[ "$editrmsf" =~ ^[Yy]$ ]]; then
      nano rmsf.xvg
    fi
  fi
fi
