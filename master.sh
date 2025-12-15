#!/bin/bash

# ==============================================================================
# Script Name: master.sh (MASTER NODE / SERIAL VERSION)
# Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi
# Email: dheeraj@scfbio-iitd.res.in
# Date: December 14, 2025
# Version: 1.3.0 (Strictly for DRUGBANK/FDA on Master Node)
#
# Description:
#   This script performs virtual screening on a provided PDB file against 
#   specific small molecule databases (DrugBank or FDA).
#   
#   USAGE GUIDELINES:
#   1. This script is designed for SERIAL execution (no parallel background jobs).
#   2. Run this directly on the master node or an interactive session.
#   3. Arguments:
#      $1 : PDB Filename (without extension)
#      $2 : Ligand Code (e.g., LIG)
#      $3 : Database Identifier ("DB" for DrugBank, "FDA" for FDA)
#
#   Example: 
#      ./master.sh 1abc.pdb LIG DB
#
# Copyright (c) 2025 Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio), IIT Delhi.
# All rights reserved.
# ==============================================================================

date
# --- Assign variables ---

JOB_DIR=`pwd`
JOBID=`date +%s`
PDB_FILE=$1
LIGAND_CODE=$2
DB=$3

# ==============================================================================
# SECTION 1: JOB SETUP AND STATUS REPORTING
# ==============================================================================

# Exit immediately if a command fails.
set -e

if [ -z "$SEARCH_ML_HOME" ]; then
    echo "Error: SEARCH_ML_HOME is not defined."
    echo "Please set the SEARCH-ML home path by executing the below commands:"
    echo "   conda env config vars set SEARCH_ML_HOME="$(pwd)" --name $ENV_NAME"
    echo "   conda deactivate $ENV_NAME"
    echo "   conda activate $ENV_NAME"
    echo ""
    exit 1
else
    echo "SEARCH_ML_HOME is set to: $SEARCH_ML_HOME"
fi

# Create the job directory
mkdir $JOBID

# CRITICAL: Navigate to the specific job directory.
cd "$JOB_DIR/$JOBID" || { echo "FATAL: Could not navigate to job directory $JOB_DIR/$JOBID"; exit 1; }

# --- Define a function to handle errors ---
handle_error() {
  local exit_code=$?
  echo "Job failed with exit code $exit_code at or near line $1." > error.log
  touch FAILED
  rm -f RUNNING
  exit $exit_code
}

# --- Set the trap ---
trap 'handle_error $LINENO' ERR

# --- Signal that the job has officially started ---
touch RUNNING
echo "$1 , $2 , $3" >>args.txt
exec > >(tee -a job.log)
exec 2> >(tee -a job.log >&2)

echo "Job started at $(date)"
echo "Job Directory: $JOB_DIR/$JOBID"

# ==============================================================================
# SECTION 2: SCIENTIFIC LOGIC AND DATABASE SELECTION
# ==============================================================================

if [ -d "$SEARCH_ML_HOME" ]; then
    export data_path="$SEARCH_ML_HOME/datasets"
    export scripts_path="$SEARCH_ML_HOME/scripts"
    export models="$SEARCH_ML_HOME/models"
    export parameters="$SEARCH_ML_HOME/parameters"
else
    echo "Error: SEARCH-ML directory not found at $SEARCH_ML_HOME."
    exit 1
fi

# --- 1. Check PDB file existence ---
echo "Checking PDB file existence: $PDB_FILE"
if [ ! -f "../$PDB_FILE" ]; then
    echo "Error: PDB file '$PDB_FILE' not found."
    exit 2
else
    echo "PDB file '$PDB_FILE' found."
    cp "../$PDB_FILE" ./
    aa_count=$(grep '^ATOM' "$PDB_FILE" | grep -E "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" | grep ' CA ' | cut -c 18-26 | sort -u | wc -l)
    if [ "$aa_count" -gt 24 ]; then
        echo "Protein length ($aa_count) is OK."
    else
        echo "Error: Protein size must be >= 25 amino acids. Found: $aa_count."
        exit 3
    fi
fi

# --- 2. Validate and select database ---
echo "Validating database selection: $DB"
database="" # Initialize the variable

case "$DB" in
  "DB") 
    database="drugbank" 
    ;;
  "FDA") 
    database="fda" 
    ;;
  *)
    echo "Error: Invalid database selection. Argument 3 must be 'DB' or 'FDA'."
    exit 6
    ;;
esac

echo "Selection '$DB' is valid. Selected Database path: $database"

# --- Define the core scientific workflow as a function ---
run_screening_logic() {
    # This function expects to be run inside a directory that has the data linked.
    echo "Preparing PDB and ligand files..."
    grep '^ATOM' "$PDB_FILE" | grep -E "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" > PROT.pdb
    grep -E "^ATOM|^HETATM" "$PDB_FILE" | grep "$LIGAND_CODE" >> PROT.pdb
    LIG=$(echo "$LIGAND_CODE" | awk '{print $1}')
    echo "$LIG" >> lig

    echo "Loading Intel module and running tleap..."
    #module load intel/2018.3.222
    tleap -f leap.cmd
    #module purge
    echo "tleap finished."

    echo "Running Python scripts for feature calculation..."
    python protein_features.py PROT_LIG.pdb "$LIG"
    python pocket_features.py PROT_LIG.pdb "$LIG"
    echo "Feature calculation finished."

    echo "Merging Protein Parameters..."
    python merge_protein_features.py
    echo "Merging finished."

    echo "Starting final virtual screening..."
    python screening.py > results.txt
    echo "Virtual screening finished."
}

# --- 3. Execute workflow (Serial Mode) ---
echo "Linking files for serial execution..."

# Link database files
ln -s "$data_path/$database"/* ./

# Link models and scripts
ln -s "$models"/* ./
ln -s "$scripts_path"/* ./

# Link the parameter file from the parent directory
ln -s "$parameters"/* .

# Execute logic
run_screening_logic

# ==============================================================================
# SECTION 3: FINAL STATUS REPORTING
# ==============================================================================

echo "Job completed successfully at $(date)."

touch COMPLETED
rm -f RUNNING
date
exit 0
