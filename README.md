# 🧬 Search-ML: Automated Virtual Screening Pipeline

![Version](https://img.shields.io/badge/version-1.3.1-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)
![Lab](https://img.shields.io/badge/SCFBio-IIT%20Delhi-orange.svg)

**Search-ML** is a streamlined, serial execution pipeline designed to automate the virtual screening of protein-ligand complexes against high-value databases like **DrugBank** and **FDA** approved compounds.

Developed at the **Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio)**, this tool integrates structure preparation, feature extraction, and machine learning predictions into a single command.

---

## 🚀 Key Features

* **Automated Pre-processing:** Extracts protein and ligand coordinates and cleans PDB files automatically.
* **AmberTools Integration:** Seamlessly runs `tleap` for system parameterization.
* **Smart Feature Extraction:** Calculates physicochemical properties, pocket volume, and SASA using custom Python modules.
* **Targeted Screening:** Optimized for **DrugBank** and **FDA** databases.
* **Error Handling:** Robust logging (`job.log`) and error trapping to ensure data integrity.

---

## 🛠️ Prerequisites

Before running the pipeline, ensure you have the following installed on your Master Node or Workstation:

1.  **Bash Shell** (Required for process substitution).
2.  **AmberTools** (specifically `tleap`).
3.  **Python 3.x** with the following dependencies:
    * `numpy`
    * `pandas`
    * `scikit-learn` (or specific ML libraries used in your models)
    * `rdkit` (if used for physicochemical properties)

---

## ⚙️ Configuration (Important)

Because this script relies on specific dataset locations, you **must** configure the paths before the first run.

1.  Open `master.sh` in your text editor.
2.  Navigate to **Section 2 (Lines 75-79)**.
3.  Update the `export HOME` variable to point to your local installation directory:

```bash
# Inside master.sh

# Ensure these subdirectories exist inside that path:
# - /datasets (containing 'drugbank' and 'fda' folders)
# - /scripts
# - /models
```

## 🪟 Windows Support

This tool runs natively on Linux and macOS. For Windows users, we support execution via **WSL2 (Windows Subsystem for Linux)**.

**Step 1: Install WSL**
1. Open PowerShell as Administrator.
2. Run the command: `wsl --install`
3. Restart your computer.
4. Open the newly installed "Ubuntu" app from your Start Menu.

**Step 2: Install in WSL**
Now that you are inside the Ubuntu terminal, the steps are identical to Linux:

1. **Install Git:** `sudo apt-get update && sudo apt-get install git`
2. **Clone the Repo:** `git clone https://github.com/dheeraj-scfbio/search-ml.git`
3. **Run Installer:** `bash install.sh`
4. **Download Model:** (Follow the "Model Download" instructions above)
5. **Run Tool:** `bash master.sh ...`
