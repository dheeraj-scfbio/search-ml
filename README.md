# 🧬 Search-ML: Automated Virtual Screening Pipeline

![Version](https://img.shields.io/badge/version-1.3.1-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
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
2.  **AmberTools** (specifically `tleap` must be in your PATH).
3.  **Python 3.8+** with the following dependencies:
    ```bash
    pip install numpy pandas scikit-learn rdkit joblib
    ```

---

## 📥 Model Download (Required)

The Random Forest model (`final_rf_model.joblib`) is too large for GitHub (>100MB) and must be downloaded separately.

1.  **[Click here to download the model](INSERT_YOUR_GOOGLE_DRIVE_LINK_HERE)**
2.  **Setup:**
    * Download the file.
    * Move it into the `models/` directory inside this repository:
        ```bash
        mv ~/Downloads/final_rf_model.joblib ./models/
        ```

---

## ⚙️ Configuration (Important)

Because this script relies on specific dataset locations, you **must** configure the paths before the first run.

1.  Open `run_virtual_screening.sh` in your text editor.
2.  Navigate to **Section 2 (Lines 75-79)**.
3.  Update the `export HOME` variable to point to your local installation directory:

```bash
# Inside run_virtual_screening.sh

# CHANGE THIS PATH to your actual directory
export HOME="/path/to/your/search-ml" 

# Ensure these subdirectories exist inside that path:
# - /datasets (containing 'drugbank' and 'fda' folders)
# - /scripts
# - /models
