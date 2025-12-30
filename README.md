# 🧬 Search-ML: Automated Virtual Screening Pipeline

![Version](https://img.shields.io/badge/version-1.3.1-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)
![License](https://img.shields.io/badge/license-AFL-green.svg)
![Lab](https://img.shields.io/badge/SCFBio-IIT%20Delhi-orange.svg)

**Search-ML** is a streamlined, serial execution pipeline designed to automate the virtual screening of protein-ligand complexes against high-value databases like **DrugBank** and **FDA** approved compounds.

Developed at the **Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio)**, this tool integrates structure preparation, feature extraction, and machine learning predictions into a single command.

---

## 🚀 Key Features

* **Automated Pre-processing:** Extracts protein and ligand coordinates and cleans PDB files automatically.
* **AmberTools Integration:** Seamlessly runs `tleap` for system parameterization and topology generation.
* **Smart Feature Extraction:** Calculates physicochemical properties, pocket volume, and Solvent Accessible Surface Area (SASA) using RDKit and custom modules.
* **Targeted Screening:** Optimized for **DrugBank** and **FDA** databases.
* **Robust Error Handling:** Includes detailed logging (`job.log`) and automated environment configuration.

---

## 🛠️ Prerequisites

Before installing, ensure you have the following:

1.  **Operating System:** Linux, macOS, or Windows (via WSL2).
2.  **Conda Manager:** **Miniconda** or **Anaconda** must be installed.
    * *If you don't have it, download it here: [Miniconda Documentation](https://docs.anaconda.com/miniconda/)*
3.  **Git:** To clone the repository.

---

## 📥 Installation Guide

We provide an automated installer (`install.sh`) that sets up the Conda environment and configures necessary paths for you.
Because this script relies on specific dataset locations, you **must** configure the paths before the first run.

1.  Open `master.sh` in your text editor.
2.  Navigate to **Section 2 (Lines 75-79)**.
3.  Update the `export HOME` variable to point to your local installation directory:

### Step 1: Clone the Repository
Open your terminal and run:
```bash
git clone https://github.com/SimulatedLife/search-ml.git
cd search-ml
```

# Inside master.sh

# Ensure these subdirectories exist inside that path:
# - /datasets (containing 'drugbank' and 'fda' folders)
# - /scripts
# - /models
```

---

## 📥 Installation Guide

We provide an automated installer (`install.sh`) that sets up the Conda environment and configures necessary paths for you.

### Step 1: Clone the Repository
Open your terminal and run:
```bash
git clone https://github.com/dheeraj-scfbio/search-ml.git
cd search-ml
```

### Step 2: Download the Model (Crucial)
The pre-trained Random Forest model (~800MB) is too large for GitHub and is hosted externally.

1.  **[Click here to download the model](https://scfbio.iitd.ac.in/search-ml/download.php)**
2.  Move the downloaded file (`final_rf_model.joblib`) into the `models/` folder:
    ```bash
    mv ~/Downloads/final_rf_model.joblib ./models/
    ```

### Step 3: Run the Automated Installer
This script will detect your Conda installation, create the `search-ml` environment, install all dependencies (AmberTools, RDKit, etc.), and configure the project path.

```bash
chmod +x install.sh
./install.sh
```

### Step 4: Activate the Environment
Once the installer finishes, activate the environment to start using the tool:

```bash
conda env config vars set SEARCH_ML="$(pwd)" --name search-ml
conda activate search-ml
```

---

## 📖 Usage Guide

This script is designed for **serial execution** on a Master Node or Workstation.

### Syntax
```bash
./master.sh <PDB_FILENAME> <LIGAND_CODE> <DATABASE>
./master.sh <PDB_FILENAME> <LIGAND_CODE> <MODE>
```

### Argument Breakdown

| Position | Argument | Description |
| :--- | :--- | :--- |
| **$1** | `PDB_FILENAME` | The name of your PDB file **with the extension**. <br>*(e.g., Use `1abc.pdb`)* |
| **$2** | `LIGAND_CODE` | The 3-letter residue name of the ligand inside the PDB. <br>*(e.g., `LIG`, `MOL`, `DRG`)* |
| **$3** | `DATABASE` | The target database identifier. <br>**Options:** `DB` (DrugBank) or `FDA` (FDA Approved) |
| **$3** | `MODE` | The screening mode. <br>**Options:** DB, FDA, SINGLE, CUSTOM |

### Example

**Screening against DrugBank**
Screening protein `9kte.pdb` containing ligand `9UM`:
```bash
bash master.sh 9kte.pdb 9UM DB
bash master.sh 9kte 9UM DB
```

For detailed usage instructions, please read [USAGE.md](USAGE.md).

---

## 🪟 Windows Installation Guide (WSL2)

Since this tool relies on **AmberTools** (Unix-based software), it **cannot** run directly in the Windows Command Prompt or PowerShell.

Instead, you must run it inside **WSL2 (Windows Subsystem for Linux)**. This allows you to run a full Ubuntu environment seamlessly inside Windows.

### Step 1: Install WSL2
1.  Open **PowerShell** as Administrator (Right-click Start Menu -> Terminal (Admin) / PowerShell (Admin)).
2.  Run the following command:
    ```powershell
    wsl --install
    ```
3.  **Restart your computer** when prompted.
4.  After restarting, open the **"Ubuntu"** app from your Start Menu.
5.  Follow the on-screen prompts to create a Username and Password for your Linux system.

### Step 2: Set Up Prerequisites
Inside your new Ubuntu terminal, run the following commands to update the system and install Git:

```bash
sudo apt-get update
sudo apt-get install git -y
```

### Step 3: Clone the Repository
Now you can download the Search-ML code:

```bash
git clone "https://github.com/dheeraj-scfbio/search-ml.git"
cd search-ml
```

### Step 4: Download the Model
*Note: You cannot drag-and-drop files easily into the command line, so we will move it from your Windows Downloads folder.*

1.  **[Click here to download the model](https://scfbio.iitd.ac.in/search-ml/download.php)** to your Windows **Downloads** folder.
2.  In the Ubuntu terminal, move the file from Windows to the project folder:
    ```bash
    # Replace 'YourWindowsUsername' with your actual Windows user name
    mv /mnt/c/Users/YourWindowsUsername/Downloads/final_rf_model.joblib ./models/
    ```

### Step 5: Run the Installer
Now run the automated setup script. This works exactly the same as on Linux.

```bash
chmod +x install.sh
./install.sh
```

### Step 6: Activate and Run
Once the installation finishes:

```bash
conda env config vars set SEARCH_ML="$(pwd)" --name search-ml
conda activate search-ml
./master.sh 9kte.pdb 9UM DB
./master.sh 9kte 9UM DB
```
