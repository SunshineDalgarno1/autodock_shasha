```markdown
# 🧬 Universal AutoDock GNINA Pipeline
**A Smart & Dynamic Molecular Docking Tool integrating P2Rank and GNINA.**

![Python](https://img.shields.io/badge/Python-3.6%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-success)

This pipeline automates the entire molecular docking process. It dynamically predicts binding pockets using **P2Rank**, evaluates multiple pockets automatically, and performs high-accuracy docking using **GNINA**. It also features a "Smart Bio-Physics Filter" to accurately distinguish between super-binders (large, highly flexible drugs like Tetracyclines) and normal binders.

---

## ✨ Key Features
1. **🤖 Dynamic Pocket Detection:** Automatically detects and evaluates binding pockets. Only pockets with a probability of `> 0.8` are selected for docking. Supports both single-pocket and multi-pocket (e.g., homodimer) targets effortlessly.
2. **🧠 Smart Bio-Physics Filter:** Intelligently filters results based on physical constraints. It tolerates higher intramolecular strain (+8.5) for "Super-binders" (Vina < -9.0 kcal/mol) while maintaining strict strain limits for standard binders, effectively removing false positives.
3. **📦 Auto-Install P2Rank:** No need to manually install P2Rank! The script checks if it exists; if not, it automatically downloads and extracts the official release.
4. **🚀 Highly Optimized Defaults:** Pre-configured with optimal parameters (`Exhaustiveness 64`, `RMSD 2.0`, `Grid Size 20`) for publication-quality results while saving CPU time.
5. **🪶 Lightweight Core:** The core script relies **entirely on Python Standard Libraries**. No `pip install` required for the script itself!

---

## 🛠️ Prerequisites & Installation

While the Python script uses standard libraries, the pipeline relies on external bioinformatics tools. We recommend using **Conda** for a clean setup.

### 1. Clone the Repository
```bash
git clone [https://github.com/SunshineDalgarno1/autodock_shasha.git](https://github.com/SunshineDalgarno1/autodock_shasha.git)
cd autodock_shasha
```

### 2. Create the Conda Environment
Use the provided `environment.yml` to install necessary dependencies (Python, OpenBabel, Wget, Tar):
```bash
conda env create -f environment.yml
conda activate gnina_pipeline
```

### 3. Install GNINA
GNINA must be installed and accessible in your system's PATH. 
- Download the pre-compiled binary from the [Official GNINA GitHub](https://github.com/gnina/gnina).
- Or install it via your preferred method depending on your GPU (CUDA) setup.

---

## 🚀 Usage Guide

The script is designed to be plug-and-play. Make sure your environment is activated before running.

### Basic Usage (Quick Run)
For a standard run using default high-accuracy parameters:
```bash
python autodock_shasha.py -r target_protein.pdb -l ligands.csv
```
*Note: Your `ligands.csv` must contain at least two columns with headers: `Name` and `SMILES`.*

### Advanced Usage
You can customize grid sizes, exhaustiveness, and pose diversity constraints:
```bash
python autodock_shasha.py -r target_protein.pdb -l ligands.csv -sx 20 -sy 20 -sz 20 -ex 64 -no 100 --rmsd 2.0
```

---

## ⚙️ Command-Line Arguments

### 📌 Required Arguments
| Argument | Flag | Description |
| :--- | :--- | :--- |
| `--receptor` | `-r` | Path to the target protein file (`.pdb`). |
| `--ligands` | `-l` | Path to the ligands list (`.csv`). |

### 📦 Grid Box & Precision Parameters
| Argument | Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `--size_x` | `-sx` | `25` | Grid box size in the X dimension (Å). |
| `--size_y` | `-sy` | `25` | Grid box size in the Y dimension (Å). |
| `--size_z` | `-sz` | `25` | Grid box size in the Z dimension (Å). |
| `--exhaustiveness` | `-ex` | `16` | Docking search depth. Use `8` for quick tests, `64` for final publications. |
| `--num_modes` | `-no` | `10` | Number of docking poses generated per pocket. |
| `--rmsd` | | `2.0` | Minimum RMSD threshold to force pose diversity. Highly recommended for large/flexible ligands. |

### 🔧 System Paths
| Argument | Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `--p2rank_path` | | `~/p2rank_2.4.1/prank` | Path to the P2Rank executable. *If not found locally, it will be downloaded automatically.* |

---

## 📊 Output Files
All results will be safely stored in a dynamically generated, timestamped directory (e.g., `results_universal_YYYYMMDD_HHMMSS/`).

1. **`p2rank_output/`**: Contains raw pocket predictions and probabilities from P2Rank.
2. **`*_docked.pdbqt`**: The 3D coordinates of the docked ligands.
3. **`summary_all_modes.csv`**: Raw data of all generated poses across all valid pockets.
4. **`summary_best_binders_smart.csv`**: The crown jewel. This file contains the best pose for each drug, rigorously filtered by the Smart Bio-Physics logic (Vina Affinity, Intramolecular Strain, CNN Pose, and CNN Affinity).

---

## 📝 Author
Developed by **Parinya Tipanyo (Shasha)**.
```
