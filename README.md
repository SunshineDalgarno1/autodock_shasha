# 🧬 Universal AutoDock GNINA Pipeline
**A Smart & Dynamic Molecular Docking Tool integrating P2Rank and GNINA.**

![Python](https://img.shields.io/badge/Python-3.6%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Active-success)

This pipeline automates the entire molecular docking process. It dynamically predicts binding pockets using **P2Rank**, evaluates multiple pockets automatically, and performs high-accuracy docking using **GNINA**. It also features a "Smart Bio-Physics Filter" to accurately distinguish between super-binders (large, highly flexible drugs) and normal binders.

---

## ✨ Key Features
1. **🤖 Dynamic Pocket Detection:** Automatically detects and evaluates binding pockets. Only pockets with a probability of `> 0.8` are selected for docking. Supports both single-pocket and multi-pocket (e.g., homodimer) targets effortlessly.
2. **🧠 Smart Bio-Physics Filter:** Intelligently filters results based on physical constraints. It tolerates higher intramolecular strain for "Super-binders" (Vina < -9.0 kcal/mol) while maintaining strict strain limits for standard binders, effectively removing false positives.
3. **📦 Auto-Install P2Rank:** No need to manually install P2Rank! The script checks if it exists; if not, it automatically downloads and extracts the official release.
4. **🚀 Highly Optimized Defaults:** Pre-configured with optimal parameters (`Exhaustiveness 16`, `RMSD 2.0`, `Grid Size 25`) for publication-quality results while saving CPU time.
5. **🪶 Lightweight:** The core script relies **entirely on Python Standard Libraries**. No `pip install` required for the script itself!

---

## 🛠️ Prerequisites & Installation

While the Python script uses standard libraries, the pipeline relies on external bioinformatics tools. We recommend using **Conda** for a clean setup.

### 1. Create a Conda Environment (Recommended)
Create a file named `environment.yml` with the following content:
```yaml
name: gnina_pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.8
  - openbabel
  - wget
  - tar