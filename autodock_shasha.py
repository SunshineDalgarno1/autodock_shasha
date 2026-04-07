#!/usr/bin/env python3
"""
autodock_shasha.py - Universal Automated Binding Pocket Prediction & Molecular Docking Pipeline
Integrates P2Rank (Auto-Install) for dynamic pocket detection and GNINA for high-accuracy docking.
"""

import argparse
import os
import subprocess
import csv
from datetime import datetime

def print_header():
    print("""
    =============================================================
     🧬 Universal AutoDock GNINA Pipeline (Smart & Dynamic Mode)
    =============================================================
    """)

def check_and_install_p2rank(p2rank_path):
    """Check if P2Rank exists, if not, auto-download and extract it."""
    expanded_path = os.path.expanduser(p2rank_path)
    
    # 1. Check p2rank
    if os.path.exists(expanded_path):
        print(f"✅ Found P2Rank at: {expanded_path}")
        return expanded_path
        
    print(f"⚠️ P2Rank not found at {expanded_path}.")
    print("⏳ Auto-downloading P2Rank 2.4.1 from GitHub...")
    
    # 2. If not dowload it
    try:
        # Download p2rank
        subprocess.run(["wget", "https://github.com/rdk/p2rank/releases/download/2.4.1/p2rank_2.4.1.tar.gz"], check=True)
        print("📦 Extracting P2Rank...")
        # Extraction
        subprocess.run(["tar", "-xzf", "p2rank_2.4.1.tar.gz"], check=True)
        
        # Clean
        if os.path.exists("p2rank_2.4.1.tar.gz"):
            os.remove("p2rank_2.4.1.tar.gz")
            
        new_path = os.path.abspath("p2rank_2.4.1/prank")
        print(f"✅ P2Rank successfully installed at: {new_path}")
        return new_path
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during P2Rank download/extraction: {e}")
        exit(1)

def main():
    print_header()
    
    # 1. Setup Arguments (Optimized Defaults for High Accuracy)
    parser = argparse.ArgumentParser(description="Universal P2Rank + GNINA Molecular Docking Pipeline")
    parser.add_argument("-r", "--receptor", required=True, help="Path to receptor .pdb file")
    parser.add_argument("-l", "--ligands", required=True, help="Path to ligands .csv file (Columns: Name, SMILES)")
    parser.add_argument("-sx", "--size_x", type=float, default=25, help="Grid size X (Default: 25)")
    parser.add_argument("-sy", "--size_y", type=float, default=25, help="Grid size Y (Default: 25)")
    parser.add_argument("-sz", "--size_z", type=float, default=25, help="Grid size Z (Default: 25)")
    parser.add_argument("-ex", "--exhaustiveness", type=int, default=16, help="Exhaustiveness (Default: 16)")
    parser.add_argument("-no", "--num_modes", type=int, default=10, help="Number of modes per pocket (Default: 10)")
    parser.add_argument("--rmsd", type=float, default=2.0, help="Min RMSD filter for pose diversity (Default: 2.0)")
    parser.add_argument("--p2rank_path", type=str, default="~/p2rank_2.4.1/prank", help="Default Path to P2Rank executable")
    
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = f"results_universal_{timestamp}"
    os.makedirs(results_dir, exist_ok=True)
    print(f"📂 Created results directory: {results_dir}")

    receptor_basename = os.path.splitext(os.path.basename(args.receptor))[0]

    # =========================================================
    # 2. Check & Install P2Rank
    # =========================================================
    p2rank_exec = check_and_install_p2rank(args.p2rank_path)

    # =========================================================
    # 3. Run P2Rank and Detect Dynamic Pockets (Prob > 0.8)
    # =========================================================
    print(f"\n🔍 Running P2Rank to map binding pockets...")
    p2rank_out_dir = os.path.join(results_dir, "p2rank_output")
    subprocess.run([p2rank_exec, "predict", "-f", args.receptor, "-o", p2rank_out_dir], capture_output=True)
    
    csv_path = os.path.join(p2rank_out_dir, f"{receptor_basename}.pdb_predictions.csv")
    valid_pockets = {}
    
    if os.path.exists(csv_path):
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            header = [h.strip() for h in next(reader)]
            try:
                name_idx = header.index('name')
                prob_idx = header.index('probability')
                cx_idx, cy_idx, cz_idx = header.index('center_x'), header.index('center_y'), header.index('center_z')
                
                for row in reader:
                    p_name = row[name_idx].strip()
                    p_prob = float(row[prob_idx].strip())
                    
                    # 💡 DYNAMIC POCKET DETECTION: Only keep pockets with > 80% probability
                    if p_prob > 0.8:
                        valid_pockets[p_name] = {
                            'coords': (float(row[cx_idx].strip()), float(row[cy_idx].strip()), float(row[cz_idx].strip())),
                            'probability': p_prob
                        }
            except ValueError as e:
                print(f"❌ Error parsing P2Rank columns: {e}")
                return
                
    if not valid_pockets:
        print("❌ No pockets found with Probability > 0.8. Aborting docking.")
        return
        
    print(f"✅ Found {len(valid_pockets)} highly probable pocket(s) (Probability > 0.8):")
    for pkt, data in valid_pockets.items():
        print(f"   - {pkt.upper()}: Prob = {data['probability']:.3f} | Coords = {data['coords']}")

    # =========================================================
    # 4. Prepare Receptor (PDBQT)
    # =========================================================
    print(f"\n⚙️ Preparing Receptor: {receptor_basename}")
    receptor_pdbqt = os.path.join(results_dir, f"{receptor_basename}.pdbqt")
    subprocess.run(["obabel", args.receptor, "-O", receptor_pdbqt, "-xr"], capture_output=True)

    # =========================================================
    # 5. Docking Loop (Iterate through Ligands & Valid Pockets)
    # =========================================================
    all_modes_data = []
    
    with open(args.ligands, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        next(reader) # Skip header
        
        for row in reader:
            if len(row) < 2: continue
            raw_name, smiles = row[0].strip(), row[1].strip()
            safe_name = raw_name.replace(" ", "_")
            
            print(f"\n=========================================")
            print(f"🚀 Processing Ligand: {raw_name}")
            
            lig_pdbqt = os.path.join(results_dir, f"{safe_name}_ligand.pdbqt")
            subprocess.run(["obabel", f"-:{smiles}", "-O", lig_pdbqt, "--gen3d", "-p", "7.4"], capture_output=True)
            
            if not os.path.exists(lig_pdbqt):
                print(f"⚠️ Failed to generate 3D structure for {raw_name}. Skipping...")
                continue
            
            for target_pocket, data in valid_pockets.items():
                coords = data['coords']
                print(f"   🧬 Docking into {target_pocket.upper()}...")
                
                out_pdbqt = os.path.join(results_dir, f"{safe_name}_{target_pocket}_docked.pdbqt")
                log_file = os.path.join(results_dir, f"{safe_name}_{target_pocket}_gnina.log")
                
                # 💡 Advanced GNINA params
                gnina_cmd = [
                    "gnina", "-r", receptor_pdbqt, "-l", lig_pdbqt,
                    "--center_x", str(coords[0]), "--center_y", str(coords[1]), "--center_z", str(coords[2]),
                    "--size_x", str(args.size_x), "--size_y", str(args.size_y), "--size_z", str(args.size_z),
                    "--exhaustiveness", str(args.exhaustiveness), 
                    "--num_modes", str(args.num_modes),
                    "--min_rmsd_filter", str(args.rmsd),
                    "--cnn_scoring", "rescore",
                    "--out", out_pdbqt
                ]
                
                with open(log_file, "w") as out_log:
                    subprocess.run(gnina_cmd, stdout=out_log, stderr=subprocess.STDOUT)
                
                # Parse GNINA output
                try:
                    with open(log_file, "r") as log:
                        lines = log.readlines()
                        for i, line in enumerate(lines):
                            if "-----+------------+------------+------------+----------" in line:
                                for j in range(i + 1, len(lines)):
                                    pose = lines[j].split()
                                    if len(pose) < 5: break
                                    
                                    all_modes_data.append({
                                        'Drug_Name': raw_name,
                                        'Target_Pocket': target_pocket.upper(),
                                        'Mode': int(pose[0]),
                                        'Vina_Affinity': float(pose[1]),
                                        'Intramol_Energy': float(pose[2]),
                                        'CNN_Pose_Score': float(pose[3]),
                                        'CNN_Affinity': float(pose[4])
                                    })
                                break
                except Exception as e:
                    print(f"⚠️ Error reading log for {raw_name} in {target_pocket}: {e}")

    # =========================================================
    # 6. Summarize & Apply Smart Bio-Physics Filter
    # =========================================================
    if not all_modes_data: return

    all_modes_csv = os.path.join(results_dir, "summary_all_modes.csv")
    smart_modes_csv = os.path.join(results_dir, "summary_best_binders_smart.csv")
    
    # 6.1 Save ALL modes raw data
    with open(all_modes_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=all_modes_data[0].keys())
        writer.writeheader()
        writer.writerows(all_modes_data)
        
    # 6.2 Smart Filter Logic (The Physicist's Approach)
    best_modes = {}
    for row in all_modes_data:
        name = row['Drug_Name']
        
        base_pass = (row['Vina_Affinity'] <= -7.0 and 
                     row['CNN_Pose_Score'] >= 0.3 and 
                     row['CNN_Affinity'] >= 5.0)
                     
        # Dynamic Intramol Penalty
        if row['Vina_Affinity'] <= -9.0:
            intramol_pass = row['Intramol_Energy'] <= 8.5  # Super-binders can fold under stress
        else:
            intramol_pass = row['Intramol_Energy'] <= 3.0  # Normal binders must be relaxed
            
        if base_pass and intramol_pass:
            # Save the highest CNN_Affinity across ALL pockets for this drug
            if name not in best_modes or row['CNN_Affinity'] > best_modes[name]['CNN_Affinity']:
                best_modes[name] = row
            
    # Sort best hits from highest CNN affinity to lowest
    sorted_best = sorted(best_modes.values(), key=lambda x: x['CNN_Affinity'], reverse=True)
    
    # 6.3 Save the Filtered list
    if sorted_best:
        with open(smart_modes_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=sorted_best[0].keys())
            writer.writeheader()
            writer.writerows(sorted_best)
        print(f"\n🥇 Smartly filtered poses successfully saved to: {smart_modes_csv}")
    else:
        print(f"\n⚠️ WARNING: No poses passed the strict criteria! Try checking {all_modes_csv}")

    print(f"\n🎉 Pipeline Complete! Results safely stored in: {results_dir}/")

if __name__ == "__main__":
    main()
