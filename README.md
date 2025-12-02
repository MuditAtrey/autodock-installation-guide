# AutoDock Vina Environment Installation Guide for M4 MacBook Air
## Step-by-step instructions for Claude Sonnet 4.5

---

## Prerequisites Check

### Step 1: Verify System Requirements
```bash
# Check macOS version
sw_vers

# Check architecture (should show arm64)
uname -m

# Check if Conda is installed
conda --version
```

**Expected Results:**
- macOS 13.0 or later
- Architecture: arm64 (Apple Silicon)
- Conda should be installed (any version 4.x+)

**If Conda is NOT installed:**
```bash
# Download Miniforge (optimized for Apple Silicon)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh

# Install Miniforge
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3

# Initialize shell
$HOME/miniforge3/bin/conda init zsh

# Restart terminal or source config
source ~/.zshrc

# Verify installation
conda --version
```

---

## Installation Process

### Step 2: Create Conda Environment
```bash
# Create new environment with Python 3.11
conda create -n autodock python=3.11 -y

# Activate the environment
conda activate autodock

# Verify Python version
python --version
```

**Expected Output:** `Python 3.11.x`

---

### Step 3: Install AutoDock Vina
```bash
# Install Vina from conda-forge
conda install -c conda-forge vina -y

# Verify installation
vina --version
```

**Expected Output:** `AutoDock Vina v1.2.5` (or later)

---

### Step 4: Install PyMOL
```bash
# Install PyMOL open-source version
conda install -c conda-forge pymol-open-source -y

# Test PyMOL command-line mode (quick test)
pymol -c -e -Q
```

**Expected Output:** PyMOL should start and exit without errors

---

### Step 5: Install Meeko (Ligand Preparation Tool)
```bash
# Install Meeko
conda install -c conda-forge meeko -y

# Verify Meeko CLI is available
mk_prepare_ligand.py --help | head -5
```

**Expected Output:** Usage information for mk_prepare_ligand.py

---

### Step 6: Install RDKit and Scientific Libraries
```bash
# Install RDKit and additional scientific packages
conda install -c conda-forge rdkit numpy pandas scipy matplotlib -y

# Verify RDKit installation
python -c "import rdkit; print('RDKit version:', rdkit.__version__)"
```

**Expected Output:** `RDKit version: 2025.03.x` (or later)

---

## Verification Tests

### Step 7: Test AutoDock Vina
```bash
conda activate autodock

# Show help menu
vina --help

# Verify it shows usage information with options like:
# --receptor, --ligand, --center_x, --size_x, etc.
```

**Success Criteria:** Help menu displays without errors

---

### Step 8: Test PyMOL Python API
```bash
conda activate autodock

python -c "
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])
cmd.fragment('gly')
natoms = cmd.count_atoms()
print(f'PyMOL test: Created glycine with {natoms} atoms')
cmd.delete('all')
pymol.cmd.quit()
print('PyMOL API test PASSED')
"
```

**Expected Output:** 
```
PyMOL test: Created glycine with 7 atoms
PyMOL API test PASSED
```

---

### Step 9: Test RDKit Functionality
```bash
conda activate autodock

python -c "
from rdkit import Chem
from rdkit.Chem import AllChem

# Create aspirin molecule
mol = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

print(f'RDKit test: Created aspirin with {mol.GetNumAtoms()} atoms')
print('3D coordinates generated successfully')
print('RDKit test PASSED')
"
```

**Expected Output:**
```
RDKit test: Created aspirin with 21 atoms
3D coordinates generated successfully
RDKit test PASSED
```

---

### Step 10: Test Meeko
```bash
conda activate autodock

python -c "
from meeko import MoleculePreparation
print('Meeko imports successfully')
print('Meeko test PASSED')
"
```

**Expected Output:**
```
Meeko imports successfully
Meeko test PASSED
```

---

### Step 11: Test Scientific Libraries
```bash
conda activate autodock

python -c "
import numpy as np
import pandas as pd
import scipy
import matplotlib

arr = np.array([1, 2, 3, 4, 5])
mean = np.mean(arr)
df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})

print(f'NumPy mean: {mean}')
print(f'Pandas DataFrame shape: {df.shape}')
print(f'SciPy version: {scipy.__version__}')
print('Scientific libraries test PASSED')
"
```

**Expected Output:**
```
NumPy mean: 3.0
Pandas DataFrame shape: (3, 2)
SciPy version: 1.16.x
Scientific libraries test PASSED
```

---

### Step 12: Launch PyMOL GUI
```bash
conda activate autodock

# Launch PyMOL with graphical interface
pymol &
```

**Success Criteria:** PyMOL window opens without crashes

**To close:** Use File > Quit in PyMOL or press Cmd+Q

---

## Comprehensive Verification

### Step 13: Run Complete Integration Test
```bash
conda activate autodock

python -c "
import sys
print('='*60)
print('AUTODOCK ENVIRONMENT VERIFICATION')
print('='*60)
print()

# Test all components
tests_passed = []
tests_failed = []

# Test 1: Vina
try:
    import subprocess
    result = subprocess.run(['vina', '--version'], capture_output=True, text=True)
    if 'AutoDock Vina' in result.stdout:
        tests_passed.append('AutoDock Vina')
    else:
        tests_failed.append('AutoDock Vina')
except:
    tests_failed.append('AutoDock Vina')

# Test 2: PyMOL
try:
    import pymol
    tests_passed.append('PyMOL')
except:
    tests_failed.append('PyMOL')

# Test 3: Meeko
try:
    from meeko import MoleculePreparation
    tests_passed.append('Meeko')
except:
    tests_failed.append('Meeko')

# Test 4: RDKit
try:
    from rdkit import Chem
    mol = Chem.MolFromSmiles('CCO')
    if mol.GetNumAtoms() == 3:
        tests_passed.append('RDKit')
    else:
        tests_failed.append('RDKit')
except:
    tests_failed.append('RDKit')

# Test 5: NumPy
try:
    import numpy as np
    arr = np.array([1, 2, 3])
    if np.mean(arr) == 2.0:
        tests_passed.append('NumPy')
    else:
        tests_failed.append('NumPy')
except:
    tests_failed.append('NumPy')

# Test 6: Pandas
try:
    import pandas as pd
    df = pd.DataFrame({'A': [1]})
    tests_passed.append('Pandas')
except:
    tests_failed.append('Pandas')

# Test 7: SciPy
try:
    import scipy
    tests_passed.append('SciPy')
except:
    tests_failed.append('SciPy')

# Test 8: Matplotlib
try:
    import matplotlib
    tests_passed.append('Matplotlib')
except:
    tests_failed.append('Matplotlib')

# Report results
print('PASSED TESTS:')
for test in tests_passed:
    print(f'  ✓ {test}')
print()

if tests_failed:
    print('FAILED TESTS:')
    for test in tests_failed:
        print(f'  ✗ {test}')
    print()
    print('STATUS: INSTALLATION INCOMPLETE')
else:
    print('STATUS: ALL TESTS PASSED ✓')
    print()
    print('Environment is ready for:')
    print('  • Molecular docking with AutoDock Vina')
    print('  • Visualization with PyMOL')
    print('  • Ligand preparation with Meeko')
    print('  • Cheminformatics with RDKit')
    print('  • Data analysis')

print('='*60)
"
```

**Expected Output:**
```
============================================================
AUTODOCK ENVIRONMENT VERIFICATION
============================================================

PASSED TESTS:
  ✓ AutoDock Vina
  ✓ PyMOL
  ✓ Meeko
  ✓ RDKit
  ✓ NumPy
  ✓ Pandas
  ✓ SciPy
  ✓ Matplotlib

STATUS: ALL TESTS PASSED ✓

Environment is ready for:
  • Molecular docking with AutoDock Vina
  • Visualization with PyMOL
  • Ligand preparation with Meeko
  • Cheminformatics with RDKit
  • Data analysis
============================================================
```

---

## Troubleshooting

### Issue: Conda not found
**Solution:** Install Miniforge as shown in Step 1

### Issue: Package installation fails
**Solution:** 
```bash
# Update conda
conda update -n base conda -y

# Clear cache and retry
conda clean --all -y
conda install -c conda-forge <package_name> -y
```

### Issue: PyMOL crashes or won't open GUI
**Solution:**
```bash
# Ensure XQuartz is NOT running (it can cause conflicts)
pkill -9 XQuartz

# Try launching PyMOL again
conda activate autodock
pymol
```

### Issue: Environment activation not working
**Solution:**
```bash
# Re-initialize conda for zsh
conda init zsh
source ~/.zshrc

# Try activation again
conda activate autodock
```

### Issue: Import errors in Python
**Solution:**
```bash
# Verify you're in the correct environment
conda activate autodock
which python  # Should show path with 'autodock' in it

# Reinstall problematic package
conda install -c conda-forge --force-reinstall <package_name> -y
```

---

## Quick Reference

### Daily Usage Commands
```bash
# Activate environment
conda activate autodock

# Launch PyMOL GUI
pymol

# Run Vina docking
vina --receptor protein.pdbqt --ligand ligand.pdbqt \
     --center_x 0 --center_y 0 --center_z 0 \
     --size_x 20 --size_y 20 --size_z 20 \
     --out result.pdbqt

# Prepare ligand with Meeko
mk_prepare_ligand.py -i molecule.sdf -o ligand.pdbqt

# Deactivate environment
conda deactivate
```

### Check Installed Versions
```bash
conda activate autodock
conda list | grep -E "vina|pymol|meeko|rdkit|numpy|pandas|scipy"
```

---

## Installation Summary

**Total Time:** ~15-20 minutes (depending on internet speed)

**Disk Space Required:** ~2-3 GB

**Components Installed:**
1. Conda environment: `autodock`
2. Python 3.11.x
3. AutoDock Vina 1.2.5+
4. PyMOL 3.1.0+
5. Meeko 0.6.0+
6. RDKit 2025.03.x+
7. NumPy, Pandas, SciPy, Matplotlib (latest versions)

**Architecture:** Native ARM64 (optimized for Apple Silicon M4)

---

## Success Criteria Checklist

Execute this final check to confirm everything is ready:

- [ ] `conda activate autodock` works without errors
- [ ] `vina --version` shows version information
- [ ] `pymol -c -e -Q` runs without errors
- [ ] `python -c "import rdkit; print('OK')"` prints OK
- [ ] `python -c "import meeko; print('OK')"` prints OK
- [ ] `python -c "import numpy, pandas, scipy; print('OK')"` prints OK
- [ ] `pymol` launches GUI successfully (window appears)
- [ ] All verification tests in Step 13 pass

**If all items are checked:** Installation is complete and verified! ✓

---

## Notes for Claude

- Execute steps sequentially in order
- Wait for each command to complete before proceeding
- If any step fails, refer to Troubleshooting section
- Report any unexpected errors with full output
- The entire process should complete without manual intervention
- All commands are designed for zsh shell on macOS
- Test results should match expected outputs exactly

---

**Document Version:** 1.0  
**Last Updated:** December 3, 2025  
**Target Platform:** M4 MacBook Air (Apple Silicon)  
**Tested On:** M3 MacBook (verified compatible)
