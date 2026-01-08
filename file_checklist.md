# SIMC Phase 2 - Complete File Checklist

## How to Download Files

Since the artifact system saves as .txt, you have two options:

### Option 1: Manual (Save and Rename)
1. Click on each artifact
2. Copy the content
3. Save to correct location with correct extension
4. Example: Save "SimcConstants.h" content to `include/simc/SimcConstants.h`

### Option 2: Copy-Paste with Script
1. Run `setup_project.sh` first to create directories
2. Use your text editor to save each file with the correct name

## Complete File List (21 Files Total)

### ✅ Root Directory (4 files)

- [ ] **CMakeLists.txt** 
  - Location: `./CMakeLists.txt`
  - Extension: `.txt` → NO extension change needed
  
- [ ] **README.md**
  - Location: `./README.md`
  - Extension: `.txt` → `.md`
  
- [ ] **setup_project.sh**
  - Location: `./setup_project.sh`
  - Extension: `.txt` → `.sh`
  - Make executable: `chmod +x setup_project.sh`
  
- [ ] **FILE_CHECKLIST.md** (this file)
  - Location: `./FILE_CHECKLIST.md`
  - Extension: `.txt` → `.md`

---

### ✅ Headers - include/simc/ (7 files)

- [ ] **SimcConstants.h**
  - Location: `include/simc/SimcConstants.h`
  - Extension: `.txt` → `.h`
  
- [ ] **SimcTypes.h**
  - Location: `include/simc/SimcTypes.h`
  - Extension: `.txt` → `.h`
  
- [ ] **SimcEvent.h**
  - Location: `include/simc/SimcEvent.h`
  - Extension: `.txt` → `.h`
  
- [ ] **SimcEventLinkDef.h**
  - Location: `include/simc/SimcEventLinkDef.h`
  - Extension: `.txt` → `.h`
  
- [ ] **ConfigManager.h**
  - Location: `include/simc/ConfigManager.h`
  - Extension: `.txt` → `.h`
  
- [ ] **RandomGenerator.h**
  - Location: `include/simc/RandomGenerator.h`
  - Extension: `.txt` → `.h`
  
- [ ] **OutputManager.h**
  - Location: `include/simc/OutputManager.h`
  - Extension: `.txt` → `.h`

---

### ✅ Source Files - src/ (6 files)

- [ ] **src/CMakeLists.txt**
  - Location: `src/CMakeLists.txt`
  - Extension: `.txt` → NO change
  
- [ ] **src/main.cpp**
  - Location: `src/main.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **src/core/SimcEvent.cpp**
  - Location: `src/core/SimcEvent.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **src/core/ConfigManager.cpp**
  - Location: `src/core/ConfigManager.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **src/core/RandomGenerator.cpp**
  - Location: `src/core/RandomGenerator.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **src/io/OutputManager.cpp**
  - Location: `src/io/OutputManager.cpp`
  - Extension: `.txt` → `.cpp`

---

### ✅ Test Files - tests/ (4 files)

- [ ] **tests/CMakeLists.txt**
  - Location: `tests/CMakeLists.txt`
  - Extension: `.txt` → NO change
  
- [ ] **tests/test_event.cpp**
  - Location: `tests/test_event.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **tests/test_config.cpp**
  - Location: `tests/test_config.cpp`
  - Extension: `.txt` → `.cpp`
  
- [ ] **tests/test_random.cpp**
  - Location: `tests/test_random.cpp`
  - Extension: `.txt` → `.cpp`

---

### ✅ Configuration - data/config/ (1 file)

- [ ] **default.json**
  - Location: `data/config/default.json`
  - Extension: `.txt` → `.json`

---

### ✅ Documentation - docs/ (1 file)

- [ ] **phase2.tex**
  - Location: `docs/phase2.tex`
  - Extension: `.txt` → `.tex`

---

## Quick Setup Procedure

### Step 1: Create Directory Structure
```bash
mkdir -p include/simc src/{core,io} data/config docs tests build
```

### Step 2: Save All Files
Download and save each file from the artifacts to its correct location.

### Step 3: Make Setup Script Executable
```bash
chmod +x setup_project.sh
```

### Step 4: Run Setup Script to Verify
```bash
./setup_project.sh
```

This will check if all files are present.

### Step 5: Build
```bash
cd build
cmake ..
make -j4
```

### Step 6: Run
```bash
./simc
```

---

## Extension Mapping Quick Reference

| Download As | Save As |
|-------------|---------|
| `.txt` (CMakeLists) | No change |
| `.txt` (other) | See specific extension |
| Headers | `.h` |
| Source | `.cpp` |
| Markdown | `.md` |
| JSON | `.json` |
| LaTeX | `.tex` |
| Shell script | `.sh` |

---

## Troubleshooting

### Problem: CMake can't find files
**Solution:** Check file names and locations exactly match the structure above.

### Problem: Test compilation errors
**Solution:** Tests require Google Test. If not installed, disable with:
```bash
cmake -DBUILD_TESTS=OFF ..
```

### Problem: ROOT not found
**Solution:** Source ROOT environment:
```bash
source /path/to/root/bin/thisroot.sh
cmake ..
```

### Problem: nlohmann/json not found
**Solution:** CMake will automatically fetch it. Ensure internet connection or install:
```bash
# Ubuntu/Debian
sudo apt-get install nlohmann-json3-dev

# macOS
brew install nlohmann-json
```

---

## Verification Checklist

After setting up, verify:

- [ ] All 21 files present
- [ ] `setup_project.sh` runs without errors
- [ ] `cmake ..` completes successfully
- [ ] `make` compiles without errors
- [ ] `./simc` runs and creates `simc_output.root`
- [ ] `root simc_output.root` can open the file

---

## File Count Summary

- **Total Files:** 21
- **Headers (.h):** 7
- **Source (.cpp):** 10
- **CMake files:** 3
- **Config (.json):** 1
- **Documentation (.tex, .md):** 3
- **Script (.sh):** 1

---

## Next Steps After Setup

1. **Compile LaTeX documentation:**
   ```bash
   cd docs
   pdflatex phase2.tex
   pdflatex phase2.tex  # Run twice for references
   ```

2. **Test the installation:**
   ```bash
   cd build
   ./simc
   root simc_output.root
   ```

3. **Run unit tests (if Google Test installed):**
   ```bash
   cd build
   make test
   ```

4. **Review the code:**
   - Read `README.md` for overview
   - Read `docs/phase2.pdf` for architecture
   - Browse headers in `include/simc/`

---

## Contact for Issues

If you encounter problems:
1. Run `./setup_project.sh` and check output
2. Verify ROOT is properly installed
3. Check CMake version: `cmake --version` (need ≥ 3.15)
4. Check compiler: `g++ --version` (need ≥ 7 for C++17)

---

**Last Updated:** January 2026  
**Phase:** 2 - Core Infrastructure  
**Status:** Complete - Ready for Phase 3
