#!/bin/bash
# setup_project.sh
# Helper script to set up SIMC C++ project structure

echo "================================================"
echo "  SIMC C++/ROOT Project Setup"
echo "================================================"
echo ""

# Create directory structure
echo "Creating directory structure..."
mkdir -p include/simc
mkdir -p src/core
mkdir -p src/io
mkdir -p data/config
mkdir -p docs
mkdir -p tests
mkdir -p build

echo "✓ Directories created"
echo ""

# List expected files
echo "Expected file structure:"
echo ""
echo "Root directory:"
echo "  - CMakeLists.txt"
echo "  - README.md"
echo "  - setup_project.sh (this script)"
echo ""
echo "include/simc/:"
echo "  - SimcConstants.h"
echo "  - SimcTypes.h"
echo "  - SimcEvent.h"
echo "  - SimcEventLinkDef.h"
echo "  - ConfigManager.h"
echo "  - RandomGenerator.h"
echo "  - OutputManager.h"
echo ""
echo "src/:"
echo "  - CMakeLists.txt"
echo "  - main.cpp"
echo "  - core/SimcEvent.cpp"
echo "  - core/ConfigManager.cpp"
echo "  - core/RandomGenerator.cpp"
echo "  - io/OutputManager.cpp"
echo ""
echo "tests/:"
echo "  - CMakeLists.txt"
echo "  - test_event.cpp"
echo "  - test_config.cpp"
echo "  - test_random.cpp"
echo ""
echo "data/config/:"
echo "  - default.json"
echo ""
echo "docs/:"
echo "  - phase2.tex"
echo ""

# Check which files exist
echo "Checking for files..."
missing_files=0

check_file() {
    if [ -f "$1" ]; then
        echo "  ✓ $1"
    else
        echo "  ✗ $1 (MISSING)"
        ((missing_files++))
    fi
}

# Check root files
check_file "CMakeLists.txt"
check_file "README.md"

# Check headers
check_file "include/simc/SimcConstants.h"
check_file "include/simc/SimcTypes.h"
check_file "include/simc/SimcEvent.h"
check_file "include/simc/SimcEventLinkDef.h"
check_file "include/simc/ConfigManager.h"
check_file "include/simc/RandomGenerator.h"
check_file "include/simc/OutputManager.h"

# Check source files
check_file "src/CMakeLists.txt"
check_file "src/main.cpp"
check_file "src/core/SimcEvent.cpp"
check_file "src/core/ConfigManager.cpp"
check_file "src/core/RandomGenerator.cpp"
check_file "src/io/OutputManager.cpp"

# Check test files
check_file "tests/CMakeLists.txt"
check_file "tests/test_event.cpp"
check_file "tests/test_config.cpp"
check_file "tests/test_random.cpp"

# Check config
check_file "data/config/default.json"

# Check docs
check_file "docs/phase2.tex"

echo ""
if [ $missing_files -eq 0 ]; then
    echo "✓ All files present!"
    echo ""
    echo "Ready to build. Run:"
    echo "  cd build"
    echo "  cmake .."
    echo "  make -j4"
else
    echo "⚠ Warning: $missing_files file(s) missing"
    echo ""
    echo "Please ensure all files are in their correct locations."
    echo "Refer to README.md for the complete file list."
fi

echo ""
echo "================================================"
