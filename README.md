# SIMC C++/ROOT Monte Carlo - Phase 2

**Physics Monte Carlo for Jefferson Lab Hall C/A Experiments**

This is Phase 2 of the SIMC Fortran to C++/ROOT conversion project, implementing the core infrastructure.

## Project Status

**Phase 2 Complete** - Core infrastructure implemented:
- ✅ Event data structures with ROOT compatibility
- ✅ Configuration management (JSON-based)
- ✅ Random number generation
- ✅ ROOT TTree output system
- ✅ Histogram management
- ✅ CMake build system

**Next:** Phase 3 will implement physics models and event generators.

## Prerequisites

- **CMake** ≥ 3.15
- **ROOT** ≥ 6.20
- **C++17** compatible compiler (g++ ≥ 7, clang++ ≥ 5)
- **nlohmann/json** ≥ 3.9 (automatically fetched if not installed)

### Installing ROOT

**On Ubuntu/Debian:**
```bash
sudo apt-get install root-system libroot-dev
```

**On macOS (with Homebrew):**
```bash
brew install root
```

**From source:**
See https://root.cern/install/build_from_source/

## Building

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake ..

# Build
make -j$(nproc)

# Run tests (optional)
make test

# Install (optional)
sudo make install
```

### Build Options

```bash
# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build (optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build without tests
cmake -DBUILD_TESTS=OFF ..

# Specify ROOT location
cmake -DROOT_DIR=/path/to/root ..
```

## Running SIMC

### Basic Usage

```bash
# Run with default configuration
./simc

# Run with custom configuration
./simc path/to/config.json
```

### Configuration File

Edit `data/config/default.json` to configure:
- Beam energy and properties
- Target material and geometry
- Spectrometer settings
- Monte Carlo parameters
- Output options

Example:
```json
{
  "beam": {
    "energy": 10.6,
    "energy_spread": 0.05
  },
  "target": {
    "type": "LH2",
    "thickness": 0.1
  },
  "monte_carlo": {
    "nevents": 100000,
    "random_seed": 12345
  }
}
```

### Output

SIMC produces a ROOT file (`simc_output.root` by default) containing:
- **T**: Main event tree with all kinematic quantities
- **T_gen**: Generated events (before cuts)
- **T_stats**: Run statistics
- **histograms/**: Diagnostic histograms

### Analyzing Output

```cpp
// Open ROOT file
TFile* f = new TFile("simc_output.root");
TTree* T = (TTree*)f->Get("T");

// Draw histogram
T->Draw("Q2/1e6 >> hQ2(100,0,5)", "success");

// Access event data
simc::SimcEvent* evt = nullptr;
T->SetBranchAddress("SimcEvent", &evt);
T->GetEntry(0);
evt->Print();
```

## Project Structure

```
simc_cpp/
├── CMakeLists.txt              # Main build configuration
├── README.md                   # This file
├── include/simc/               # Public headers
│   ├── SimcConstants.h         # Physical constants
│   ├── SimcTypes.h             # Common type definitions
│   ├── SimcEvent.h             # Event structure
│   ├── ConfigManager.h         # Configuration management
│   ├── RandomGenerator.h       # Random number generation
│   └── OutputManager.h         # ROOT I/O
├── src/                        # Implementation
│   ├── core/                   # Core functionality
│   │   ├── SimcEvent.cpp
│   │   ├── ConfigManager.cpp
│   │   └── RandomGenerator.cpp
│   ├── io/                     # Input/output
│   │   └── OutputManager.cpp
│   └── main.cpp                # Main program
├── data/config/                # Configuration files
│   └── default.json
├── tests/                      # Unit tests
└── docs/                       # Documentation
    └── phase2.tex              # LaTeX documentation
```

## File Locations

Save files in the following locations:

### Headers (`include/simc/`)
- `SimcConstants.h`
- `SimcTypes.h`
- `SimcEvent.h`
- `SimcEventLinkDef.h` (for ROOT dictionary)
- `ConfigManager.h`
- `RandomGenerator.h`
- `OutputManager.h`

### Source Files (`src/`)
- `core/SimcEvent.cpp`
- `core/ConfigManager.cpp`
- `core/RandomGenerator.cpp`
- `io/OutputManager.cpp`
- `main.cpp`

### Build Files
- `CMakeLists.txt` (root directory)
- `src/CMakeLists.txt`

### Data Files
- `data/config/default.json`

## Documentation

Compile the LaTeX documentation:

```bash
cd docs
pdflatex phase2.tex
pdflatex phase2.tex  # Run twice for references
```

This generates `phase2.pdf` with complete architecture documentation.

## Development

### Code Style

- **C++17 standard** throughout
- **Header/source separation** for all classes
- **Doxygen-style** comments for all public APIs
- **snake_case** for variables and functions
- **PascalCase** for classes and types

### Adding New Features

1. Add header to `include/simc/`
2. Add implementation to `src/core/` or `src/io/`
3. Update `src/CMakeLists.txt`
4. Add unit tests to `tests/`
5. Update documentation

## Testing

Run unit tests:
```bash
cd build
make test
```

Or run directly:
```bash
./simc_tests
```

## Troubleshooting

### ROOT Not Found

```bash
# Set ROOT environment
source /path/to/root/bin/thisroot.sh

# Or specify ROOT location
cmake -DROOT_DIR=/path/to/root ..
```

### Compilation Errors

Make sure you have C++17 support:
```bash
g++ --version  # Should be ≥ 7.0
```

### Runtime Errors

Check ROOT file permissions:
```bash
ls -l simc_output.root
```

## Contributing

This is an active development project. Phase 3 will add:
- Physics event generators
- Cross section calculations
- Spectrometer Monte Carlos
- Complete physics models

## License

SIMC is developed at Jefferson Lab for nuclear physics research.

## Contact

For questions about SIMC C++/ROOT:
- Jefferson Lab Hall C: https://hallcweb.jlab.org/
- SIMC Documentation: https://hallcweb.jlab.org/wiki/index.php/SIMC_Monte_Carlo

## Acknowledgments

Based on the original SIMC Fortran code developed by the Jefferson Lab Hall C collaboration.

---

**Version:** 2.0.0-phase2  
**Last Updated:** January 2026









# SIMC C++/ROOT Monte Carlo - Phase 3

**Physics Monte Carlo for Jefferson Lab Hall C/A Experiments**

This is Phase 3 of the SIMC Fortran to C++/ROOT conversion project, implementing physics calculations and event generation.

## Project Status

**Phase 2 Complete** ✅ - Core infrastructure
**Phase 3 In Progress** 🚧 - Physics implementation

- ✅ Event data structures with ROOT compatibility
- ✅ Configuration management (JSON-based)
- ✅ Random number generation
- ✅ ROOT TTree output system
- ✅ Histogram management
- ✅ CMake build system
- ✅ Standalone test framework (no Google Test dependency)

### In Progress (Phase 3):
- 🚧 Kinematic calculations
- 🚧 Cross section framework
- 🚧 Event generators
- 🚧 Radiative corrections
- 🚧 Multiple scattering
- 🚧 Energy loss

### Upcoming:
- ⏳ Spectrometer optics (COSY models)
- ⏳ Spectral functions
- ⏳ Complete physics models

## Prerequisites

- **CMake** ≥ 3.15
- **ROOT** ≥ 6.20
- **C++17** compatible compiler (g++ ≥ 7, clang++ ≥ 5)
- **nlohmann/json** ≥ 3.9 (automatically fetched if not installed)

**Note:** Google Test is NOT required. We use a standalone test framework.

### Installing ROOT

**On Ubuntu/Debian:**
```bash
sudo apt-get install root-system libroot-dev
```

**On macOS (with Homebrew):**
```bash
brew install root
```

**From source:**
See https://root.cern/install/build_from_source/

## Building

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake ..

# Build
make -j$(nproc)

# Run tests (optional, no external dependencies)
make test
# OR
./simc_tests

# Install (optional)
sudo make install
```

### Build Options

```bash
# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build (optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build without tests
cmake -DBUILD_TESTS=OFF ..

# Specify ROOT location
cmake -DROOT_DIR=/path/to/root ..
```

## Running SIMC

### Basic Usage

```bash
# Run with default configuration
./simc

# Run with custom configuration
./simc path/to/config.json
```

### Configuration File

Edit `data/config/default.json` to configure:
- Beam energy and properties
- Target material and geometry
- Spectrometer settings
- Monte Carlo parameters
- Reaction type and physics settings
- Output options

Example:
```json
{
  "beam": {
    "energy": 10.6,
    "energy_spread": 0.05
  },
  "target": {
    "type": "LH2",
    "thickness": 0.1
  },
  "reaction": {
    "type": "elastic"
  },
  "monte_carlo": {
    "nevents": 100000,
    "random_seed": 12345,
    "radiative_corrections": true,
    "energy_loss": true,
    "multiple_scattering": true
  }
}
```

### Supported Reactions

Phase 3 implements:
- **Elastic**: e + p → e' + p
- **Quasi-elastic**: e + A(p) → e' + p + (A-1) [In progress]
- **Pion production**: e + p → e' + π + X [In progress]
- **Kaon production**: e + p → e' + K + X [In progress]

### Output

SIMC produces a ROOT file (`simc_output.root` by default) containing:
- **T**: Main event tree with all kinematic quantities
- **T_gen**: Generated events (before cuts)
- **T_stats**: Run statistics
- **histograms/**: Diagnostic histograms

### Analyzing Output

```cpp
// Open ROOT file
TFile* f = new TFile("simc_output.root");
TTree* T = (TTree*)f->Get("T");

// Draw histogram
T->Draw("Q2/1e6 >> hQ2(100,0,5)", "success");

// Access event data
simc::SimcEvent* evt = nullptr;
T->SetBranchAddress("SimcEvent", &evt);
T->GetEntry(0);
evt->Print();

// Analyze kinematics
T->Draw("Em:Pm >> h(100,0,500,100,-100,200)", "success");
```

## Project Structure

```
simc_cpp/
├── CMakeLists.txt              # Main build configuration
├── README.md                   # This file
├── include/simc/               # Public headers
│   ├── SimcConstants.h         # Physical constants
│   ├── SimcTypes.h             # Common type definitions
│   ├── SimcEvent.h             # Event structure
│   ├── ConfigManager.h         # Configuration management
│   ├── RandomGenerator.h       # Random number generation
│   ├── OutputManager.h         # ROOT I/O
│   ├── Kinematics.h           # Kinematic calculations (NEW)
│   └── CrossSection.h         # Cross sections (NEW)
├── src/                        # Implementation
│   ├── core/                   # Core functionality
│   │   ├── SimcEvent.cpp
│   │   ├── ConfigManager.cpp
│   │   └── RandomGenerator.cpp
│   ├── physics/               # Physics calculations (NEW)
│   │   ├── Kinematics.cpp
│   │   └── CrossSection.cpp
│   ├── io/                     # Input/output
│   │   └── OutputManager.cpp
│   └── main.cpp                # Main program
├── data/config/                # Configuration files
│   └── default.json
├── tests/                      # Unit tests (standalone, no Google Test)
│   ├── test_main.cpp          # Test framework
│   ├── test_event.cpp
│   ├── test_config.cpp
│   └── test_random.cpp
└── docs/                       # Documentation
    └── phase2.tex              # LaTeX documentation
```

## Testing

### Running Tests

```bash
cd build

# Run all tests
./simc_tests

# Or use CTest
make test
```

**No external test framework required!** Tests use a simple standalone framework built into the project.

### Test Coverage

Current tests cover:
- ✅ Event data structures
- ✅ Configuration loading and parsing
- ✅ Random number generation
- 🚧 Kinematic calculations
- 🚧 Cross section calculations

## Physics Implementation

### Kinematic Calculations

The `Kinematics` class provides:
- Q², W, x_Bjorken, ε calculations
- Missing energy and momentum
- Angle transformations
- Coordinate system conversions

Example:
```cpp
using namespace simc;

SimcEvent evt;
Kinematics::Calculate(evt, 
    Ein=10600.0,    // Beam energy (MeV)
    Ee=8000.0,      // Electron energy (MeV)
    theta_e=0.218,  // Electron angle (rad)
    Ep=5000.0,      // Hadron energy (MeV)
    theta_p=0.523,  // Hadron angle (rad)
    phi_pq=0.0,     // Out-of-plane angle (rad)
    M_target=938.27, // Proton mass (MeV)
    M_hadron=938.27  // Proton mass (MeV)
);

std::cout << "Q² = " << evt.Q2/1e6 << " (GeV/c)²\n";
std::cout << "W = " << evt.W << " MeV\n";
std::cout << "x = " << evt.xbj << "\n";
```

### Cross Sections

The `CrossSection` framework provides:
- Elastic: Rosenbluth formula with form factors
- Quasi-elastic: Spectral function approach
- Pion/Kaon production: Hadronic models

Example:
```cpp
// Create cross section calculator
auto xs = CrossSectionFactory::Create(ReactionType::ELASTIC);

// Calculate for event
double sigma = xs->Calculate(evt);  // microbarns
```

## Development

### Code Style

- **C++17 standard** throughout
- **Header/source separation** for all classes
- **Doxygen-style** comments for all public APIs
- **snake_case** for variables and functions
- **PascalCase** for classes and types

### Adding New Physics

1. Create header in `include/simc/`
2. Add implementation to `src/physics/`
3. Update `src/CMakeLists.txt`
4. Add tests to `tests/`
5. Update documentation

Example for new reaction:
```cpp
// MyReaction.h
class MyReactionCrossSection : public CrossSectionBase {
public:
    double Calculate(const SimcEvent& evt) const override;
    // ... implement required methods
};
```

## Documentation

### Compile LaTeX Documentation

```bash
cd docs
pdflatex phase2.tex
pdflatex phase2.tex  # Run twice for references
```

This generates `phase2.pdf` with complete Phase 2 architecture documentation.

### Code Documentation

All public APIs are documented with Doxygen-style comments. To build HTML documentation:

```bash
cmake -DBUILD_DOCUMENTATION=ON ..
make doc
```

## Troubleshooting

### ROOT Not Found

```bash
# Set ROOT environment
source /path/to/root/bin/thisroot.sh

# Or specify ROOT location
cmake -DROOT_DIR=/path/to/root ..
```

### Compilation Errors

Make sure you have C++17 support:
```bash
g++ --version  # Should be ≥ 7.0
clang++ --version  # Should be ≥ 5.0
```

### Runtime Errors

Check ROOT file permissions:
```bash
ls -l simc_output.root
```

Verify configuration file format:
```bash
cat data/config/default.json | python -m json.tool
```

### Test Failures

Tests are standalone and don't require Google Test. If tests fail:

1. Check that core libraries compiled successfully
2. Verify ROOT is properly linked
3. Run with verbose output: `./simc_tests`

## Contributing

This is an active development project for Jefferson Lab. 

### Current Priorities (Phase 3):
1. Complete cross section implementations
2. Add event generators
3. Implement radiative corrections
4. Add multiple scattering
5. Energy loss calculations

### Future Phases:
- **Phase 4**: Spectrometer optics and apertures
- **Phase 5**: Complete detector simulation
- **Phase 6**: Validation and optimization

## Performance

Phase 3 targets:
- **Event generation**: >10,000 events/sec
- **Memory usage**: <100 MB for typical runs
- **ROOT I/O**: Efficient ntuple structure

## License

SIMC is developed at Jefferson Lab for nuclear physics research.

## Contact

For questions about SIMC C++/ROOT:
- Jefferson Lab Hall C: https://hallcweb.jlab.org/
- SIMC Documentation: https://hallcweb.jlab.org/wiki/index.php/SIMC_Monte_Carlo

## Acknowledgments

Based on the original SIMC Fortran code developed by the Jefferson Lab Hall C collaboration.

## Version History

- **2.0.0-phase3** (Current) - Physics implementation in progress
- **2.0.0-phase2** - Core infrastructure complete
- **1.x** - Original Fortran SIMC

---

**Version:** 2.0.0-phase3  
**Last Updated:** January 2026  
**Status:** Active Development
