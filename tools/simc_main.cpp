// tools/simc_main.cpp — the `simc++` command line (task P7.2). The implementation is simc::driver::CommandLineMain
// (src/driver/CommandLine.cpp, P15.3), so that programs with additional plug-ins can reuse it.
#include "simc/driver/CommandLine.h"

int main(int argc, char** argv) { return simc::driver::CommandLineMain(argc, argv); }
