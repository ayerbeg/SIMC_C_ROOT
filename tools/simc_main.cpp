// tools/simc_main.cpp — the `simc++` program. The implementation is simc::driver::CommandLineMain
// (src/driver/CommandLine.cpp), so that programs linked with additional plug-ins can reuse it.
#include "simc/driver/CommandLine.h"

int main(int argc, char** argv) { return simc::driver::CommandLineMain(argc, argv); }
