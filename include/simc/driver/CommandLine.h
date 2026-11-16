/**
 * @file
 * @ingroup driver
 * @brief The simc++ command line as a library function (task P15.3).
 *
 * the `simc++` command line as a library function (task P15.3).
 * tools/simc_main.cpp is `return simc::driver::CommandLineMain(argc, argv);`. A program that links extra plug-ins
 * (classes registered with SIMC_REGISTER, docs/manual/markdown/reference/EXTENDING.md) calls the same function and gets the full simc++ command
 * line with its plug-ins selectable from the JSON pipeline (examples/out_of_tree_plugin).
 * Usage, options and exit status: src/driver/CommandLine.cpp (or `simc++` without arguments).
 */
#ifndef SIMC_DRIVER_COMMANDLINE_H
#define SIMC_DRIVER_COMMANDLINE_H

namespace simc {
namespace driver {

/// @addtogroup driver
/// @{

/// `simc++ <file.inp|file.json> [-o OUTDIR] [--seed N] [--ngen N] [--pipeline FILE] ...`; returns the exit status.
/// @param argc argument count, as given to main().
/// @param argv argument vector, as given to main().
/// @return the process exit status.
int CommandLineMain(int argc, char** argv);

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_COMMANDLINE_H
