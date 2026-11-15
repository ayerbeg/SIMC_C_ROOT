/**
 * @file
 * @ingroup core
 * @brief The SIMC++ exception hierarchy.
 *
 * | exception           | meaning |
 * |---------------------|---------|
 * | simc::Error         | base of every error that SIMC++ reports as a run error (a std::runtime_error) |
 * | simc::FortranStop   | the run reached a condition where the original program stops (a Fortran `STOP`); the message names the routine and carries the original text. The run cannot continue. |
 * | simc::NotPorted     | a feature of the original program that SIMC++ does not provide (reaction family, arm, option); the message names the feature |
 *
 * Programming errors (a stage used before `init()`, a missing collaborator) are std::logic_error; bad user input in
 * the JSON pipeline is std::invalid_argument. The `simc++` command line turns every std::exception into exit status 3.
 * Task P7.4: one hierarchy instead of one exception type per module.
 */
#ifndef SIMC_CORE_ERRORS_H
#define SIMC_CORE_ERRORS_H

#include <stdexcept>
#include <string>

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief Base class of every run error reported by SIMC++.
 *
 * Catch `simc::Error` to handle all SIMC++ run errors (and only those); the message (what()) is meant for the user.
 * Constructed with the std::runtime_error constructors (a message string).
 */
class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief The run reached a condition where the original Fortran SIMC executes `STOP`.
 *
 * Thrown by the set-up and the kernels when the input or the event leads to a state the original program refuses
 * (e.g. an unknown arm phi, an inconsistent target). The message names the routine (file:line) and repeats the
 * original text; the run cannot continue.
 */
class FortranStop : public Error {
public:
    using Error::Error;
};

/**
 * @brief A feature of the original program that SIMC++ does not provide.
 *
 * Thrown when an input selects a reaction family, spectrometer or option that is not available; the message names the
 * feature.
 */
class NotPorted : public Error {
public:
    using Error::Error;
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_ERRORS_H
