/**
 * @file
 * @ingroup core
 * @brief Name-keyed factory registry of the plug-in stages, and the #SIMC_REGISTER macro.
 *
 * Every stage implementation registers a factory under a name for its interface; the driver creates stages by the
 * names given in the JSON `pipeline` section:
 * @code
 * SIMC_REGISTER(ISpectrometer, "hms", HmsSpectrometer);          // in HmsSpectrometer.cpp
 * auto s = simc::Registry<simc::ISpectrometer>::create("hms");   // in the driver
 * @endcode
 * Registration happens during static initialisation. Plug-ins must therefore be linked as object files (CMake OBJECT
 * libraries or whole-archive), otherwise the linker may drop a translation unit that nothing references.
 */
#ifndef SIMC_CORE_REGISTRY_H
#define SIMC_CORE_REGISTRY_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief Factory registry of the implementations of one stage interface.
 *
 * One registry exists per interface type (a function-local static map, so it is safe to use during static
 * initialisation in any order). Registration is meant to happen during static initialisation (#SIMC_REGISTER);
 * afterwards the registry is only read, and concurrent reads (create(), contains(), names()) are thread-safe.
 * Registering while other threads read is not.
 * @tparam Interface The stage interface (e.g. simc::IReaction); every factory returns a new instance of an
 *         implementation of it.
 */
template <class Interface>
class Registry {
public:
    /// A factory: returns a new, default-constructed (not yet configured) instance; the caller owns it.
    using Factory = std::function<std::unique_ptr<Interface>()>;

    /**
     * @brief Registers a factory under a name.
     * @param name The plug-in name used in the JSON pipeline (case-sensitive), e.g. `"hms"`.
     * @param f The factory.
     * @return Always true, so that the call can initialise a static variable (see #SIMC_REGISTER).
     * @throws std::logic_error if @p name is already registered for this interface.
     */
    static bool add(const std::string& name, Factory f) {
        auto& m = map();
        if (m.count(name)) throw std::logic_error("simc::Registry: duplicate plug-in name '" + name + "'");
        m.emplace(name, std::move(f));
        return true;
    }

    /**
     * @brief Creates a new instance of the plug-in registered as @p name.
     * @param name The plug-in name.
     * @return The new instance (owned by the caller), not yet configured or initialised.
     * @throws std::out_of_range if no plug-in has this name; the message lists the known names.
     */
    static std::unique_ptr<Interface> create(const std::string& name) {
        auto& m = map();
        auto it = m.find(name);
        if (it == m.end()) {
            std::string known;
            for (const auto& kv : m) known += (known.empty() ? "" : ", ") + kv.first;
            throw std::out_of_range("simc::Registry: no plug-in named '" + name + "' (known: " + known + ")");
        }
        return it->second();
    }

    /**
     * @brief Tells whether a plug-in name is registered.
     * @param name The plug-in name.
     * @return True if create(@p name) would succeed.
     */
    static bool contains(const std::string& name) { return map().count(name) != 0; }

    /**
     * @brief The registered plug-in names.
     * @return All names, in alphabetical order.
     */
    static std::vector<std::string> names() {
        std::vector<std::string> v;
        for (const auto& kv : map()) v.push_back(kv.first);
        return v;
    }

private:
    static std::map<std::string, Factory>& map() {
        static std::map<std::string, Factory> m;   // function-local: safe static-init order
        return m;
    }
};

/** @} */

}  // namespace simc

/**
 * @ingroup core
 * @brief Token pasting helper of #SIMC_REGISTER (implementation detail).
 * @param a First token.
 * @param b Second token.
 */
#define SIMC_REGISTRY_CONCAT_(a, b) a##b
/**
 * @ingroup core
 * @brief Pastes two tokens after macro expansion of both (implementation detail of #SIMC_REGISTER).
 * @param a First token (expanded, e.g. `__LINE__`).
 * @param b Second token (expanded).
 */
#define SIMC_REGISTRY_CONCAT(a, b) SIMC_REGISTRY_CONCAT_(a, b)

/**
 * @ingroup core
 * @brief Registers a stage implementation under a name, at static initialisation.
 *
 * Use it once at namespace scope in the implementation file of the plug-in. It defines a variable in an anonymous
 * namespace whose initialiser calls simc::Registry<Iface>::add(); a duplicate name throws std::logic_error during static
 * initialisation (the program terminates).
 * @param Iface The stage interface (e.g. `ISpectrometer`, resolved in the enclosing namespace).
 * @param Name The plug-in name (a string literal, e.g. `"hms"`).
 * @param Class The implementation: default-constructible and derived from @p Iface.
 */
#define SIMC_REGISTER(Iface, Name, Class)                                                      \
    namespace {                                                                                \
    const bool SIMC_REGISTRY_CONCAT(simc_registered_, __LINE__) =                              \
        ::simc::Registry<Iface>::add(Name, [] { return std::unique_ptr<Iface>(std::make_unique<Class>()); }); \
    }

#endif  // SIMC_CORE_REGISTRY_H
