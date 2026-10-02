// include/simc/core/Registry.h
//
// Name-keyed factory registry for plug-in stages (docs/ARCHITECTURE.md §3).
//
//   SIMC_REGISTER(ISpectrometer, "hms", HMS);                 // in HMS.cpp
//   auto s = Registry<ISpectrometer>::create("hms");          // in the driver
//
// Registration happens during static initialisation. Plug-ins must therefore be linked
// as object files (CMake OBJECT libraries or whole-archive), otherwise the linker may drop
// a translation unit that nothing references.
#ifndef SIMC_CORE_REGISTRY_H
#define SIMC_CORE_REGISTRY_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace simc {

template <class Interface>
class Registry {
public:
    using Factory = std::function<std::unique_ptr<Interface>()>;

    /// Register a factory under `name`. Returns true (usable in static initialisers).
    /// Registering the same name twice is an error (throws).
    static bool add(const std::string& name, Factory f) {
        auto& m = map();
        if (m.count(name)) throw std::logic_error("simc::Registry: duplicate plug-in name '" + name + "'");
        m.emplace(name, std::move(f));
        return true;
    }

    /// Create the plug-in registered as `name`; throws std::out_of_range listing known names.
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

    static bool contains(const std::string& name) { return map().count(name) != 0; }

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

}  // namespace simc

#define SIMC_REGISTRY_CONCAT_(a, b) a##b
#define SIMC_REGISTRY_CONCAT(a, b) SIMC_REGISTRY_CONCAT_(a, b)

/// Register `Class` (default-constructible, derived from `Iface`) under the string `Name`.
#define SIMC_REGISTER(Iface, Name, Class)                                                      \
    namespace {                                                                                \
    const bool SIMC_REGISTRY_CONCAT(simc_registered_, __LINE__) =                              \
        ::simc::Registry<Iface>::add(Name, [] { return std::unique_ptr<Iface>(std::make_unique<Class>()); }); \
    }

#endif  // SIMC_CORE_REGISTRY_H
