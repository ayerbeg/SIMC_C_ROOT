#include "simc/EventGenerator.h"
#include "simc/ConfigManager.h"
#include "simc/RandomGenerator.h"
#include <iostream>

int main() {
    try {
        simc::ConfigManager config("../data/config/default.json");
        auto random = std::make_shared<simc::RandomGenerator>(12345);
        simc::EventGenerator gen(config, random);
        
        if (gen.Initialize()) {
            std::cout << "EventGenerator initialized successfully!\n";
            std::cout << "Reaction type: " 
                      << (int)gen.GetReactionType() << "\n";
            return 0;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 1;
}
