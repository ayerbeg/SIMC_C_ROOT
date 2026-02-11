// tests/test_event.cpp
// Unit tests for SimcEvent (standalone version)

#include "simc/core/SimcEvent.h"
#include "simc/core/SimcConstants.h"
#include <iostream>
#include <cmath>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_EQ(double a, double b, const std::string& msg = "");
template<typename T> void ASSERT_EQ_T(T a, T b, const std::string& msg = "");

using namespace simc;

void TestEvent() {
    SimcEvent evt;
    
    // Test 1: Default initialization
    evt.Clear();
    ASSERT_EQ(evt.Ein, 0.0, "Ein should be 0");
    ASSERT_EQ(evt.Q2, 0.0, "Q2 should be 0");
    ASSERT_EQ(evt.W, 0.0, "W should be 0");
    
    // Test 2: Clear method
    evt.Ein = 10000.0;
    evt.Q2 = 5000.0;
    evt.e_delta = 5.0;
    evt.Clear();
    ASSERT_EQ(evt.Ein, 0.0, "Ein should be cleared");
    ASSERT_EQ(evt.Q2, 0.0, "Q2 should be cleared");
    
    // Test 3: Set/Get electron state
    ArmState arm;
    arm.delta = 5.0;
    arm.xptar = 0.01;
    arm.yptar = 0.02;
    arm.E = 8000.0;
    arm.P = 8000.0;
    
    evt.SetElectron(arm);
    ASSERT_EQ(evt.e_delta, 5.0, "e_delta not set correctly");
    ASSERT_EQ(evt.e_xptar, 0.01, "e_xptar not set correctly");
    ASSERT_EQ(evt.e_E, 8000.0, "e_E not set correctly");
    
    ArmState retrieved = evt.GetElectronState();
    ASSERT_EQ(retrieved.delta, 5.0, "retrieved delta wrong");
    
    // Test 4: Unit vectors
    Vector3D ue{1.0, 0.0, 0.0};
    Vector3D up{0.0, 1.0, 0.0};
    Vector3D uq{0.0, 0.0, 1.0};
    
    evt.SetUnitVectors(ue, up, uq);
    ASSERT_EQ(evt.ue_x, 1.0, "ue_x wrong");
    ASSERT_EQ(evt.up_y, 1.0, "up_y wrong");
    ASSERT_EQ(evt.uq_z, 1.0, "uq_z wrong");
    
    std::cout << "  All event tests passed\n";
}
