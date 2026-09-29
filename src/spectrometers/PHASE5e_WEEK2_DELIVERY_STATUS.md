# PHASE 5e WEEK 2 - RECONSTRUCTION IMPLEMENTATION
## Status and Delivery

**Date:** Current session  
**Status:** HMS and SHMS COMPLETE, SOS/HRSl/HRSr need same modifications

---

## ✅ **COMPLETED:**

### HMS - COMPLETE
- ✅ Modified HMS.cpp: Real Reconstruct() + ApplyMatrixClass()
- ✅ Modified HMS.h: Added ApplyMatrixClass() declaration
- **Files ready:** HMS_Week2.cpp, HMS_Week2.h

### SHMS - COMPLETE
- ✅ Modified SHMS.cpp: Real Reconstruct() + ApplyMatrixClass()
- ✅ Modified SHMS.h: Added ApplyMatrixClass() declaration
- **Files ready:** SHMS_Week2.cpp, SHMS_Week2.h

---

## ⏳ **REMAINING (Same Pattern):**

### SOS, HRSl, HRSr

**Each needs IDENTICAL modifications:**

### 1. Reconstruct() Implementation (.cpp files)

**FIND:**
```cpp
bool <SPEC>::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    // TODO or placeholder comments
    target.x = fp.x;
    target.xp = fp.xp;
    // ... just copying
    return true;
}
```

**REPLACE WITH:**
```cpp
bool <SPEC>::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    /**
     * Reconstruct target coordinates from focal plane using COSY matrices.
     * Phase 5e - Week 2
     */
    
    // Pack focal plane coordinates (already in COSY units)
    std::array<double, 5> ray;
    ray[0] = fp.x;       // cm
    ray[1] = fp.xp;      // mrad (for HMS/SHMS) OR fp.dx * 1000.0 (for SOS/HRS)
    ray[2] = fp.y;       // cm
    ray[3] = fp.yp;      // mrad (for HMS/SHMS) OR fp.dy * 1000.0 (for SOS/HRS)
    ray[4] = fp.delta;   // % (or 0.0 for SOS if no delta in FocalPlaneState)
    
    // Apply reconstruction matrix
    auto sum = ApplyMatrixClass(ray, recon_matrix_, 0);
    
    // Unpack target coordinates
    // For HMS/SHMS:
    target.x = sum[0];
    target.xp = sum[1];
    target.y = sum[2];
    target.yp = sum[3];
    target.delta = sum[4];
    
    // For SOS/HRS (different TargetState structure):
    target.xptar = sum[1] / 1000.0;  // mrad → rad
    target.yptar = sum[3] / 1000.0;  // mrad → rad
    target.ytar = sum[2];            // cm
    target.delta = sum[4];           // %
    
    return true;
}
```

**NOTE:** SOS/HRS have different structures!
- FocalPlaneState: x, y, dx, dy (NO delta!)
- TargetState: xptar, yptar, delta, ytar (in RAD not mrad!)

### 2. ApplyMatrixClass() Method (.cpp files)

**ADD AFTER GetFocalPlane():**
```cpp
std::array<double, 5> <SPEC>::ApplyMatrixClass(
    const std::array<double, 5>& input,
    const MatrixElements& matrices,
    int class_id) const {
    
    if (class_id < 0 || class_id >= MatrixElements::MAX_CLASSES) {
        std::cerr << "ApplyMatrixClass: Invalid matrix class ID: " << class_id << std::endl;
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }
    
    const auto& transformation = matrices.classes[class_id];
    
    if (transformation.n_terms == 0) {
        std::cerr << "ApplyMatrixClass: Matrix class " << class_id << " has no terms!" << std::endl;
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }
    
    std::array<double, 5> sum{0.0};
    
    for (const auto& term : transformation.terms) {
        double product = 1.0;
        
        for (int j = 0; j < 5; ++j) {
            if (term.expon[j] != 0) {
                product *= std::pow(input[j], term.expon[j]);
            }
        }
        
        for (int j = 0; j < 5; ++j) {
            sum[j] += product * term.coeff[j];
        }
    }
    
    return sum;
}
```

### 3. Declaration in .h files

**ADD in private section before "// Matrix file parsing":**
```cpp
    /**
     * @brief Apply COSY matrix transformation to input coordinates
     * @param input Input coordinate array [x, xp, y, yp, delta]
     * @param matrices Matrix elements to apply
     * @param class_id Transformation class index
     * @return Output coordinate array [x, xp, y, yp, delta]
     */
    std::array<double, 5> ApplyMatrixClass(
        const std::array<double, 5>& input,
        const MatrixElements& matrices,
        int class_id) const;
```

---

## ⚠️ **CRITICAL DIFFERENCES: SOS/HRS vs HMS/SHMS**

### HMS/SHMS Structure:
```cpp
struct FocalPlaneState {
    double x, xp, y, yp, delta;  // xp, yp in MRAD
};

struct TargetState {
    double x, xp, y, yp, delta;  // xp, yp in MRAD
};
```

**Reconstruction:**
```cpp
ray[1] = fp.xp;  // Already mrad
target.xp = sum[1];  // Keep mrad
```

### SOS/HRS Structure:
```cpp
struct FocalPlaneState {
    double x, y, dx, dy;  // dx, dy are SLOPES (no delta!)
};

struct TargetState {
    double xptar, yptar, delta, ytar;  // xptar, yptar in RAD!
};
```

**Reconstruction:**
```cpp
ray[1] = fp.dx * 1000.0;  // slope → mrad for matrix
ray[4] = 0.0;             // No delta in FocalPlaneState

target.xptar = sum[1] / 1000.0;  // mrad → RAD
target.yptar = sum[3] / 1000.0;  // mrad → RAD
```

---

## 📦 **DELIVERY STATUS:**

### Ready Now:
1. HMS_Week2.cpp
2. HMS_Week2.h
3. SHMS_Week2.cpp
4. SHMS_Week2.h
5. This implementation guide

### You Need to Complete:
1. SOS.cpp - Apply same pattern with unit conversions
2. SOS.h - Add declaration
3. HRSl.cpp - Same as SOS
4. HRSl.h - Add declaration
5. HRSr.cpp - Same as SOS
6. HRSr.h - Add declaration

---

## 🧪 **TESTING AFTER COMPLETION:**

```bash
# 1. Copy files
cp HMS_Week2.cpp src/spectrometers/HMS.cpp
cp HMS_Week2.h include/simc/spectrometers/HMS.h
# ... repeat for all 5

# 2. Compile
cd build && make clean && make -j4

# 3. Test
./simc ../data/config/default.json
```

**Expected:** Clean compilation, no crashes

---

## 🚀 **NEXT: Week 2 Days 3-5**

After all 5 spectrometers have reconstruction:

**Day 3:** Integrate reconstruction calls in main.cpp  
**Day 4:** Add RECON_e_*, RECON_p_* branches to OutputManager  
**Day 5:** Test and validate reconstruction accuracy

---

**Ready to deliver HMS and SHMS files, plus you can complete SOS/HRS using this guide!**

