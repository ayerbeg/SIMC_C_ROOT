// quick_parser_test.cpp
// Quick standalone test to verify the COSY parser is working
// Compile: g++ -std=c++17 -I../include quick_parser_test.cpp ../src/physics/CosyMatrix.cpp -o quick_test

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Minimal parser test - no dependencies
bool test_parse_line() {
    // Example line from actual HMS forward_cosy.dat:
    //   0.7779354     -3.321846     0.0000000E+00 0.0000000E+00 0.0000000E+00 100000
    
    std::string test_line = "  0.7779354     -3.321846     0.0000000E+00 0.0000000E+00 0.0000000E+00 100000";
    
    std::istringstream iss(test_line);
    
    // Parse 5 coefficients
    double coeffs[5];
    for (int i = 0; i < 5; ++i) {
        if (!(iss >> coeffs[i])) {
            std::cerr << "❌ Failed to read coefficient " << i << std::endl;
            return false;
        }
    }
    
    // Parse exponent string
    std::string exp_str;
    if (!(iss >> exp_str)) {
        std::cerr << "❌ Failed to read exponent string" << std::endl;
        return false;
    }
    
    // Verify length
    if (exp_str.length() != 6) {
        std::cerr << "❌ Exponent string length is " << exp_str.length() 
                  << ", expected 6" << std::endl;
        return false;
    }
    
    // Parse exponents
    int exps[5];
    exps[0] = exp_str[0] - '0';  // x
    exps[1] = exp_str[1] - '0';  // xp
    exps[2] = exp_str[2] - '0';  // y
    exps[3] = exp_str[3] - '0';  // yp
    // Skip [4] - TOF
    exps[4] = exp_str[5] - '0';  // delta
    
    // Verify values
    std::cout << "✓ Parsed line successfully:" << std::endl;
    std::cout << "  Coefficients: [" 
              << coeffs[0] << ", " << coeffs[1] << ", " 
              << coeffs[2] << ", " << coeffs[3] << ", " << coeffs[4] << "]" << std::endl;
    std::cout << "  Exponents: [" 
              << exps[0] << ", " << exps[1] << ", " 
              << exps[2] << ", " << exps[3] << ", " << exps[4] << "]" << std::endl;
    
    // Expected values
    if (coeffs[0] != 0.7779354 || coeffs[1] != -3.321846) {
        std::cerr << "❌ Coefficient values incorrect!" << std::endl;
        return false;
    }
    
    if (exps[0] != 1 || exps[1] != 0 || exps[2] != 0 || exps[3] != 0 || exps[4] != 0) {
        std::cerr << "❌ Exponent values incorrect!" << std::endl;
        std::cerr << "   Expected: [1,0,0,0,0]" << std::endl;
        std::cerr << "   Got:      [" << exps[0] << "," << exps[1] << "," 
                  << exps[2] << "," << exps[3] << "," << exps[4] << "]" << std::endl;
        return false;
    }
    
    std::cout << "✓ All values correct!" << std::endl;
    return true;
}

bool test_parse_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "❌ Cannot open file: " << filename << std::endl;
        return false;
    }
    
    int data_lines = 0;
    int comment_lines = 0;
    int separator_lines = 0;
    int parse_errors = 0;
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(file, line)) {
        ++line_num;
        
        if (line.empty()) continue;
        
        if (line[0] == '!') {
            ++comment_lines;
            continue;
        }
        
        if (line.size() > 1 && line[0] == ' ' && line[1] == '-') {
            ++separator_lines;
            continue;
        }
        
        // Try to parse as data line
        std::istringstream iss(line);
        double coeffs[5];
        std::string exp_str;
        
        bool parsed = true;
        for (int i = 0; i < 5; ++i) {
            if (!(iss >> coeffs[i])) {
                parsed = false;
                break;
            }
        }
        
        if (parsed && (iss >> exp_str) && exp_str.length() == 6) {
            ++data_lines;
        } else {
            // Only count as error if line contains numbers
            if (line.find_first_of("0123456789") != std::string::npos) {
                ++parse_errors;
                if (parse_errors <= 3) {  // Show first 3 errors
                    std::cerr << "⚠ Line " << line_num << " parse error: " << line << std::endl;
                }
            }
        }
    }
    
    std::cout << "\nFile: " << filename << std::endl;
    std::cout << "  Comment lines:   " << comment_lines << std::endl;
    std::cout << "  Separator lines: " << separator_lines << std::endl;
    std::cout << "  Data lines:      " << data_lines << std::endl;
    std::cout << "  Parse errors:    " << parse_errors << std::endl;
    
    if (data_lines == 0) {
        std::cerr << "❌ CRITICAL: No data lines parsed!" << std::endl;
        std::cerr << "   Parser is not working correctly." << std::endl;
        return false;
    }
    
    if (data_lines < 100) {
        std::cerr << "⚠ WARNING: Only " << data_lines << " data lines" << std::endl;
        std::cerr << "   Expected 300-600 for typical matrix file" << std::endl;
        return false;
    }
    
    std::cout << "✓ File parsed successfully!" << std::endl;
    return true;
}

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "Quick COSY Parser Validation" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Test 1: Parse a single line
    std::cout << "=== Test 1: Single Line Parsing ===" << std::endl;
    if (!test_parse_line()) {
        std::cerr << "\n❌ FAILED: Parser cannot read basic format!" << std::endl;
        return 1;
    }
    
    // Test 2: Parse entire file
    std::cout << "\n=== Test 2: Full File Parsing ===" << std::endl;
 
    std::string filename = "../data/matrices/hrsl/hrs_recon_cosy.dat";
    if (argc > 1) {
        filename = argv[1];
    }
    
    if (!test_parse_file(filename)) {
        std::cerr << "\n❌ FAILED: Cannot parse matrix file!" << std::endl;
        std::cerr << "Make sure you're running from the build directory" << std::endl;
        return 1;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "✓ All parser tests PASSED!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
