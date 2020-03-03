#include <iostream>
#include <cmath>
#include <random>

#define NUM_ELEMENTS 100

static inline uint8_t sw_compress_one(float input) {
    uint8_t result;
    static const uint32_t mask1 = 0x3F800000;
    static const uint32_t mask2 = 0x3FFFFFFF;
    static union ieee {
        float f;
        uint32_t i;
    } num;
    num.f = input;
    num.i |= mask1;
    num.i &= mask2;
    int8_t exp = ilogbf(input);
    if (1.5 < num.f) exp += 1; // add 1
    if (exp >= 17) result = 127;
    else if (exp <= -110 || input == 0) result = 0;
    else result = uint8_t(exp + 110);
    if (input < 0) result |= 0x80u;
    return result;
}

static inline float sw_decompress_one(uint8_t input) {
    if (input == 0) {
        return 0.;
    }
    static union ieee {
        float f;
        uint32_t i;
    } num;
    uint8_t sign = (input & 0x80u);
    num.i = uint32_t((input & 0x7Fu) + 17) << 23u;
    return sign ? -num.f : num.f;
}

int main(void) {
    float testinputs[NUM_ELEMENTS];
    uint8_t swresults[NUM_ELEMENTS];
//    ap_int<8> hwcomp[NUM_ELEMENTS];

    // random number generator
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<float> dist(0., 5.);

    unsigned error_count = 0;

    // Generate test vectors and expected results
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        // also test 0
        if (i == 0) testinputs[i] = 0;
        else testinputs[i] = dist(e2);
        swresults[i] = sw_compress_one(testinputs[i]);
    }

    // Run the accelerator

    // Check accelerator outputs
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        auto decompressed = sw_decompress_one(swresults[i]);

        // compute correct answer
        int exp;
        float frac = frexpf(testinputs[i], &exp);
        if (frac > 0) {
            if (frac > 0.75) frac = 1;
            else frac = 0.5;
        } else if (frac < 0) {
            if (frac < -0.75) frac = -1;
            else frac = -0.5;
        }
        float answer = ldexpf(frac, exp);
        if (decompressed != answer) {
            // printf("original %f, should get %f, but got %f\n", testinputs[i], answer, decompressed);
            error_count += 1;
        }
    }

    // Print final test status
    if (error_count) {
        std::cout << "!!! TEST FAILED - " << error_count;
        std::cout << " results incorrect." << std::endl;
    } else
        std::cout << "Test Passed" << std::endl;

    return error_count;
}