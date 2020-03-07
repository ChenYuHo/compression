#include <iostream>
#include <cmath>
#include <random>

#include "compress.hpp"

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


int main(void) {
	float testinputs[NUM_ELEMENTS];
	uint8_t swresults[NUM_ELEMENTS];
	uint8_t hwcomp[NUM_ELEMENTS];

	hls::stream<apuint32_t> floatstream;
	hls::stream<uint8_t> compstream;

	std::random_device rd;
	std::mt19937 e2(rd());
	std::normal_distribution<float> dist(0., 5.);

	unsigned error_count = 0;

	// Generate test vectors and expected results
	for (int i = 0; i < NUM_ELEMENTS; i++) {
		// also test 0
		if (i == 0) testinputs[i] = 0;
		else testinputs[i] = dist(e2);
		printf("%f  ",testinputs[i]);
		swresults[i] = sw_compress_one(testinputs[i]);
		printf("%u\n",swresults[i]);

	}

	// Load test vectors into stream object
	for (int i = 0; i < NUM_ELEMENTS; i++) {
		float inval = testinputs[i];
		ap_uint<32> intVal = *(ap_uint<32>*) &inval;
		floatstream.write(intVal);
	}

	// Run the accelerator
	for (int i = 0; i < NUM_ELEMENTS; i++) {
		compress_one(floatstream, compstream);
	}

	// Unload results from stream object
	for (int i = 0; i < NUM_ELEMENTS; i++) {
			compstream.read(hwcomp[i]);
	}

	// Check accelerator outputs
	for (int i = 0; i < NUM_ELEMENTS; i++) {
		printf("SW:%d    HW:%d\n", swresults[i], (uint8_t) hwcomp[i]);
		if (swresults[i] != (uint8_t) hwcomp[i])
			error_count += 1;
	}


	// Print final test status
	if (error_count) {
		std::cout << "!!! TEST FAILED - " << error_count;
		std::cout << " results incorrect." << std::endl;
	} else
		std::cout << "Test Passed" << std::endl;

	return error_count;
}


//void generate_data() {
//	std::random_device rd;
//	std::mt19937 e2(rd());
//	std::normal_distribution<float> dist(0., 0.001);
//	auto gen = [&dist, &e2]() {
//		return dist(e2);
//    };
//	original_data.assign(NUM_ELEMENTS, 0);
//	init_compressed_data();
//	generate(original_data.begin(), original_data.end(), gen);
//}
