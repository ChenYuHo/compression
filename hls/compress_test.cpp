#include <iostream>
#include <cmath>
#include <random.h>

#include "compress.hpp"

#define NUM_ELEMENTS 100

int main(void) {
	float testinputs[NUM_ELEMENTS];
	uint8_t swresults[NUM_ELEMENTS];
	ap_int<8> hwcomp[NUM_ELEMENTS];

	unsigned error_count = 0;

	// Generate test vectors and expected results
	for (int i = 0; i < NUM_ELEMENTS; i++) {

	}

	// Run the accelerator

	// Check accelerator outputs

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
