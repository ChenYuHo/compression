#include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <vector>
using namespace std;
using namespace std::chrono;

float get_prob(float input) {
    static const uint32_t mask1 = 0x3F800000;
    static const uint32_t mask2 = 0x3FFFFFFF;
    thread_local static union ieee {
        float f;
        uint32_t i;
    } num;
    num.f = input;
    num.i |= mask1;
    num.i &= mask2;
    return num.f;
}

uint8_t compress_one_pure(const float value) {
    int8_t exp = ilogbf(value);
    uint8_t compressed;
    //if (1.5 < get_prob(value)) exp+=1; // add 1
    if (rand() / float(RAND_MAX) + 1.0 < get_prob(value)) exp+=1;
    if (exp >= 17) compressed = 127;
    else if (exp <= -110 || value == 0) compressed = 0;
    else compressed = uint8_t(exp + 110);
    if (value < 0) compressed |= 0x80u;
    return compressed;
}

void compress(float *input, uint8_t *output, int num_elements) {
    for (unsigned i = 0; i < num_elements; ++i) {
        output[i] = compress_one_pure(input[i]);
    }
}

float decompress_one_pure(const uint8_t exponent) {
    if (exponent == 0) {
        return 0.;
    }
    thread_local static union ieee {
        float f;
        uint32_t i;
    } num;
    auto sign = (exponent & 0x80u);
    num.i = uint32_t((exponent & 0x7Fu) + 17) << 23u;
    return sign ? -num.f : num.f;
}

void decompress(uint8_t *input, float *output, int num_elements) {
    for (unsigned i = 0; i < num_elements; ++i) {
        output[i] = decompress_one_pure(input[i]);
    }
}


int main(int argc, char *argv[]) {
    srand(time(NULL));
    int num_elements = atoi(argv[1]);
    int count = atoi(argv[2]);
    vector<float> data;
    data.reserve(num_elements);
    vector<uint8_t> compressed(num_elements);
    vector<float> decompressed(num_elements);
    for(unsigned i = 0; i < num_elements; ++i) {
        data.push_back(float(i));
    }
    //for (const auto &v : data) cout << v << " ";
    //cout << endl;
    for (unsigned i = 1; i <= count; ++i) {
        high_resolution_clock::time_point begin = high_resolution_clock::now();
        compress(data.data(), compressed.data(), num_elements);
        high_resolution_clock::time_point end = high_resolution_clock::now();
        auto ns_taken = duration_cast<nanoseconds>(end - begin).count();
        cout << "compress round " << i << ": " << ns_taken << " nanoseconds elapsed, "
             << double(num_elements) / ns_taken * 1e9 << " element/s\n";
        //for (const auto &v : compressed) cout << (unsigned) v << " ";
        //cout << endl;
    }
    for (unsigned i = 1; i <= count; ++i) {
        high_resolution_clock::time_point begin_decompress = high_resolution_clock::now();
        decompress(compressed.data(), decompressed.data(), num_elements);
        high_resolution_clock::time_point end_decompress = high_resolution_clock::now();
        auto ns_taken_decompress = duration_cast<nanoseconds>(end_decompress - begin_decompress).count();
        cout << "decompress round " << i << ": " << ns_taken_decompress << " nanoseconds elapsed, "
             << double(num_elements) / ns_taken_decompress * 1e9 << " element/s\n";
        //for (const auto &v : decompressed) cout << v << " ";
        //cout << endl;
    }
//    int r = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.
    return 0;
}






 //   void print() override {
 //       cout << "original data:" << endl;
 //       for (const auto &v : original_data) cout << v << " ";
 //       cout << endl << "compressed data:" << endl;
 //       for (const auto &v : compressed_data) cout << (unsigned) v << " ";
 //       cout << endl << "decompressed data:" << endl;
 //       for (const auto &v : decompressed_data) cout << v << " ";
 //       cout << endl;

 //   }

