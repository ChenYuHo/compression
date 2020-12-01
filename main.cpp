#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <cmath>
#include <random>
#include <chrono>
#include "include/cxxopts.h"

using namespace std;
using namespace std::chrono;
using namespace cxxopts;

ParseResult parse(int argc, char *argv[]) {
    Options options(argv[0], "compression");
    options.add_options()
            ("i,input", "input file name, randomly generate elements if not provided", value<string>())
            ("o,data-output", "output file name to save original_data (no-op if given input file)",
             value<string>())
            ("r,compressed-output", "output file name to save compressed compressed_data", value<string>())
            ("decompressed-output", "output file name to save decompressed compressed_data", value<string>())
            ("m,method", "compress method", value<string>()->default_value("intml"))
#ifdef DEBUG
    ("s,size", "number of float32 elements",
value<uint32_t>()->default_value("100"))
    ("repeat", "repeat compression and/or decompression", value<uint32_t>()->default_value("1"))
#else
            ("s,size", "number of float32 elements",
             value<uint32_t>()->default_value("26214400")) // default: 100 MB
            ("repeat", "repeat compression and/or decompression", value<uint32_t>()->default_value("10"))
#endif
            ("p,print", "print original, compressed, and decompressed data")
            ("c,compress", "do and measure compression")
            ("d,decompress", "do and measure decompression")
            ("h,help", "print usage");
    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cout << options.help() << endl;
            exit(0);
        }
        return result;
    } catch (const OptionException &e) {
        cout << "error parsing options: " << e.what() << endl;
        cout << options.help() << endl;
        exit(1);
    }
}

template<typename T>
class Compressor {
protected:
    uint32_t num_elements = 0;
    vector<T> original_data;
    vector<T> decompressed_data;
public:

    void read_data(const string &filename) {
        ifstream input(filename);
        if (!input.good()) {
            cerr << filename << " does not exist" << endl;
            exit(1);
        }
        for (istream_iterator<T> p{input}, e; p != e; ++p) {
            original_data.push_back(*p);
        }
        num_elements = original_data.size();
        init_compressed_data();
    }

    virtual void init_compressed_data() = 0;

    virtual void init_decompressed_data() = 0;

    virtual void generate_data() = 0;

    virtual void compress() = 0;

    virtual void decompress() = 0;

    virtual void print() = 0;

    void write_original_data(const string &filename) {
        ofstream file(filename);
        for (const T &v : original_data) file << v << " ";
        file << endl;
    }

    void write_decompressed_data(const string &filename) {
        ofstream file(filename);
        for (const T &v : decompressed_data) file << v << " ";
        file << endl;
    }

    virtual void write_compressed_data(const string &) = 0;

};

class IntMLCompressor :
        public Compressor<float> {
public:
    IntMLCompressor() : Compressor() {};

    explicit IntMLCompressor(uint32_t num_elements) : IntMLCompressor() {
        this->num_elements = num_elements;
    };

    void generate_data() override {
        random_device rd;
        mt19937 e2(rd());
        normal_distribution<float> dist(0., 0.001);
        auto gen = [&dist, &e2]() {
            return dist(e2);
        };
        original_data.assign(num_elements, 0);
        init_compressed_data();
        generate(original_data.begin(), original_data.end(), gen);
    }

    inline void init_compressed_data() override {
        compressed_data.assign(num_elements, 0);
    }

    inline void init_decompressed_data() override {
        decompressed_data.assign(num_elements, 0);
    }


    void write_compressed_data(const string &filename) override {
        ofstream file(filename);
        for (const auto &v : compressed_data) file << (unsigned) v << " ";
        file << endl;
    }

    void print() override {
        cout << "original data:" << endl;
        for (const auto &v : original_data) cout << v << " ";
        cout << endl << "compressed data:" << endl;
        for (const auto &v : compressed_data) cout << (unsigned) v << " ";
        cout << endl << "decompressed data:" << endl;
        for (const auto &v : decompressed_data) cout << v << " ";
        cout << endl;

    }

    void compress() override {
        float *current_data_ptr = original_data.data();
        uint8_t *current_result_ptr = compressed_data.data();

        for (unsigned i = 0; i < num_elements; ++i) {
            compress_one(current_data_ptr + i, current_result_ptr + i);
        }
    }

private:
    vector<uint8_t> compressed_data;

    static inline float get_prob(float input) {
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

    static inline void compress_one(const float *value, uint8_t *compressed) {
        thread_local static random_device rd;
        thread_local static mt19937 e2(rd());
        thread_local static uniform_real_distribution<float> dist(1, 2);
        int8_t exp = ilogbf(*value);
#ifndef NORANDOM
        if (dist(e2) < get_prob(*value)) exp += 1; // add 1
#else
        if (1.5 < get_prob(*value)) exp+=1; // add 1
#endif
        if (exp >= 17) *compressed = 127;
        else if (exp <= -110 || *value == 0) *compressed = 0;
        else *compressed = uint8_t(exp + 110);
        if (*value < 0) *compressed |= 0x80u;
    }

    void decompress() override {
        float *result_ptr = decompressed_data.data();
        uint8_t *exp_ptr = compressed_data.data();

        for (unsigned i = 0; i < num_elements; ++i) {
            decompress_one(result_ptr + i, exp_ptr + i);
        }
    }

    static inline void decompress_one(float *result, const uint8_t *exponent) {
        if (*exponent == 0) {
            *result = 0;
            return;
        }
        thread_local static union ieee {
            float f;
            uint32_t i;
        } num;
        auto sign = (*exponent & 0x80u);
        num.i = uint32_t((*exponent & 0x7Fu) + 17) << 23u;
        *result = sign ? -num.f : num.f;
    }
};

void do_work(const ParseResult &result, Compressor<float> &compressor) {
    if (result.count("input")) {
        compressor.read_data(result["input"].as<string>());
    } else {
        compressor.generate_data();
    }
    uint32_t count = result["repeat"].as<uint32_t>();
    uint32_t num_elements = result["size"].as<uint32_t>();
    if (result.count("compress")) {
        for (unsigned i = 1; i <= count; ++i) {
            compressor.init_compressed_data();
            high_resolution_clock::time_point begin = high_resolution_clock::now();
            compressor.compress();
            high_resolution_clock::time_point end = high_resolution_clock::now();
            auto ns_taken = duration_cast<nanoseconds>(end - begin).count();
            cout << "compress round " << i << ": " << ns_taken << " nanoseconds elapsed, "
                 << double(num_elements) / ns_taken * 1e9 << " element/s\n";
        }
    }
    if (result.count("decompress")) {
        for (unsigned i = 1; i <= count; ++i) {
            compressor.init_decompressed_data();
            high_resolution_clock::time_point begin = high_resolution_clock::now();
            compressor.decompress();
            high_resolution_clock::time_point end = high_resolution_clock::now();
            auto ns_taken = duration_cast<nanoseconds>(end - begin).count();
            cout << "decompress round " << i << ": " << ns_taken << " nanoseconds elapsed, "
                 << double(num_elements) / ns_taken * 1e9 << " element/s\n";
        }
    }
    if (result["print"].as<bool>() && !result.count("no-print")) {
        compressor.print();
    }
    if (!result.count("input") && result.count("data-output")) {
        compressor.write_original_data(result["data-output"].as<string>());
    }
    if (result.count("compressed-output")) {
        compressor.write_compressed_data(result["compressed-output"].as<string>());
    }
    if (result.count("decompressed-output")) {
        compressor.write_decompressed_data(result["decompressed-output"].as<string>());
    }
}

int main(int argc, char *argv[]) {
    auto result = parse(argc, argv);
    auto method = result["method"].as<string>();
    if (method == "intml") {
        IntMLCompressor compressor(result["size"].as<uint32_t>());
        do_work(result, compressor);
    } else {
        cerr << method << " not implemented" << endl;
        return 1;
    }
    return 0;
}
