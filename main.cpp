#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <immintrin.h>
#include <cmath>
#include <random>
#include <chrono>
#include "include/simdxorshift128plus.h"
#include "include/aligned_allocator.h"
#include "include/cxxopts.h"

using namespace std;
using namespace std::chrono;
using namespace cxxopts;

ParseResult parse(int argc, char *argv[]) {
    Options options(argv[0], "compression");
    options.add_options()
#ifdef DEBUG
            ("s,size", "number of float32 elements",
             value<uint32_t>()->default_value("100"))
#else
            ("s,size", "number of float32 elements",
             value<uint32_t>()->default_value("26214400")) // default: 100 MB
#endif
            ("i,input", "input file name, randomly generate elements if not provided", value<string>())
            ("o,data-output", "output file name to save data (no-op if given input file)", value<string>())
            ("r,result-output", "output file name to save compressed result", value<string>())
            ("m,method", "compress method", value<string>()->default_value("intml"))
#ifdef DEBUG
            ("repeat", "repeat ", value<uint32_t>()->default_value("1"))
#else
            ("repeat", "repeat ", value<uint32_t>()->default_value("10"))
#endif
            ("c,compress", "do and measure compression")
            ("d,decompress", "do and measure decompression")
            ("h,help", "Print help");
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
    vector<T, aligned_allocator<T, 64>> data;
public:
    Compressor() = default;

    explicit Compressor(uint32_t num_elements) {
        this->num_elements = num_elements;
        data.assign(num_elements, 0);
    };

    void read_data(const string &filename) {
        ifstream input(filename);
        if (!input.good()) {
            cerr << filename << " does not exist" << endl;
            exit(1);
        }
        copy(istream_iterator<T>(input), istream_iterator<T>(), data.data());
        this->num_elements = data.size();
    }

    virtual void generate_data() = 0;

    virtual void compress() = 0;

    virtual void decompress() = 0;

    void write_data(const string &filename) {
        ofstream file(filename);
        for (const auto &v : this->data) file << v << "\n";
    }

    virtual void write_result(const string &) = 0;

};

class IntMLCompressor :
        public Compressor<float> {
protected:
    vector<uint8_t, aligned_allocator<uint8_t, 64>> result;
public:

    IntMLCompressor() = default;

    explicit IntMLCompressor(uint32_t num_elements) : Compressor(num_elements) {
        result.assign(num_elements, 0);
    };

    void generate_data() override {
        random_device rd;
        mt19937 e2(rd());
        normal_distribution<float> dist(0., 0.001);
        auto gen = [&dist, &e2]() {
            return dist(e2);
        };
        generate(this->data.begin(), this->data.end(), gen);
    }

    void write_result(const string &filename) override {
        ofstream file(filename);
        for (const auto &v : this->result) file << (unsigned) v << "\n";
    }

    void compress() override {
#ifdef DEBUG
        cout<<"before compression:"<<endl;
        for (const float &v : this->data) cout << v << " ";
        cout<<endl;
#endif
        float *current_data_ptr = this->data.data();
        uint8_t *current_result_ptr = this->result.data();

#pragma omp parallel for default(none) shared(current_data_ptr, current_result_ptr)
        for (unsigned i = 0; i < num_elements / 32; ++i) {
            fill_avx(current_data_ptr + i * 32, current_result_ptr + i * 32);
        }

        unsigned offset = 32 * (this->num_elements / 32);
        current_data_ptr += offset;
        current_result_ptr += offset;

#pragma omp parallel for default(none) shared(current_data_ptr, current_result_ptr)
        for (unsigned i = 0; i < num_elements % 32; ++i) {
            fill_one(current_data_ptr + i, current_result_ptr + i);
        }
#ifdef DEBUG
        cout<<"after compression:"<<endl;
        for (const unsigned &v : this->result) cout << v << " ";
        cout<<endl;
#endif
    }

private:

    static inline float get_prob(float input) {
        static const uint32_t mask1 = 0x3F800000;
        static const uint32_t mask2 = 0x3FFFFFFF;
        static union ieee {
            float f;
            uint32_t i;
        } num;
        num.f = input;
        num.i |= mask1;
        num.i &= mask2;
        return num.f;
    }

    static inline void fill_one(const float *value, uint8_t *compressed) {
        thread_local static random_device rd;
        thread_local static mt19937 e2(rd());
        thread_local static uniform_real_distribution<float> dist(1, 2);
        if (*value == 0) {
            *compressed = 0x40;
            return;
        }
        int8_t exp = ilogbf(*value);
#ifndef NORANDOM
        if (dist(e2) < get_prob(*value)) exp += 1; // add 1
#else
        if (1.5 < get_prob(*value)) exp+=1; // add 1
#endif
        if (exp >= 10) *compressed = 60;
        else if (exp <= -50) *compressed = 0;
        else *compressed = uint8_t(exp + 50);
        if (*value < 0) *compressed |= 0x80u;
    }

    static inline void fill_avx(const float *p, uint8_t *exp) {
        thread_local static avx_xorshift128plus_key_t key;
        // load 32 floats
        __m256i a = _mm256_castps_si256(_mm256_loadu_ps(p));
        __m256i b = _mm256_castps_si256(_mm256_loadu_ps(p + 8));
        __m256i c = _mm256_castps_si256(_mm256_loadu_ps(p + 16));
        __m256i d = _mm256_castps_si256(_mm256_loadu_ps(p + 24));

        // mask 23bits
        const static __m256i m = _mm256_set1_epi32(0x007fffff);

        // mask out sign and exponent, only leave mantissa
        __m256i am = _mm256_and_si256(a, m);
        __m256i bm = _mm256_and_si256(b, m);
        __m256i cm = _mm256_and_si256(c, m);
        __m256i dm = _mm256_and_si256(d, m);

#ifndef NORANDOM
        // if am larger than random numbers
        __m256i aa = _mm256_cmpgt_epi32(am, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i bb = _mm256_cmpgt_epi32(bm, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i cc = _mm256_cmpgt_epi32(cm, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i dd = _mm256_cmpgt_epi32(dm, _mm256_and_si256(avx_xorshift128plus(&key), m));
#else
        // rounding
        const static __m256i round = _mm256_set1_epi32(0x00400000);
        // comparing only mantissa doesn't raise sign issues
        __m256i aa = _mm256_cmpgt_epi32(am, round);
        __m256i bb = _mm256_cmpgt_epi32(bm, round);
        __m256i cc = _mm256_cmpgt_epi32(cm, round);
        __m256i dd = _mm256_cmpgt_epi32(dm, round);
#endif

        // shift right (zero extending), only sign and exponent left, add one if am larger than ar (undefined behavior if float is Inf)
        const static __m256i one = _mm256_set1_epi32(1);
        __m256i aexp = _mm256_add_epi32(_mm256_srli_epi32(a, 23), _mm256_and_si256(aa, one));
        __m256i bexp = _mm256_add_epi32(_mm256_srli_epi32(b, 23), _mm256_and_si256(bb, one));
        __m256i cexp = _mm256_add_epi32(_mm256_srli_epi32(c, 23), _mm256_and_si256(cc, one));
        __m256i dexp = _mm256_add_epi32(_mm256_srli_epi32(d, 23), _mm256_and_si256(dd, one));

        // pack as 16 uint16_t
        __m256i ab = _mm256_packus_epi32(aexp, bexp);
        __m256i cd = _mm256_packus_epi32(cexp, dexp);

        // save sign and put to the 8th bit
        __m256i absign = _mm256_slli_epi16(_mm256_srli_epi16(ab, 8), 7);
        __m256i cdsign = _mm256_slli_epi16(_mm256_srli_epi16(cd, 8), 7);
        __m256i abcdsign = _mm256_packus_epi16(absign, cdsign);

        // mask out sign
        const static __m256i m2 = _mm256_set1_epi16(0x00FF);
        __m256i abmasksign = _mm256_and_si256(ab, m2);
        __m256i cdmasksign = _mm256_and_si256(cd, m2);

        // pack as 32 uint8_t
        __m256i abcd = _mm256_packus_epi16(abmasksign, cdmasksign);

        // check 0 for later use
        __m256i abcdiszero = _mm256_cmpeq_epi8(abcd, _mm256_setzero_si256());

        // -50 + 127, smallest to represent 2^-50
        const static __m256i seventyseven = _mm256_set1_epi8(77);
        // biggest to represent 2^10 is 60 (10+127-77), so to saturate, 255-60=195
        const static __m256i onenintyfive = _mm256_set1_epi8(195);

        __m256i abcdshift = _mm256_subs_epu8(abcd, seventyseven); // saturate lower
        __m256i abcdsaturate = _mm256_subs_epu8(_mm256_adds_epu8(abcdshift, onenintyfive),
                                                onenintyfive); // saturate upper
        __m256i abcdsigned = _mm256_or_si256(abcdsign, abcdsaturate); // append sign

        // set 7th bit if value is zero
        const static __m256i seventhbit = _mm256_set1_epi8(0x40);
        __m256i abcdsigned_with_zero = _mm256_or_si256(abcdsigned, _mm256_and_si256(abcdiszero, seventhbit));

        // packing will result in different order, but we unpack when storing so no need to shuffle
        // const static __m256i shuffleidx = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
        // __m256i abcdordered = _mm256_permutevar8x32_epi32(abcdsigned_with_zero, shuffleidx);

        _mm256_storeu_si256((__m256i *) exp, abcdsigned_with_zero);
    }

    void decompress() override {}
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
        high_resolution_clock::time_point begin = high_resolution_clock::now();
        for (unsigned i = 0; i < count; ++i)
            compressor.compress();
        high_resolution_clock::time_point end = high_resolution_clock::now();
        auto us_taken = duration_cast<microseconds>(end - begin).count() / count;
        cout << "compress: " << us_taken << " microseconds "
             << num_elements * 1e6 * count / us_taken  << " CTE/s" << endl;
    }
    if (result.count("decompress")) {
        high_resolution_clock::time_point begin = high_resolution_clock::now();
        for (unsigned i = 0; i < count; ++i)
            compressor.decompress();
        high_resolution_clock::time_point end = high_resolution_clock::now();
        auto us_taken = duration_cast<microseconds>(end - begin).count() / count;
        cout << "decompress: " << us_taken << " microseconds "
             << num_elements * 1e6 * count / us_taken << " CTE/s" << endl;
    }
    if (!result.count("input") && result.count("output")) {
        compressor.write_data(result["output"].as<string>());
    }
    if (result.count("compressed-output")) {
        compressor.write_result(result["compressed-output"].as<string>());
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