#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <immintrin.h>
#include <cmath>
#include <random>
#include <chrono>
#include "include/aligned_allocator.h"
#include "include/cxxopts.h"
#ifndef NORANDOM
#include "include/simdxorshift128plus.h"
#endif

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
    vector<T, aligned_allocator<T, 64>> original_data;
    vector<T, aligned_allocator<T, 64>> decompressed_data;
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

#pragma omp parallel for default(none) shared(current_data_ptr, current_result_ptr)
        for (unsigned i = 0; i < num_elements / 32; ++i) {
            compress_avx(current_data_ptr + i * 32, current_result_ptr + i * 32);
        }

        unsigned offset = 32 * (num_elements / 32);
        current_data_ptr += offset;
        current_result_ptr += offset;

#pragma omp parallel for default(none) shared(current_data_ptr, current_result_ptr)
        for (unsigned i = 0; i < num_elements % 32; ++i) {
            compress_one(current_data_ptr + i, current_result_ptr + i);
        }
    }

private:
    vector<uint8_t, aligned_allocator<uint8_t, 64>> compressed_data;

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

    static inline void compress_avx(const float *p, uint8_t *exp) {
#ifndef NORANDOM
        thread_local static random_device rd;
        thread_local static avx_xorshift128plus_key_t key;
        avx_xorshift128plus_init(rd(), rd(), &key);
#endif
        // load 32 floats
        __m256i a = _mm256_castps_si256(_mm256_loadu_ps(p));
        __m256i b = _mm256_castps_si256(_mm256_loadu_ps(p + 8));
        __m256i c = _mm256_castps_si256(_mm256_loadu_ps(p + 16));
        __m256i d = _mm256_castps_si256(_mm256_loadu_ps(p + 24));

        // mask out sign and exponent, only leave mantissa
        const static __m256i m = _mm256_set1_epi32(0x007fffff);
        __m256i am = _mm256_and_si256(a, m);
        __m256i bm = _mm256_and_si256(b, m);
        __m256i cm = _mm256_and_si256(c, m);
        __m256i dm = _mm256_and_si256(d, m);

#ifndef NORANDOM
        // stochastic rounding
        // true if mantissa larger than a random number (exp should +1)
        __m256i aa = _mm256_cmpgt_epi32(am, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i bb = _mm256_cmpgt_epi32(bm, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i cc = _mm256_cmpgt_epi32(cm, _mm256_and_si256(avx_xorshift128plus(&key), m));
        __m256i dd = _mm256_cmpgt_epi32(dm, _mm256_and_si256(avx_xorshift128plus(&key), m));
#else
        // deterministic rounding
        const static __m256i round = _mm256_set1_epi32(0x00400000);
        __m256i aa = _mm256_cmpgt_epi32(am, round);
        __m256i bb = _mm256_cmpgt_epi32(bm, round);
        __m256i cc = _mm256_cmpgt_epi32(cm, round);
        __m256i dd = _mm256_cmpgt_epi32(dm, round);
#endif

        // shift right (fill zeros) so that only sign and exponent left
        // add 1 to exp accordingly (undefined behavior if float is Inf)
        const static __m256i one = _mm256_set1_epi32(1);
        __m256i aexp = _mm256_add_epi32(_mm256_srli_epi32(a, 23), _mm256_and_si256(aa, one));
        __m256i bexp = _mm256_add_epi32(_mm256_srli_epi32(b, 23), _mm256_and_si256(bb, one));
        __m256i cexp = _mm256_add_epi32(_mm256_srli_epi32(c, 23), _mm256_and_si256(cc, one));
        __m256i dexp = _mm256_add_epi32(_mm256_srli_epi32(d, 23), _mm256_and_si256(dd, one));

        // pack as 16 uint16_t
        __m256i ab = _mm256_packus_epi32(aexp, bexp);
        __m256i cd = _mm256_packus_epi32(cexp, dexp);

        // save sign and put to the MSB (8th)
        __m256i absign = _mm256_slli_epi16(_mm256_srli_epi16(ab, 8), 7);
        __m256i cdsign = _mm256_slli_epi16(_mm256_srli_epi16(cd, 8), 7);
        __m256i abcdsign = _mm256_packus_epi16(absign, cdsign);

        // mask out sign (becuase _mm256_packus_epi16 saturates signed word to unsigned byte)
        const static __m256i m2 = _mm256_set1_epi16(0x00FF);
        __m256i abmasksign = _mm256_and_si256(ab, m2);
        __m256i cdmasksign = _mm256_and_si256(cd, m2);

        // pack as 32 uint8_t
        __m256i abcd = _mm256_packus_epi16(abmasksign, cdmasksign);

        // we want to use MSB for value sign, so choose to preserve exponent -109 ~ 17
        // and treat exponent -110 (and below) as value 0, 18 (and above) as value 2^17
        //        exponent:    below    -110      -109      0      17      18    above
        // view as uint8_t:      <17      17        18    127     144     145     >145
        //  +110(saturate):     <127     127       128    237     254     255      255
        //  -127(saturate):        0       0         1    110     127     127      127
        //    decompressed:        0       0    2^-109    2^0    2^17    2^17     2^17

        // saturate to 17~145
        const static __m256i oneten = _mm256_set1_epi8(110); // 255-145
        const static __m256i onetwentyseven = _mm256_set1_epi8(127);
        __m256i abcdsaturate = _mm256_subs_epu8(_mm256_adds_epu8(abcd, oneten), onetwentyseven);
        // put value sign to MSB
        __m256i abcdsigned = _mm256_or_si256(abcdsign, abcdsaturate);

        // packing will compressed_data in different order,
        // if we do corresponding unpack when decompressing then no need to shuffle
        // const static __m256i shuffleidx = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
        // __m256i abcdordered = _mm256_permutevar8x32_epi32(abcdsigned, shuffleidx);

        _mm256_storeu_si256((__m256i *) exp, abcdsigned);
    }

    void decompress() override {
        float *result_ptr = decompressed_data.data();
        uint8_t *exp_ptr = compressed_data.data();

#pragma omp parallel for default(none) shared(result_ptr, exp_ptr)
        for (unsigned i = 0; i < num_elements / 32; ++i) {
            decompress_avx(result_ptr + i * 32, exp_ptr + i * 32);
        }

        unsigned offset = 32 * (num_elements / 32);
        result_ptr += offset;
        exp_ptr += offset;

#pragma omp parallel for default(none) shared(result_ptr, exp_ptr)
        for (unsigned i = 0; i < num_elements % 32; ++i) {
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


    static void decompress_avx(float *result, const uint8_t *exponent) {
        // load 32 uint8_t
        __m256i abcdwsign = _mm256_loadu_si256((__m256i const *) (exponent));

        // unpack to 2 * 16 uint16_t
        __m256i abwsign = _mm256_unpacklo_epi8(abcdwsign, _mm256_setzero_si256());
        __m256i cdwsign = _mm256_unpackhi_epi8(abcdwsign, _mm256_setzero_si256());

        // mask out value sign
        const static __m256i m = _mm256_set1_epi16(0x007F);
        __m256i abwosign = _mm256_and_si256(abwsign, m);
        __m256i cdwosign = _mm256_and_si256(cdwsign, m);

        // if exponent is 0, the decompressed value should be 0
        __m256i abiszero = _mm256_cmpeq_epi16(abwosign, _mm256_setzero_si256());
        __m256i cdiszero = _mm256_cmpeq_epi16(cdwosign, _mm256_setzero_si256());

        // add 17 (127-110) to get the original exponent
        const static __m256i seventeen = _mm256_set1_epi16(17);
        __m256i abexp = _mm256_add_epi16(abwosign, seventeen);
        __m256i cdexp = _mm256_add_epi16(cdwosign, seventeen);

        // move sign to 9th bit, mask out others
        __m256i absign = _mm256_slli_epi16(_mm256_srli_epi16(abwsign, 7), 8);
        __m256i cdsign = _mm256_slli_epi16(_mm256_srli_epi16(cdwsign, 7), 8);

        // put value sign back to 9th bit
        __m256i absignexp = _mm256_or_si256(abexp, absign);
        __m256i cdsignexp = _mm256_or_si256(cdexp, cdsign);

        // set zero, _mm256_andnot_si256(A, B) performs bitwise ((not A) and B)
        __m256i absignexp_with_zero = _mm256_andnot_si256(abiszero, absignexp);
        __m256i cdsignexp_with_zero = _mm256_andnot_si256(cdiszero, cdsignexp);

        // shift 23 left, cast back to float
        __m256 a = _mm256_castsi256_ps(
                _mm256_slli_epi32(_mm256_unpacklo_epi16(absignexp_with_zero, _mm256_setzero_si256()), 23));
        __m256 b = _mm256_castsi256_ps(
                _mm256_slli_epi32(_mm256_unpackhi_epi16(absignexp_with_zero, _mm256_setzero_si256()), 23));
        __m256 c = _mm256_castsi256_ps(
                _mm256_slli_epi32(_mm256_unpacklo_epi16(cdsignexp_with_zero, _mm256_setzero_si256()), 23));
        __m256 d = _mm256_castsi256_ps(
                _mm256_slli_epi32(_mm256_unpackhi_epi16(cdsignexp_with_zero, _mm256_setzero_si256()), 23));

        // store floats back
        _mm256_storeu_ps(result, a);
        _mm256_storeu_ps(result + 8, b);
        _mm256_storeu_ps(result + 16, c);
        _mm256_storeu_ps(result + 24, d);
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