#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

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
        original_data.assign(num_elements, 0);
        init_compressed_data();
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
        int8_t exp = ilogbf(*value);
        if (1.5 < get_prob(*value)) exp+=1; // add 1
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

void do_work(Compressor<float> &compressor) {
    compressor.generate_data();
    uint32_t count = 10;
    uint32_t num_elements = 26214400;
        for (unsigned i = 1; i <= count; ++i) {
            compressor.init_compressed_data();
            high_resolution_clock::time_point begin = high_resolution_clock::now();
            compressor.compress();
            high_resolution_clock::time_point end = high_resolution_clock::now();
            auto ns_taken = duration_cast<nanoseconds>(end - begin).count();
            cout << "compress round " << i << ": " << ns_taken << " nanoseconds elapsed, "
                 << double(num_elements) / ns_taken * 1e9 << " element/s\n";
        }
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

int main(int argc, char *argv[]) {
    IntMLCompressor compressor(26214400);
    do_work(compressor);
    return 0;
}
