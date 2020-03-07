#include "ap_int.h"
#include "hls_stream.h"

typedef ap_uint<32> apuint32_t;
//typedef ap_uint<8> uint8_t;

//void compress_one(float *value, ap_int<8> *compressed);
void compress_one(hls::stream<apuint32_t> &in_stream, hls::stream<uint8_t> &out_stream);
