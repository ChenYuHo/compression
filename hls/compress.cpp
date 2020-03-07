#include "compress.hpp"

//void compress_one(float *value, ap_int<8> *compressed) {
void compress_one(hls::stream<apuint32_t> &in_stream, hls::stream<apuint8_t> &out_stream) {
#pragma HLS interface axis port=in_stream
#pragma HLS interface axis port=out_stream


#ifndef __SYNTHESIS__
//	printf("%f  ",*value);
#endif
//	float inval = *value;
//	ap_uint<32> intVal = *(ap_uint<32>*) &inval;
	apuint32_t intVal;
	apuint8_t compressed;

	in_stream.read(intVal);
	ap_uint<23> mant = intVal(22,0);
	ap_int<8> exponent = intVal(30,23)-127;
	bool signbit = intVal[31];
#ifndef __SYNTHESIS__
//	printf("%d %d %d\n",signbit, (int8_t) exponent, (uint32_t) mant);
#endif

	if (intVal[22]==1) exponent+=1; // add 1

	if (exponent >= 17) compressed = 127;
	else if (exponent <= -110 || ((mant == 0) && (exponent == 0))) compressed = 0;
	else compressed = (apuint8_t)(exponent + 110);
	if (signbit) compressed |= (apuint8_t) 0x80u;

	out_stream.write(compressed);

}
