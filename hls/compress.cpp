#include "compress.hpp"

void compress_one(float *value, ap_int<8> *compressed) {
#ifndef __SYNTHESIS__
//	printf("%f  ",*value);
#endif
	float inval = *value;
	ap_uint<32> intVal = *(ap_uint<32>*) &inval;
	ap_uint<23> mant = intVal(22,0);
	ap_int<8> exponent = intVal(30,23)-127;
	bool signbit = intVal[31];
#ifndef __SYNTHESIS__
//	printf("%d %d %d\n",signbit, (int8_t) exponent, (int32_t) mant);
#endif

	if (intVal[22]==1) exponent+=1; // add 1

	if (exponent >= 17) *compressed = 127;
	else if (exponent <= -110 || *value == 0) *compressed = 0;
	else *compressed = (ap_int<8>)(exponent + 110);
	if (signbit) *compressed |= (ap_int<8>) 0x80u;

}
