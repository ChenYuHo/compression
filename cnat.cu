#include <iostream>
#include <cuda_runtime.h>
#include <stdint.h>
#include <math.h>
#include <curand.h>

using namespace std;
#define CUDA_CALL(x) do { \
      cudaError_t _m_cudaStat = x; \
      if((_m_cudaStat) != cudaSuccess) { \
      fprintf(stderr, "Error %s at line %d in file %s", \
      cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); \
      exit(1);}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
      printf("CURAND rrror at %s:%d\n",__FILE__,__LINE__);            \
      exit(1);}} while(0)

__global__ void cnat_compress_cuda_kernel(
    float* __restrict__ input,
    uint8_t* __restrict__ output,
    int len) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if(index < len) {
    if (input[index] == 0)
      output[index] = 0;
    else {
      int rand_bits = output[index];
      rand_bits <<= 15;
      int exp;
      int prob = reinterpret_cast<int &>(input[index]) & 0b00000000011111111000000000000000;
      frexpf(input[index], &exp);
      if (rand_bits >= prob) exp -= 1;
      exp += 127;
      uint8_t encode;
      if (exp<=17) encode = 0;
      else if (exp<=143) encode = uint8_t(exp-17);
      else encode = 127;
      if (input[index] < 0) encode += 128;
      output[index] = encode;
    }
  }
}

/*
  generate random bits using the output buffer (8 bits for each element)
  and compare with the 8 most significant bits of mantissa
*/
void cnat_compress(float* input, uint8_t* output, int count, cudaStream_t stream, curandGenerator_t* gen) {
    const int threads = 1024;
    auto blocks = count/threads;
    if (count%threads || !blocks) blocks++;
    CURAND_CALL(curandSetStream(*gen, stream));
    CURAND_CALL(curandGenerate(*gen, (uint32_t*)output, count/4)); // (0, 1]
    cnat_compress_cuda_kernel<<<blocks, threads, 0, stream>>>(
            input,
            output,
            count);
}

__global__ void cnat_decompress_cuda_kernel(
    uint8_t* __restrict__ input,
    float* __restrict__ output,
    int len) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if(index < len) {
      uint32_t decode;
      if (input[index]<=127) decode = input[index]+17;
      else decode = input[index]+145;
      if (!input[index] % 128) decode -= 17;
      uint32_t sign_and_exp = decode << 23;
      output[index] = reinterpret_cast<float &>(sign_and_exp);
  }
}

void cnat_decompress(uint8_t* input, float* output, int count, cudaStream_t stream) {
  const int threads = 1024;
  auto blocks = count/threads;
  if (count%threads || !blocks) blocks++;
  cnat_decompress_cuda_kernel<<<blocks, threads, 0, stream>>>(
    input,
    output,
    count);
}

int main(){
    curandGenerator_t gen;
    cudaStream_t stream;
    CUDA_CALL(cudaStreamCreate(&stream));
    CURAND_CALL(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, time(NULL)));
    float host[] = {3., 3., 3., 3., 3., 3., 3., 3.};
    float* ptr;
    uint8_t* output;
    uint8_t host_output[8];
    for(int i=0; i<8; ++i) cout<<host[i]<<", ";
    cout<<endl;
    CUDA_CALL(cudaMalloc((void **)&ptr, 8*sizeof(float)));
    CUDA_CALL(cudaMalloc((void **)&output, 8));
    cudaMemcpy(ptr, host, 8*sizeof(float), cudaMemcpyHostToDevice);
    cnat_compress(ptr, output, 8, stream, &gen);
    cnat_decompress(output, ptr, 8, stream);
    cudaStreamSynchronize(stream);
    cudaMemcpy(host_output, output, 8, cudaMemcpyDeviceToHost);
    cudaMemcpy(host, ptr, 32, cudaMemcpyDeviceToHost);
    cudaStreamSynchronize(stream);
    for(int i=0; i<8; ++i) cout<<int(host_output[i])<<", ";
    cout<<endl;
    for(int i=0; i<8; ++i) cout<<(host[i])<<", ";
    cout<<endl;
    CURAND_CALL(curandDestroyGenerator(gen));
    return 0;
}
