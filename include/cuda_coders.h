#ifndef _CUDA_CODER_LIB_H
#define _CUDA_CODER_LIB_H

#define THREAD_PER_BLOCK 256

#ifndef PI
#define PI 3.1416
#endif

void CudaGpuInitialize();
int CudaG711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int CudaG711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int CudaG722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int CudaPcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr);
int CudaPcmToG711u(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* ulaw_data_ptr);
int CudaPcmToG722(const short int* pcm_data_ptr, short int* band_data_ptr, unsigned int no_of_data, unsigned char* g722_data_ptr);
void CudaSinusoidalSynthesis(const float* parameter_ptr, int max_sin, int number_of_packets, int samples_per_packet, short int* output_ptr);
bool DetectCudaDevice();
#endif
