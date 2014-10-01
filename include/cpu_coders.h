#ifndef _CPU_CODER_LIB_H
#define _CPU_CODER_LIB_H

#ifndef PI
#define PI 3.1416
#endif

int G711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int G711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int G722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);

int PcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr);
int PcmToG711u(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* ulaw_data_ptr);
int PcmToG722(const short int* pcm_data_ptr, short int* band_data_ptr, unsigned int no_of_data, unsigned char* g722_data_ptr);

void SinusoidalSynthesis(const float* parameter_ptr, int max_sin, int number_of_packets, int samples_per_packet, short int* output_ptr);

#endif
