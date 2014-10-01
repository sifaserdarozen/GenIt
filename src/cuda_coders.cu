#include <iostream>
#include <vector>
#include <cuda.h>
#include "cuda_coders.h"

#include "g722coder.h"

#define checkCudaErrors(val) CheckErrors( (val), __FILE__, __LINE__)

__device__ short int* q6_ptr;
__device__ short int * misil_ptr;
__device__ short int * coef_qmf_ptr;
__device__ short int * ril4_ptr;
__device__ short int * risil_ptr;
__device__ short int * oq4_ptr;
__device__ short int * wl_ptr;
__device__ short int * ila_ptr;
__device__ short int * misih_ptr;
__device__ short int * ih2_ptr;
__device__ short int * sih_ptr;
__device__ short int * oq2_ptr;
__device__ short int * wh_ptr;


void CheckErrors(cudaError_t cuda_error, const char* const file, const int line)
{
	if (cuda_error != cudaSuccess)
	{
		std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
		std::cerr << cudaGetErrorString(cuda_error) << std::endl;
		exit(1);
	}		
}

__global__ void CudaKernelG711aToPcm(unsigned char* d_alaw_data_ptr, short int* d_pcm_data_ptr)
{
	unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);
	
	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;

	for (int k=0; k<160; k++)
	{
		alaw_data = d_alaw_data_ptr[idx+k];	
		alaw_data^=0x55;

		quantization_value= (alaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)alaw_data & (0x70)) >> (4);
		switch (quantization_segment)
		{
		case 0: 
			quantization_value+=(0x0008);
			break;
		case 1:
			quantization_value+=(0x0108);
			break;
		default:
			quantization_value+=(0x0108);
			quantization_value <<= (quantization_segment-1);
		};

		d_pcm_data_ptr[idx+k]=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

// memory coalesced version of alaw to pcm conversion
__global__ void CudaKernelG711aToPcmCM(unsigned char* d_alaw_data_ptr, short int* d_pcm_data_ptr)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;

	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
		alaw_data = d_alaw_data_ptr[idx];	
		alaw_data^=0x55;

		quantization_value= (alaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)alaw_data & (0x70)) >> (4);
		switch (quantization_segment)
		{
		case 0: 
			quantization_value+=(0x0008);
			break;
		case 1:
			quantization_value+=(0x0108);
			break;
		default:
			quantization_value+=(0x0108);
			quantization_value <<= (quantization_segment-1);
		};

		d_pcm_data_ptr[idx]=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

// memory coalesced version of pcm to alaw conversion
__global__ void CudaKernelPcmToG711aCM(short int* d_pcm_data_ptr, unsigned char* d_alaw_data_ptr)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;

	short int quantization_value;
	short int quantization_segment;
	short int pcm_data;
	unsigned char alaw_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
		alaw_data = 0;
		pcm_data = d_pcm_data_ptr[idx];
		quantization_value=(pcm_data<0) ? ((~pcm_data)>>4) : (pcm_data>>4);
		
		if(quantization_value>15)
		{
			quantization_segment=1;
			while(quantization_value>(16+15))
			{
				quantization_value>>=1;
				quantization_segment++;
			}
			quantization_value-=16;

			alaw_data=quantization_value + (quantization_segment << 4);
		}

		if(pcm_data>=0)
			alaw_data |= 0x80;

		alaw_data^=0x55;
		d_alaw_data_ptr[idx] = alaw_data;
	}
}

__global__ void CudaKernelG711uToPcm(unsigned char* d_ulaw_data_ptr, short int* d_pcm_data_ptr)
{
	unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);

	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;

	for (int k=0; k<160; k++)
	{
		ulaw_data=~(d_ulaw_data_ptr[idx+k]);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);

		d_pcm_data_ptr[idx+k]=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

__global__ void CudaKernelG711uToPcmCM(unsigned char* d_ulaw_data_ptr, short int* d_pcm_data_ptr)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;

	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
		ulaw_data=~(d_ulaw_data_ptr[idx]);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);

		d_pcm_data_ptr[idx]=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

__global__ void CudaKernelPcmToG711uCM(short int* d_pcm_data_ptr, unsigned char* d_ulaw_data_ptr)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;
		
	short int quantization_value;
	short int quantization_segment = 1;
	unsigned char ulaw_data;
	short int pcm_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
		pcm_data = d_pcm_data_ptr[idx];

		quantization_value=(pcm_data<0) ? (((~pcm_data)>>2)+33) : ((pcm_data>>2)+33);

		if (quantization_value > (0x1FFF))	// clip to 8192
			quantization_value = (0x1FFF);

		quantization_segment = 1;
		// Determination of quantization segment
		for (short int i = (quantization_value >> 6); i; i>>= 1)
			quantization_segment++;

		ulaw_data =  (((0x08 - quantization_segment) << 4) | (0x000F - ((quantization_value  > quantization_segment) & 0x000F)));

		if (pcm_data >= 0)
			ulaw_data |= 0x80;

		d_ulaw_data_ptr[idx] = ulaw_data;
	}
}



/*
__global__ void CudaKernelG722ToPcm(unsigned char* d_g722_data_ptr, short int* d_pcm_data_ptr)
{
	unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);
	
	unsigned char g722_data;

	for (int k=0; k<160; k++)
	{
		g722_data=d_g722_data_ptr[idx+k];

		d_pcm_data_ptr[idx+k]=0;
	}
}
*/

__device__ short int CudaConvertLongToShort(int in_value)
{
	if (in_value > 32767)
		return 32767;
	else if (in_value < -32768)
		return -32768;
	else
		return (short)in_value;
}

__global__ void CudaKernelG722ToPcmCM(unsigned char* d_g722_data_ptr, short int* d_pcm_data_ptr, int* d_band_data_ptr, int* d_g722_consts_ptr, unsigned int no_of_data)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
//	unsigned int total_threads = blockDim.x * gridDim.x;

	unsigned char g722_data;
	
	int number_of_chunks = no_of_data/160;
	if (idx >= number_of_chunks)
		return;

	// pointers for constants, maybe copy these to shared mem
	int* wl_dec = d_g722_consts_ptr;
	int* rl42 = d_g722_consts_ptr + 8;
	int* ilb = rl42 + 16;
	int* qm4 = ilb + 32;
	int* qm6 = qm4 + 16;

	// copy band data to local variables
	int band_s = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_sp = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_sz = d_band_data_ptr[idx];
	idx += number_of_chunks;

	int band_r[3], band_a[3], band_ap[3], band_p[3];
	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
		band_r[k] = d_band_data_ptr[idx];
		band_a[k] = d_band_data_ptr[idx + 3*number_of_chunks];
		band_ap[k] = d_band_data_ptr[idx + 6*number_of_chunks];
		band_p[k] = d_band_data_ptr[idx + 9*number_of_chunks];
	}
	idx += 9*number_of_chunks;

	int band_d[7], band_b[7], band_bp[7], band_sg[7];
	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
		band_d[k] = d_band_data_ptr[idx];
		band_b[k] = d_band_data_ptr[idx + 7*number_of_chunks];
		band_bp[k] = d_band_data_ptr[idx + 14*number_of_chunks];
		band_sg[k] = d_band_data_ptr[idx + 21*number_of_chunks];
	}
	idx += 21*number_of_chunks;

	int band_nb = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_det = d_band_data_ptr[idx];
	//band_det=32;

	int dlowt;
	int rlow;
	int wd1;
	int wd2;
	int wd3;

	idx = threadIdx.x + blockDim.x * blockIdx.x;
	for (int k=0; k<160; k++, idx+=number_of_chunks)
	{
		g722_data=d_g722_data_ptr[idx];

		wd1 = g722_data & 0x3F;
		wd2 = qm6[wd1];
		wd1 >>= 2;

		/********************** Block 5 *******************/
		 // INVQBL (ITU page 43), compute quantized difference signal for the decoder output in the lower sub-band
		 wd2 = (band_det * wd2) >> 15;
		 // RECONS ( ITU page 41), compute reconstructed signal for the adaptive predictor
		 rlow = band_s + wd2;

		/********************** Block 6 ********************/
		// LIMIT (ITU page 44), limit the output reconstructed signal
		if (rlow > 16383)
			rlow = 16383;
		else if (rlow < -16384)
			 rlow = -16384;

		/********************** Block 2 ***********************/	
		// INVQAL (ITU page 37), compute the quantized differences signal for the adaptive predictor in the lower sub-band
		wd2 = qm4[wd1];
		dlowt = (band_det * wd2) >> 15;

		/********************** Block 3 ************************/
		// LOGSCL (ITU page 38), update the logarithmic quantizer scale factor in the lower sub-band
		wd2 = rl42[wd1];
		 wd1 = (band_nb * 127) >> 7;
		wd1 += wl_dec[wd2];
		if (wd1 < 0)
			wd1 = 0;
		else if (wd1 > 18432)
			wd1 = 18432;
		band_nb = wd1;

		// SCALEL (ITU page 38), compute the quantizer scale factor in the lower sub-band 
		wd1 = (band_nb >> 6) & 31;
		wd2 = 8 - (band_nb >> 11);
		wd3 = (wd2 < 0)	 ?  (ilb[wd1] << -wd2)	:  (ilb[wd1] >> wd2);
		band_det = wd3 << 2;

		/********************** Block 4 **************************/

		// RECONS (ITU page 41), compute reconstructed signal for the adaptive predictor
		band_d[0] = dlowt;
		 band_r[0] = CudaConvertLongToShort(band_s + dlowt);

		 // PARREC (ITU page 40), compute partially reconstructed signal
		band_p[0] = CudaConvertLongToShort(band_sz + dlowt);

		// UPPOL2 (ITU page 41), update second predictor coefficient
		int i;  // loop variable
		for (i = 0;	 i < 3;	 i++)
			band_sg[i] = band_p[i] >> 15;
		wd1 = CudaConvertLongToShort(band_a[1] << 2);

		wd2 = (band_sg[0] == band_sg[1])	?  -wd1	 :  wd1;
		if (wd2 > 32767)
			wd2 = 32767;
		wd3 = (band_sg[0] == band_sg[2])	?  128	:  -128;
		wd3 += (wd2 >> 7);
		wd3 += (band_a[2]*32512) >> 15;
		if (wd3 > 12288)
			wd3 = 12288;
		else if (wd3 < -12288)
			wd3 = -12288;
		band_ap[2] = wd3;

		// UPPOL1 (ITU page 42), update first predictor coefficient
		band_sg[0] = band_p[0] >> 15;
		band_sg[1] = band_p[1] >> 15;
		wd1 = (band_sg[0] == band_sg[1])	?  192	:  -192;
		wd2 = (band_a[1]*32640) >> 15;

		band_ap[1] = CudaConvertLongToShort(wd1 + wd2);
		 wd3 = CudaConvertLongToShort(15360 - band_ap[2]);
		if (band_ap[1] > wd3)
			band_ap[1] = wd3;
		else if (band_ap[1] < -wd3)
			 band_ap[1] = -wd3;

		// UPZERO (ITU page 41), update sixth order predictor coefficients
		wd1 = (dlowt == 0)  ?  0  :  128;
		band_sg[0] = dlowt >> 15;
		for (i = 1;	 i < 7;	 i++)
		{
			band_sg[i] = band_d[i] >> 15;
			wd2 = (band_sg[i] == band_sg[0])  ?  wd1  :  -wd1;
			wd3 = (band_b[i]*32640) >> 15;
			band_bp[i] = CudaConvertLongToShort(wd2 + wd3);
		}

		// DELAYA (ITU page 38), memory block delay 
		for (i = 6;	 i > 0;	 i--)
		{
			 band_d[i] = band_d[i - 1];
			band_b[i] = band_bp[i];
		}

		for (i = 2;	 i > 0;	 i--)
		{
			band_r[i] = band_r[i - 1];
			band_p[i] = band_p[i - 1];
			band_a[i] = band_ap[i];
		}

		// FILTEP (ITU page 43), compute predictor output signal, poles
		wd1 = CudaConvertLongToShort(band_r[1] + band_r[1]);
		wd1 = (band_a[1]*wd1) >> 15;
		wd2 = CudaConvertLongToShort(band_r[2] + band_r[2]);
		wd2 = (band_a[2]*wd2) >> 15;
		band_sp = CudaConvertLongToShort(wd1 + wd2);

		// FILTEZ (ITU page 42), compute predictor output signal, zeros
		band_sz = 0;
		for (i = 6;	 i > 0;	 i--)
		{
			wd1 = CudaConvertLongToShort(band_d[i] + band_d[i]);
			band_sz += (band_b[i]*wd1) >> 15;
		}
		band_sz = CudaConvertLongToShort(band_sz);

		// PREDIC (ITU page 43), compute predictor output value
		band_s = CudaConvertLongToShort(band_sp + band_sz);

		d_pcm_data_ptr[idx]=(short int)rlow;
	}

	// copy local variables back to band data
	idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	d_band_data_ptr[idx] = band_s;
	idx += number_of_chunks;
	band_sp = d_band_data_ptr[idx]= band_sp;
	idx += number_of_chunks;
	d_band_data_ptr[idx] = band_sz;
	idx += number_of_chunks;
	
	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
		d_band_data_ptr[idx] = band_r[k];
		d_band_data_ptr[idx + 3*number_of_chunks] = band_a[k];
		d_band_data_ptr[idx + 6*number_of_chunks] = band_ap[k];
		d_band_data_ptr[idx + 9*number_of_chunks] = band_p[k];
	}
	idx += 9*number_of_chunks;

	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
		d_band_data_ptr[idx] = band_d[k];
		d_band_data_ptr[idx + 7*number_of_chunks] = band_b[k];
		d_band_data_ptr[idx + 14*number_of_chunks] = band_bp[k];
		d_band_data_ptr[idx + 21*number_of_chunks] = band_sg[k];
	}
	idx += 21*number_of_chunks;

	d_band_data_ptr[idx] = band_nb;
	idx += number_of_chunks;
	d_band_data_ptr[idx] = band_det;
	idx += number_of_chunks;
}

// *******************************************************************************************************



__device__ int CudaSaturateAdd(int op1, int op2)
{
	int out = op1 + op2;
	if ((((op1 ^ op2) & MIN_32) == 0) && ((out ^ op1) & MIN_32))
		out = (op1 < 0) ? MIN_32 : MAX_32;
	return out;
}

__device__ int CudaSaturateSubtract(int op1, int op2)
{
	int out = op1 - op2;
	if ((((op1 ^ op2) & MIN_32) != 0) && ((out ^ op1) & MIN_32))
		out = (op1 < 0L) ? MIN_32 : MAX_32;
	return out;
}

__device__ int CudaShiftRight(int op1, short int op2)
{
	if (op2 > 0)
	{
		if (op2 >= 31)
			return (op1 < 0) ? -1 : 0;
		else
			return op1 >> op2;
	}
	return op1;
}

__device__ int CudaShiftLeft(int op1, short int op2)
{
	if (op2 > 0)
		for (; op2 > 0; op2--)
		{
			if (op1 > 0X3fffffff)
				return MAX_32;
			else if (op1 < (int)0xc0000000)
				return MIN_32;
			op1 *= 2;
		}
	return op1;
}


__device__ short int CudaShiftLeftShort(short int op1, short int op2)
{
	if (op1 > 0)
	{
		int result = ((int)op1) * ((int) 1 << op2);

		if ((op2 > 15 && op1 != 0) || (result != (int) ((short int) result)))
			return (op1 > 0) ? MAX_16 : MIN_16;
		else
			return (short int)result;
	}
	return op1;
}

__device__ short int CudaShiftRightShort(short int op1, short int op2)
{
	if (op2 > 0)
	{
		if (op2 >= 15)
			return (op1 < 0) ? -1 : 0;
		else
			if (op1 < 0)
				return ~((~op1) >> op2);
			else
				return op1 >> op2;
	}
	return op1;
}

__device__ int CudaClamp15ToBits(int op)
{
	if (op > 16383)
		return 16383;
	else if (op < -16384)
		return -16384;
	return op;
}

__device__ int CudaMultiplyAdd(int add_op, short int mul_op1, short int mul_op2)
{
	return CudaSaturateAdd(add_op, ((int)mul_op1 * (int)mul_op2));
}

__device__ short int CudaSaturate(int op)
{
	if (op > MAX_16)
		return MAX_16;
	else if (op < MIN_16)
		return MIN_16;
	return op;
}

__device__ short int CudaSaturateSubtractShort(short int op1, short int op2)
{
	return CudaSaturate (((int)op1 - op2));
}

__device__ short int CudaSaturateAddShort(short int op1, short int op2)
{
	return CudaSaturate (((int)op1 + op2));
}

__device__ short int CudaScaledMult(short int op1, short int op2)
{
	int product = (((int)op1 * (int)op2) & (int)(0xffff8000)) >> 15;

	if (product & (int)0x00010000)
		product |= (int)0xffff0000;

	return (short int)product;
}

__device__ short int CudaQuantl(short int el, short int detl)
{
	short int sil = CudaShiftRightShort(el, 15);
	short int wd = CudaSaturateSubtractShort(MAX_16,(el & MAX_16));
	short int mil = 0;

	if (sil == 0)
		wd = el;

	short int val = CudaScaledMult(CudaShiftLeftShort(q6_ptr[mil], 3), detl);
	while (CudaSaturateSubtractShort(val,wd) <= 0)
	{
		if (CudaSaturateSubtractShort(mil, 30) == 0)
			break;
		else
		{
			mil = CudaSaturateAddShort(mil, 1);
			val = CudaScaledMult(CudaShiftLeftShort(q6_ptr[mil], 3), detl);
		}
	}

	sil = CudaSaturateAddShort(sil, 1);

	return misil_ptr[sil*32 + mil];
}

__device__ short int CudaQuanth(short int eh, short int deth)
{
	short int sih = CudaShiftRightShort(eh, 15);
	short int wd = CudaSaturateSubtractShort(MAX_16, (eh & MAX_16));

	if (sih == 0)
		wd = eh;

	short int mih = 1;

	if (CudaSaturateSubtractShort(wd, CudaScaledMult(CudaShiftLeftShort(564, 3), deth)) >= 0)
		mih = 2;

	sih = CudaSaturateAddShort(sih, 1);

	return misih_ptr[sih*3 + mih];
}

__device__ short int CudaInvqal(short int il, short int detl)
{
	short int ril = CudaShiftRightShort(il, 2);
	short int wd1 = CudaShiftLeftShort(oq4_ptr[ril4_ptr[ril]], 3);
	short int wd2 = -wd1;

	if (risil_ptr[ril] == 0)
		wd2 = wd1;

	return CudaScaledMult(detl, wd2);
}

__device__ short int CudaInvqah(short int ih, short int deth)
{
	short int wd1 = CudaShiftLeftShort(oq2_ptr[ih2_ptr[ih]], 3);
	short int wd2 = -wd1;

	if (sih_ptr[ih] == 0)
		wd2 = wd1;

	return CudaScaledMult(wd2, deth);
}

__device__ short int CudaLogscl(short int il, short int nbl)
{
	short int ril = CudaShiftRightShort(il, 2);
	short int wd = CudaScaledMult(nbl, 32512);
	short int il4 = ril4_ptr[ril];
	short int nbpl = CudaSaturateAddShort (wd, wl_ptr[il4]);

	if (nbpl < 0)
		nbpl = 0;

	if (CudaSaturateSubtractShort(nbpl, 18432) > 0)
		nbpl = 18432;

	return nbpl;
}

__device__ short int CudaLogsch(short int ih, short int nbh)
{
	short int wd = CudaScaledMult(nbh, 32512);
	short int nbph = CudaSaturateAddShort(wd, wh_ptr[ih2_ptr[ih]]);

	if(nbph < 0)
		nbph = 0;

	if(CudaSaturateSubtractShort(nbph, 22528) > 0)
		nbph = 22528;

	return nbph;
}

__device__ short int CudaScalel(short int nbpl)
{
	short int wd1 = CudaShiftRightShort(nbpl, 6) & 511;
	short int wd2 = CudaSaturateAddShort(wd1, 64);
	return (CudaShiftLeftShort(CudaSaturateAddShort(ila_ptr[wd2], 1), 2));
}

__device__ short int CudaScaleh(short int nbph)
{
	short int wd = CudaShiftRightShort(nbph, 6) & 511;
	return CudaShiftLeftShort(CudaSaturateAddShort(ila_ptr[wd], 1), 2);
}


__device__ void CudaUpzero(short int* dlt_ptr, short int* bl_ptr)
{
	short int wd1 = 128;

	if (dlt_ptr[0] == 0)
		 wd1 = 0;

	short int sg0 = CudaShiftRightShort(dlt_ptr[0], 15);

	for (short int i = 6; i > 0; i--)
	{
		short int wd2 = CudaSaturateSubtractShort (0, wd1);
		if(sg0 == CudaShiftRightShort(dlt_ptr[i], 15))
			wd2 = CudaSaturateAddShort (0, wd1);

		bl_ptr[i] = CudaSaturateAddShort(wd2, CudaScaledMult(bl_ptr[i], 32640));
		dlt_ptr[i] = dlt_ptr[i - 1];
	}
}

__device__ void CudaUppol1(short int* al_ptr, short int* plt_ptr)
{
	short int sg0 = CudaShiftRightShort(plt_ptr[0], 15);
	short int sg1 = CudaShiftRightShort(plt_ptr[1], 15);
	short int wd1 = -192;

	if (CudaSaturateSubtractShort(sg0, sg1) == 0)
		wd1 = 192;

	short int wd2 = CudaScaledMult (al_ptr[1], 32640);
	short int apl1 = CudaSaturateAddShort(wd1, wd2);
	short int wd3 = CudaSaturateSubtractShort(15360, al_ptr[2]);

	if (CudaSaturateSubtractShort(apl1, wd3) > 0)
		apl1 = wd3;
	else if (CudaSaturateAddShort(apl1, wd3) < 0)
		apl1 = -wd3;

	/* Shift of the plt signals */
	plt_ptr[2] = plt_ptr[1];
	plt_ptr[1] = plt_ptr[0];
	al_ptr[1] = apl1;
}

__device__ void CudaUppol2(short int* al_ptr, short int* plt_ptr)
{
	short int sg0 = CudaShiftRightShort(plt_ptr[0], 15);
	short int sg1 = CudaShiftRightShort(plt_ptr[1], 15);
	short int sg2 = CudaShiftRightShort(plt_ptr[2], 15);
	short int wd1 = CudaShiftLeftShort(al_ptr[1], 2);
	short int wd2 = CudaSaturateAddShort(0, wd1);

	if (CudaSaturateSubtractShort(sg0, sg1) == 0)
		wd2 = CudaSaturateSubtractShort(0, wd1);

	wd2 = CudaShiftRightShort(wd2, 7);
	short int wd3 = -128;

	if (CudaSaturateSubtractShort(sg0, sg2) == 0)
		wd3 = 128;

	short int wd4 = CudaSaturateAddShort (wd2, wd3);
	short int wd5 = CudaScaledMult(al_ptr[2], 32512);
	short int apl2 = CudaSaturateAddShort(wd4, wd5);

	if (CudaSaturateSubtractShort(apl2, 12288) > 0)
		apl2 = 12288;

	if (CudaSaturateSubtractShort(apl2, -12288) < 0)
		apl2 = -12288;

	al_ptr[2] = apl2;
}

__device__ short int CudaFiltez(short int* dlt_ptr, short int* bl_ptr)
{
	short int szl = 0;

	for (short int i = 6; i > 0; i--)
	{
		short int wd = CudaSaturateAddShort(dlt_ptr[i], dlt_ptr[i]);
		wd = CudaScaledMult(wd, bl_ptr[i]);
		szl = CudaSaturateAddShort(szl, wd);
	}
	return szl;
}

__device__ short int CudaFiltep(short int* rlt_ptr, short int* al_ptr)
{
	// shift of rlt
	rlt_ptr[2] = rlt_ptr[1];		
	rlt_ptr[1] = rlt_ptr[0];		

	short int wd1 = CudaSaturateAddShort(rlt_ptr[1], rlt_ptr[1]);
	wd1 = CudaScaledMult(al_ptr[1], wd1);
	short int wd2 = CudaSaturateAddShort(rlt_ptr[2], rlt_ptr[2]);
	wd2 = CudaScaledMult(al_ptr[2], wd2);
	return CudaSaturateAddShort(wd1, wd2);
}

__device__ void CudaQmfTx(short int xin0, short int xin1, short int& xl, short int& xh, short int* band_qmf_tx_delayx)
{
	int accuma;
	int accumb;
	int comp_low;
	int comp_high;

	const short int* pcoef = coef_qmf_ptr;
	short int* pdelayx = band_qmf_tx_delayx;

	/* Saving past samples in delay line */
	band_qmf_tx_delayx[1] = xin1;
	band_qmf_tx_delayx[0] = xin0;

	accuma = (int)*pcoef++, (int)*pdelayx++;
	accumb = (int)*pcoef++, (int)*pdelayx++;

	for(short int i = 1; i < 12; i++)
	{
		accuma = CudaMultiplyAdd(accuma, *pcoef++, *pdelayx++);
		accumb = CudaMultiplyAdd(accumb, *pcoef++, *pdelayx++);
	}

	/* Descaling and shift of the delay line */
	for (short int i = 0; i < 22; i++)
		band_qmf_tx_delayx[23 - i] = band_qmf_tx_delayx[21 - i];

	comp_low = CudaSaturateAdd (accuma, accumb);
	comp_low = CudaSaturateAdd (comp_low, comp_low);
	comp_high = CudaSaturateSubtract (accuma, accumb);
	comp_high = CudaSaturateAdd (comp_high, comp_high);
	xl = CudaClamp15ToBits (CudaShiftRight(comp_low, 16));
	xh = CudaClamp15ToBits (CudaShiftRight(comp_high, 16));
}

__device__ short int CudaLsbCod(short int xl, short int* band_dlt, short int* band_plt, short int* band_rlt, 
		short int& band_sl, short int& band_detl, short int& band_nbl, short int& band_szl, short int* band_bl, 
		short int* band_al, short int& band_spl)
{
	short int il = CudaQuantl (CudaSaturateSubtractShort (xl, band_sl), band_detl);
	band_dlt[0] = CudaInvqal (il, band_detl);
	short int nbpl = CudaLogscl (il, band_nbl);
	band_nbl = nbpl;
	band_detl = CudaScalel (nbpl);
	band_plt[0] = CudaSaturateAddShort (band_dlt[0], band_szl);   /* parrec */
	band_rlt[0] = CudaSaturateAddShort (band_sl, band_dlt[0]);    /* recons */
	CudaUpzero (band_dlt, band_bl);
	CudaUppol2 (band_al, band_plt);
	CudaUppol1 (band_al, band_plt);
	band_szl = CudaFiltez(band_dlt, band_bl);
	band_spl = CudaFiltep(band_rlt, band_al);
	band_sl = CudaSaturateAddShort (band_spl, band_szl);          /* predic */

	/* Return encoded sample */
	return il;
}

__device__ short int CudaHsbCod(short int xh, short int* band_dh, short int* band_ph, short int* band_rh,
			short int& band_sh, short int& band_deth, short int& band_nbh, short int& band_szh,
			short int* band_bh, short int* band_ah, short int& band_sph)
{
	short int ih = CudaQuanth (CudaSaturateSubtractShort(xh, band_sh), band_deth);
	band_dh[0] = CudaInvqah (ih, band_deth);
	short int nbph = CudaLogsch (ih, band_nbh);
	band_nbh = nbph;
	band_deth = CudaScaleh (nbph);
	band_ph[0] = CudaSaturateAddShort(band_dh[0], band_szh);   /* parrec */
	band_rh[0] = CudaSaturateAddShort(band_sh, band_dh[0]);    /* recons */
	CudaUpzero (band_dh, band_bh);
	CudaUppol2 (band_ah, band_ph);
	CudaUppol1 (band_ah, band_ph);
	band_szh = CudaFiltez (band_dh, band_bh);
	band_sph = CudaFiltep (band_rh, band_ah);
	band_sh = CudaSaturateAddShort(band_sph, band_szh);        /* predic */

	return ih;
}


// ******************************************************************************************************

__global__ void CudaKernelPcmToG722CM(short int* d_pcm_data_ptr, unsigned char* d_g722_data_ptr, short int* d_band_data_ptr, 
short int* d_g722_consts_ptr, unsigned int no_of_data)
{
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
//	unsigned int total_threads = blockDim.x * gridDim.x;

	int number_of_chunks = no_of_data/320;
	if (idx >= number_of_chunks)
		return;

	// pointers for constants, maybe copy these to shared mem
	coef_qmf_ptr = d_g722_consts_ptr;
	misil_ptr = d_g722_consts_ptr + 24;
	q6_ptr = d_g722_consts_ptr + 88;
	ril4_ptr = d_g722_consts_ptr + 119;
	risil_ptr = d_g722_consts_ptr + 135;
	oq4_ptr = d_g722_consts_ptr + 151;
	wl_ptr = d_g722_consts_ptr + 159;
	ila_ptr = d_g722_consts_ptr + 167;
	misih_ptr = d_g722_consts_ptr + 520;
	ih2_ptr = d_g722_consts_ptr + 526;
	sih_ptr = d_g722_consts_ptr + 530;
	oq2_ptr = d_g722_consts_ptr + 534;
	wh_ptr = d_g722_consts_ptr + 537;

	// copy band data to local variables
	// detl
	short int band_detl = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// deth
	short int band_deth = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// nbl
	short int band_nbl = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// sl
	short int band_sl = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// spl
	short int band_spl = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// szl
	short int band_szl = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// nbh
	short int band_nbh = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// sh
	short int band_sh = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// sph
	short int band_sph = d_band_data_ptr[idx];
	idx += number_of_chunks;
	// szh
	short int band_szh = d_band_data_ptr[idx];
	idx += number_of_chunks;


	short int band_al[3], band_plt[3], band_rlt[3], band_ah[3], band_ph[3], band_rh[3];
	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
		band_al[k] = d_band_data_ptr[idx];
		band_plt[k] = d_band_data_ptr[idx + 3*number_of_chunks];
		band_rlt[k] = d_band_data_ptr[idx + 6*number_of_chunks];
		band_ah[k] = d_band_data_ptr[idx + 9*number_of_chunks];
		band_ph[k] = d_band_data_ptr[idx + 12*number_of_chunks];
		band_rh[k] = d_band_data_ptr[idx + 15*number_of_chunks];
	}
	idx += 15*number_of_chunks;

	short int band_bl[7], band_dlt[7], band_bh[7], band_dh[7];
	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
		band_bl[k] = d_band_data_ptr[idx];
		band_dlt[k] = d_band_data_ptr[idx + 7*number_of_chunks];
		band_bh[k] = d_band_data_ptr[idx + 14*number_of_chunks];
		band_dh[k] = d_band_data_ptr[idx + 21*number_of_chunks];
	}
	idx += 21*number_of_chunks;

	short int band_qmf_tx_delayx[24], band_qmf_rx_delayx[24];
	for(int k=0; k<24; k++, idx+=number_of_chunks)
	{
		band_qmf_tx_delayx[k] = d_band_data_ptr[idx];
		band_qmf_rx_delayx[k] = d_band_data_ptr[idx + 24*number_of_chunks];
	}
	idx += 24*number_of_chunks;

	short int xin1;
	short int xin0;
	short int xl;
	short int il;
	short int xh;
	short int ih;

	idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int idx1 = threadIdx.x + blockDim.x * blockIdx.x;
	for (int k=0; k<160; k++, idx += number_of_chunks, idx1 +=(2*number_of_chunks))
	{
		xin1 = d_pcm_data_ptr[idx1];
		xin0 = d_pcm_data_ptr[idx1 + number_of_chunks];

		// Calculation of the synthesis QMF samples 
		// qmf_tx (xin0, xin1, &xl, &xh, encoder);
		CudaQmfTx(xin0, xin1, xl, xh, band_qmf_tx_delayx);

		// Call the upper and lower band ADPCM encoders
		// il = lsbcod (xl, 0, encoder);
		il = CudaLsbCod(xl, band_dlt, band_plt, band_rlt, band_sl, band_detl, band_nbl, band_szl, band_bl, band_al, band_spl);
		// ih = hsbcod (xh, 0, encoder);
		ih = CudaHsbCod(xh, band_dh, band_ph, band_rh, band_sh, band_deth, band_nbh, band_szh, band_bh, band_ah, band_sph);

		// Mount the output G722 codeword: bits 0 to 5 are the lower-band
		// portion of the encoding, and bits 6 and 7 are the upper-band
		// portion of the encoding 
		// code[i] = s_and(add(shl(ih, 6), il), 0xFF);
		d_g722_data_ptr[idx] = (unsigned char) CudaSaturateAddShort(CudaShiftLeftShort(ih, 6), il);

//		d_g722_data_ptr[idx] = (unsigned char) ((xl + xh) / 2);
	}

	// copy local variables back to band data
	idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	// detl	
	d_band_data_ptr[idx] = band_detl;
	idx += number_of_chunks;
	// deth
	d_band_data_ptr[idx] = band_deth;
	idx += number_of_chunks;
	// nbl
	d_band_data_ptr[idx]= band_nbl;
	idx += number_of_chunks;
	// sl	
	d_band_data_ptr[idx] = band_sl;
	idx += number_of_chunks;
	// spl
	d_band_data_ptr[idx]= band_spl;
	idx += number_of_chunks;
	// szl
	d_band_data_ptr[idx] = band_szl;
	idx += number_of_chunks;
	// nbh
	d_band_data_ptr[idx] = band_nbh;
	idx += number_of_chunks;
	// sh	
	d_band_data_ptr[idx] = band_sh;
	idx += number_of_chunks;
	// sph
	d_band_data_ptr[idx]= band_sph;
	idx += number_of_chunks;
	// szh
	d_band_data_ptr[idx] = band_szh;
	idx += number_of_chunks;

	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
		d_band_data_ptr[idx] = band_al[k];
		d_band_data_ptr[idx + 3*number_of_chunks] = band_plt[k];
		d_band_data_ptr[idx + 6*number_of_chunks] = band_rlt[k];
		d_band_data_ptr[idx + 9*number_of_chunks] = band_ah[k];
		d_band_data_ptr[idx + 6*number_of_chunks] = band_ph[k];
		d_band_data_ptr[idx + 9*number_of_chunks] = band_rh[k];
	}
	idx += 15*number_of_chunks;

	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
		d_band_data_ptr[idx] = band_bl[k];
		d_band_data_ptr[idx + 7*number_of_chunks] = band_dlt[k];
		d_band_data_ptr[idx + 14*number_of_chunks] = band_bh[k];
		d_band_data_ptr[idx + 21*number_of_chunks] = band_dh[k];
	}
	idx += 21*number_of_chunks;

	for(int k=0; k<24; k++, idx+=number_of_chunks)
	{
		d_band_data_ptr[idx] = band_qmf_tx_delayx[k];
		d_band_data_ptr[idx + 24*number_of_chunks] = band_qmf_rx_delayx[k];
	}
	idx += 24*number_of_chunks;
}

__global__ void CudaKernelSinusoidalSynthesis(const float* parameter_ptr, int max_sin, int number_of_packets, int samples_per_packet, short int* output_ptr)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (idx >= number_of_packets)
        return;

    int size_of_each_row = 6 * number_of_packets;

    float c = 0;
    float d = 0;


        // initialize output ptr
        for (int l = 1; l <= samples_per_packet; l++)
            output_ptr[(l-1)*number_of_packets + idx] = 0;

        // each packet has max_sin sinusoidals to process
        for (int i = 0; i < max_sin; i++)
        {
            // get necessary data
            int offset = i*size_of_each_row + idx*6;
            float Ap = parameter_ptr[offset];
            float  An = parameter_ptr[offset+1];
            float Op = parameter_ptr[offset+2];
            float On = parameter_ptr[offset+3];
            float Wn = parameter_ptr[offset+4];
            float Wp = parameter_ptr[offset+5];
            float epsilon = 0;

            // classify operation
            if (Ap && An)    // interpolation
            {
                // calculate epsilon
                float total = Op - On + (Wp + Wn)*(samples_per_packet/2);
                float var_n = -2*(PI) - total;
                float var_o = -total;
                float var_p = 2*(PI) - total;
                if ((var_n <= var_o) && ((var_o <= var_p) || (var_n <= var_p)))
                {
                    // search through negative side
                    epsilon = var_n;
                    for (int k = 2; k < 100; k++)
                    {
                        if ((-2*PI*k - total) >= epsilon)
                            break;
                        else
                            epsilon = (-2*PI*k - total);
                    }
                }
                else if (((var_n >= var_o) || (var_n >= var_p)) && (var_o >= var_p))
                {
                    // search through positive side
                    epsilon = var_p;
                    for (int k = 2; k < 100; k++)
                    {
                        if (2*PI*k - total >= epsilon)
                            break;
                        else
                            epsilon = 2*PI*k - total;
                    }
              
                }
                else
                {
                    // ((var_n >= var_o) && (var_o <= var_p))
                    epsilon = var_o;
                }

                // calculate coefficients
                c = (Wn - Wp)/(2*samples_per_packet) + 3*epsilon/(samples_per_packet * samples_per_packet);
                d = -2*epsilon/(samples_per_packet*samples_per_packet*samples_per_packet);
                
            }
            else if (An)    // birth
            {
                Wp = Wn;
                Op = On - (Wp * samples_per_packet);
            }
            else if (Ap)    // death
            {
                
            }
            else    // null
            {
                continue;
            }
           
            for (int l = 1; l <= samples_per_packet; l++)
                output_ptr[(l-1)*number_of_packets + idx] += (((An - Ap)*l/samples_per_packet + Ap)*cos(Op + Wp*l + c*l*l + d*l*l*l))*32000;
        }

}


int CudaG711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);

	unsigned int size_of_alaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_alaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);

	unsigned char* d_alaw_data_ptr = NULL;
	cudaMalloc((void**)&d_alaw_data_ptr, size_of_d_alaw_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_alaw_data_ptr, alaw_data_ptr, size_of_alaw_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());

	// launch kernel here
	CudaKernelG711aToPcmCM <<< grid_dim, block_dim >>> (d_alaw_data_ptr, d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_alaw_data_ptr);
	checkCudaErrors(cudaGetLastError());
   
	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

int CudaG711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);

	unsigned int size_of_ulaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_ulaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);

	unsigned char* d_ulaw_data_ptr = NULL;
	cudaMalloc((void**)&d_ulaw_data_ptr, size_of_d_ulaw_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_ulaw_data_ptr, ulaw_data_ptr, size_of_ulaw_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());

	// launch kernel here
	CudaKernelG711uToPcmCM <<< grid_dim, block_dim >>> (d_ulaw_data_ptr, d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_ulaw_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

int CudaG722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);

	unsigned int size_of_g722_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_g722_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);

	unsigned char* d_g722_data_ptr = NULL;
	cudaMalloc((void**)&d_g722_data_ptr, size_of_d_g722_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_g722_data_ptr, g722_data_ptr, size_of_g722_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());

	// calculate space for band, 45 integers per thread
	unsigned int number_of_d_band_data = grid_dim.x * block_dim.x;
	unsigned int size_of_d_band_data = number_of_d_band_data * 45 * sizeof(int);
	unsigned int number_of_band_data = no_of_data/160;
	unsigned int size_of_band_data = number_of_band_data * 45 * sizeof(int);
	
	//std::cout << "size of band data " << size_of_d_band_data << std::endl;
	
	int* d_band_data_ptr = NULL;
	cudaMalloc((void**)&d_band_data_ptr, size_of_d_band_data);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_band_data_ptr, band_data_ptr, size_of_band_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	

	unsigned int size_of_d_g722_consts = sizeof(g722_consts);
	int* d_g722_consts_ptr = NULL;
	cudaMalloc((void**)&d_g722_consts_ptr, size_of_d_g722_consts);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_g722_consts_ptr, g722_consts, size_of_d_g722_consts, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	
 

	// launch kernel here
	CudaKernelG722ToPcmCM <<< grid_dim, block_dim >>> (d_g722_data_ptr, d_pcm_data_ptr, d_band_data_ptr, d_g722_consts_ptr, no_of_data);
	checkCudaErrors(cudaGetLastError());

	//std::cout << "host pcm data size   : " << size_of_pcm_data << std::endl;
	//std::cout << "device pcm data size : " << size_of_d_pcm_data << std::endl;
	//std::cout << "no of data           : " << no_of_data << std::endl;
	//std::cout << "number of threads    : " << no_of_d_data << std::endl;

	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(band_data_ptr, d_band_data_ptr, size_of_band_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_g722_consts_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_band_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_g722_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

int CudaPcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);
	
	unsigned int size_of_alaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);
	
	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_alaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);
	
	unsigned char* d_alaw_data_ptr = NULL;
	cudaMalloc((void**)&d_alaw_data_ptr, size_of_d_alaw_data);
	checkCudaErrors(cudaGetLastError());
	
	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());
	
	cudaMemcpy(d_pcm_data_ptr, pcm_data_ptr, size_of_pcm_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(d_alaw_data_ptr, alaw_data_ptr, size_of_alaw_data, cudaMemcpyHostToDevice);
//	checkCudaErrors(cudaGetLastError());
	
	// launch kernel here
	CudaKernelPcmToG711aCM <<< grid_dim, block_dim >>> (d_pcm_data_ptr, d_alaw_data_ptr);
	checkCudaErrors(cudaGetLastError());
//	CudaKernelG711aToPcmCM <<< grid_dim, block_dim >>> (d_alaw_data_ptr, d_pcm_data_ptr);
// 	checkCudaErrors(cudaGetLastError());
	
	cudaMemcpy(alaw_data_ptr, d_alaw_data_ptr, size_of_alaw_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
//	checkCudaErrors(cudaGetLastError());
	
	cudaFree(d_alaw_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

int CudaPcmToG711u(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* ulaw_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);

	unsigned int size_of_ulaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_ulaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);

	unsigned char* d_ulaw_data_ptr = NULL;
	cudaMalloc((void**)&d_ulaw_data_ptr, size_of_d_ulaw_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_pcm_data_ptr, pcm_data_ptr, size_of_pcm_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(d_ulaw_data_ptr, ulaw_data_ptr, size_of_ulaw_data, cudaMemcpyHostToDevice);
//	checkCudaErrors(cudaGetLastError());

	// launch kernel here
	CudaKernelPcmToG711uCM <<< grid_dim, block_dim >>> (d_pcm_data_ptr, d_ulaw_data_ptr);
	checkCudaErrors(cudaGetLastError());
//	CudaKernelG711uToPcmCM <<< grid_dim, block_dim >>> (d_ulaw_data_ptr, d_pcm_data_ptr);
//	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(ulaw_data_ptr, d_ulaw_data_ptr, size_of_ulaw_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
//	checkCudaErrors(cudaGetLastError());

	cudaFree(d_ulaw_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

int CudaPcmToG722(const short int* pcm_data_ptr, short int* band_data_ptr, unsigned int no_of_data, unsigned char* g722_data_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 320))), 1, 1);

	unsigned int size_of_g722_data = (no_of_data / 2)* sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 320;
	unsigned int size_of_d_g722_data = (no_of_d_data / 2) * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);

	unsigned char* d_g722_data_ptr = NULL;
	cudaMalloc((void**)&d_g722_data_ptr, size_of_d_g722_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_pcm_data_ptr = NULL;
	cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_pcm_data_ptr, pcm_data_ptr, size_of_pcm_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(d_g722_data_ptr, g722_data_ptr, size_of_g722_data, cudaMemcpyHostToDevice);
//	checkCudaErrors(cudaGetLastError());

	// calculate space for band, 104 short integers per thread
	unsigned int number_of_d_band_data = grid_dim.x * block_dim.x;
	unsigned int size_of_d_band_data = number_of_d_band_data * 104 * sizeof(short int);
	unsigned int number_of_band_data = no_of_data/360;
	unsigned int size_of_band_data = number_of_band_data * 104 * sizeof(short int);
	
	//std::cout << "size of band data " << size_of_d_band_data << std::endl;
	
	short int* d_band_data_ptr = NULL;
	cudaMalloc((void**)&d_band_data_ptr, size_of_d_band_data);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_band_data_ptr, band_data_ptr, size_of_band_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	

	unsigned int size_of_d_g722_consts = sizeof(full_g722_consts);
	short int* d_g722_consts_ptr = NULL;
	cudaMalloc((void**)&d_g722_consts_ptr, size_of_d_g722_consts);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_g722_consts_ptr, full_g722_consts, size_of_d_g722_consts, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	
 

	// launch kernel here
	CudaKernelPcmToG722CM <<< grid_dim, block_dim >>> (d_pcm_data_ptr, d_g722_data_ptr, d_band_data_ptr, d_g722_consts_ptr, no_of_data);
	checkCudaErrors(cudaGetLastError());

	//std::cout << "host pcm data size   : " << size_of_pcm_data << std::endl;
	//std::cout << "device pcm data size : " << size_of_d_pcm_data << std::endl;
	//std::cout << "no of data           : " << no_of_data << std::endl;
	//std::cout << "number of threads    : " << no_of_d_data << std::endl;


	cudaMemcpy(g722_data_ptr, d_g722_data_ptr, size_of_g722_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());
//	cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
//	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(band_data_ptr, d_band_data_ptr, size_of_band_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_g722_consts_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_band_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_g722_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_pcm_data_ptr);
	checkCudaErrors(cudaGetLastError());	

	return 0;
}

void CudaSinusoidalSynthesis(const float* parameter_ptr, int max_sin, int number_of_packets, int samples_per_packet, short int* output_ptr)
{
	dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(number_of_packets/((float)(block_dim.x))), 1, 1);

	unsigned int size_of_parameter_ptr = (number_of_packets * max_sin * 6)* sizeof(float);
	unsigned int size_of_output_ptr = number_of_packets * samples_per_packet * sizeof(short int);

	unsigned int no_of_d_data = grid_dim.x * block_dim.x;
	unsigned int size_of_d_parameter_data = (no_of_d_data ) * sizeof(float) * max_sin * 6;
	unsigned int size_of_d_output_data = no_of_d_data * sizeof(short int) * samples_per_packet;

	float* d_parameter_data_ptr = NULL;
	cudaMalloc((void**)&d_parameter_data_ptr, size_of_d_parameter_data);
	checkCudaErrors(cudaGetLastError());

	short int* d_output_data_ptr = NULL;
	cudaMalloc((void**)&d_output_data_ptr, size_of_d_output_data);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(d_parameter_data_ptr, parameter_ptr, size_of_parameter_ptr, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());

	// launch kernel here
	CudaKernelSinusoidalSynthesis <<< grid_dim, block_dim >>> (d_parameter_data_ptr, max_sin, number_of_packets, samples_per_packet, d_output_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaMemcpy(output_ptr, d_output_data_ptr, size_of_output_ptr, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_parameter_data_ptr);
	checkCudaErrors(cudaGetLastError());

	cudaFree(d_output_data_ptr);
	checkCudaErrors(cudaGetLastError());
}

void CudaGpuInitialize()
{
    unsigned int size_of_d_dummy_data = 1000000;
    unsigned char* d_dummy_data_ptr = NULL;

    //unsigned int device_count = cudaGetDeviceCount();
    //checkCudaErrors(cudaGetLastError()); 
    
    //std::cout << "device count : " << device_count << std::endl
	
    cudaMalloc((void**)&d_dummy_data_ptr, size_of_d_dummy_data);
    checkCudaErrors(cudaGetLastError());
	
    cudaMemset((void*)d_dummy_data_ptr, 0, size_of_d_dummy_data);
    checkCudaErrors(cudaGetLastError());
	
    cudaFree(d_dummy_data_ptr);
    checkCudaErrors(cudaGetLastError());
}

bool DetectCudaDevice()
{
    int number_of_cuda_devices = 0;
    std::vector<cudaDeviceProp> cuda_device_list;

    cudaGetDeviceCount(&number_of_cuda_devices);

    // reset the device list
    cuda_device_list.clear();
    for (int i = 0; i < number_of_cuda_devices; i++)
    {
        cudaDeviceProp cuda_device;
        cudaGetDeviceProperties(&cuda_device, i);
        cuda_device_list.push_back(cuda_device);

        std::cout << "Device : " << i << " " << cuda_device.name << std::endl;
        std::cout << "Compute number of device  : " << cuda_device.major << "." << cuda_device.minor << std::endl;
        std::cout << "Concurrent kernels        : " << cuda_device.concurrentKernels << std::endl;
        std::cout << "Number of sm              : " << cuda_device.multiProcessorCount << std::endl;
        std::cout << "Maximum threads per block : " << cuda_device.maxThreadsPerBlock << std::endl;
        std::cout << "Maximum threads per sm    : " << cuda_device.maxThreadsPerMultiProcessor << std::endl;
        std::cout << "Memory pitch              : " << cuda_device.memPitch << std::endl;
        std::cout << "Registers per block       : " << cuda_device.regsPerBlock << std::endl;
        std::cout << "Shared memory per block   : " << cuda_device.sharedMemPerBlock << std::endl;
        std::cout << std::endl;
    }

    return ((number_of_cuda_devices) ? true : false);
}
