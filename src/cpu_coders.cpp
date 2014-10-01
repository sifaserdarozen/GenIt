#include "cpu_coders.h"
#include "g722coder.h"
#include <iostream>
#include <cmath>

int G711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!alaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	//double normalizing_ratio=((double)(0x9ffff))/((double)(0x1ffff)); // this will map 13 bit pcm to pseudo 16 bit pcm
	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;
//	short int pcm_data;
	for (unsigned int k=0; k<no_of_data; k++)
	{
		alaw_data=*alaw_data_ptr++;
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
		
		//*pcm_vector++=((alaw_data & (0x80))?quantization_value:-quantization_value)*normalizing_ratio;
		*pcm_data_ptr++=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}

	return 0;
}

int PcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!alaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	short int quantization_segment;
	short int quantization_value;

	short int pcm_data;
	unsigned char alaw_data;

	for(unsigned int k=0; k<no_of_data; k++)
	{
		alaw_data = 0;
		pcm_data=*pcm_data_ptr++;
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

		*alaw_data_ptr++=alaw_data;
	}

	return 0;
}

int G711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!ulaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	//double normalizing_ratio=((double)(0x9ffff))/((double)(0x1ffff)); // this will map 13 bit pcm to pseudo 16 bit pcm
	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;
//	short int pcm_data;
	for (unsigned int k=0; k<no_of_data; k++)
	{
		ulaw_data=~(*ulaw_data_ptr++);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);
	
		//*pcm_vector++=((alaw_data & (0x80))?quantization_value:-quantization_value)*normalizing_ratio;
		*pcm_data_ptr++=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}

	return 0;
}

int PcmToG711u(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* ulaw_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!ulaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	short int quantization_segment = 1;
	short int quantization_value = 0;

	short int pcm_data;
	unsigned char ulaw_data;

	for(unsigned int k=0; k<no_of_data; k++)
	{
		pcm_data=*pcm_data_ptr++;
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

		*ulaw_data_ptr++ = ulaw_data;
	}

	return 0;
}

int G722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!g722_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	if (!band_data_ptr)
	{
		std::cerr << "band data is null" << std::endl;
		return -1;
	}
	
	unsigned int number_of_chunks = no_of_data/160;
	
	for(unsigned int k=0; k<number_of_chunks; k++)
	{
	
		unsigned char g722_data_arr[160];
		unsigned short pcm_data_arr[160];
		int band_data[45];
		
		for(int j=0; j<160; j++)
			g722_data_arr[j] = g722_data_ptr[k + number_of_chunks*j];

		for(int j=0; j<45; j++)
			band_data[j]= band_data_ptr[k + number_of_chunks*j];		
			
		G722DecoderType g722_decoder;
		g722_decoder.Decode(g722_data_arr, band_data, pcm_data_arr, 160);
	
		for(int j=0; j<160; j++)
			pcm_data_ptr[k + number_of_chunks*j] = pcm_data_arr[j];

		for(int j=0; j<45; j++)
			band_data_ptr[k+ number_of_chunks*j] = band_data[j];		
	}

	return 0;
}

int PcmToG722(const short int* pcm_data_ptr, short int* band_data_ptr, unsigned int no_of_data, unsigned char* g722_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!g722_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	if (!band_data_ptr)
	{
		std::cerr << "band data is null" << std::endl;
		return -1;
	}

	unsigned int number_of_chunks = no_of_data/320;
	
	for(unsigned int k=0; k<number_of_chunks; k++)
	{
	
		unsigned char g722_data_arr[160];
		short int pcm_data_arr[320];
		short int band_data[104];
		
		for(int j=0; j<320; j++)
			pcm_data_arr[j] = pcm_data_ptr[k + number_of_chunks*j];

		for(int j=0; j<104; j++)
			band_data[j]= band_data_ptr[k + number_of_chunks*j];		
			
		G722EncoderType g722_encoder;
		g722_encoder.Encode(pcm_data_arr, band_data, g722_data_arr, 320);
	
		for(int j=0; j<160; j++)
			g722_data_ptr[k + number_of_chunks*j] = g722_data_arr[j];

		for(int j=0; j<104; j++)
			band_data_ptr[k+ number_of_chunks*j] = band_data[j];		
	}

	return 0;
}

void SinusoidalSynthesis(const float* parameter_ptr, int max_sin, int number_of_packets, int samples_per_packet, short int* output_ptr)
{
    int row_size_of_each_packet = 6;
    int size_of_each_row = row_size_of_each_packet * number_of_packets;

    float c = 0;
    float d = 0;

    // there are number_of_packets packets to process
    for (int k = 0; k < number_of_packets; k++)
    {
        // initialize output ptr
        for (int l = 1; l <= samples_per_packet; l++)
            output_ptr[(l-1)*number_of_packets + k] = 0;

        // each packet has max_sin sinusoidals to process
        for (int i = 0; i < max_sin; i++)
        {
            // get necessary data
            int offset = i*size_of_each_row + k*row_size_of_each_packet;
            float Ap = parameter_ptr[offset];
            float An = parameter_ptr[offset+1];
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
                output_ptr[(l-1)*number_of_packets + k] += (((An - Ap)*l/samples_per_packet + Ap)*cos(Op + Wp*l + c*l*l + d*l*l*l))*32000;
        }
    }


}

