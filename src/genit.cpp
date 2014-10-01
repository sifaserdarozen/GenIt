#include <iostream>
#include <cstdio>
#include <ctime>
#include <ctime>
#include <cstdlib>
#include <cstring>    // for g++ std::memcpy
#include <cmath>
#include <iomanip>    // for std::setprecision
#include <vector>
#include <fstream>

#include "sample_data.h"
#include "cuda_coders.h"
#include "cpu_coders.h"

#define NO_OF_TEST_PACKET 1000000

#define DEFAULT_START_OF_ITERATIONS 250
#define DEFAULT_END_OF_ITERATIONS 128000
#define DEFAULT_NUMBER_OF_ITERATIONS 10
#define DEFAULT_NUMBER_OF_AVERAGING 100

#define MAX_SIN_8KHZ 25
#define MAX_SIN_16KHZ 50

enum DataType {
    NAIVE,
    COALESCED
};

int GenerateG711aData(unsigned char* alaw_data_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!alaw_data_ptr)
        return -1;
		
    unsigned char* dst_ptr = alaw_data_ptr;
    const unsigned char* src_ptr = NULL;
    
    unsigned int no_of_source_data = sizeof(sample_g711a_data) / sizeof(sample_g711a_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711a_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
//          *dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;
        } 	
	
        // copy remeinder datas
        src_ptr = sample_g711a_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
//          *dst_ptr++ = *src_ptr++;	
            *dst_ptr++ = rand() % 256;
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g711u_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
    	
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }				
    }
			
    return 0;
}

int GenerateG711uData(unsigned char* ulaw_data_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!ulaw_data_ptr)
        return -1;
		
    unsigned char* dst_ptr = ulaw_data_ptr;
    const unsigned char* src_ptr = NULL;
    
    unsigned int no_of_source_data = sizeof(sample_g711u_data) / sizeof(sample_g711u_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711u_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
        }
	
        // copy remeinder datas
        src_ptr = sample_g711a_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
            //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;			
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g711u_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }				
    }
	
    return 0;
}

int GenerateG722Data(unsigned char* g722_data_ptr, int* g722_band_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!g722_data_ptr)
        return -1;

    unsigned char* dst_ptr = g722_data_ptr;
    const unsigned char* src_ptr = NULL;

    unsigned int no_of_source_data = sizeof(sample_g722_data) / sizeof(sample_g722_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;
    unsigned int no_of_total_repetitions = ceil(no_of_data/((float)no_of_source_data));
    unsigned int band_size_of_total_repetitions = 45 * sizeof(int) * no_of_total_repetitions;

    std::memset(g722_band_ptr, 0, band_size_of_total_repetitions);

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711u_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
                //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;
        }	

        // copy remeinder datas
        src_ptr = sample_g722_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
            //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;

        for(unsigned int k=0; k<no_of_total_repetitions; k++)
            g722_band_ptr[k*45 + 44] = 32;		
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g722_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }

         for(unsigned int k=0; k<no_of_total_repetitions; k++)
            g722_band_ptr[k + 44*no_of_total_repetitions] = 32;				
    }
	
    return 0;
}

int GeneratePcmData(short int* pcm_8khz_data_ptr, short int* pcm_16khz_data_ptr, short int* g722_full_band_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
	if ((!pcm_8khz_data_ptr) || (!pcm_16khz_data_ptr))
		return -1;
		
	for(unsigned int j=0; j<no_of_data; j++)
	{
		*pcm_16khz_data_ptr++ = *pcm_8khz_data_ptr++ = (rand() % ((unsigned int) 0xffff)) - ((unsigned int) 0x8000);
		*pcm_16khz_data_ptr++ = (rand() % ((unsigned int) 0xffff)) - ((unsigned int) 0x8000);
	}

	std::memset(g722_full_band_ptr, 0, (no_of_data/160)*104*2);
	for(unsigned int k=0; k<(no_of_data / 160); k++)
	{
		g722_full_band_ptr[k] = 32;
		g722_full_band_ptr[k + (no_of_data/160)] = 8;
	}	

	return 0;
}

void GenerateParameterData(float* parameter_ptr, int no_of_test_packet, int max_sin)
{
//		*pcm_16khz_data_ptr++ = *pcm_8khz_data_ptr++ = (rand() % ((unsigned int) 0xffff)) - ((unsigned int) 0x8000);
//		*pcm_16khz_data_ptr++ = (rand() % ((unsigned int) 0xffff)) - ((unsigned int) 0x8000);
    int row_size_of_each_packet = 6;
    int size_of_each_row = row_size_of_each_packet * no_of_test_packet;

    // there are number_of_packets packets to process
    for (int k = 0; k < no_of_test_packet; k++)
    {
        // each packet has max_sin sinusoidals to process
        for (int i = 0; i < max_sin; i++)
        {
            // get necessary data
            int offset = i*size_of_each_row + k*row_size_of_each_packet;
            parameter_ptr[offset] = (float((rand() % 100)) / 100) / max_sin;
            parameter_ptr[offset+1] = (float((rand() % 100)) / 1000) * (parameter_ptr[offset]);
            parameter_ptr[offset+2] = (float((rand() % 100)) / 100) * 2 * PI;
            parameter_ptr[offset+3] = (float((rand() % 100)) / 100) * 2 * PI;
            parameter_ptr[offset+4] = ((i - (float((rand() % 10)) / 100)) / max_sin) * PI;
            parameter_ptr[offset+5] = (1 - (float((rand() % 100)) / 1000)) * parameter_ptr[offset+4];
        }
    }
}

template <class T>
bool CompareArrays(const T* first_data_ptr, const T* second_data_ptr, unsigned int no_of_data)
{
	for(unsigned int k=0; k<no_of_data; k++)
	{
		if(*first_data_ptr++ != *second_data_ptr++)
		{
			std::cout << "match failed at index " << k << std::endl;
			std::cout << (int) (*(--first_data_ptr)) << " - " << (int) *(--second_data_ptr) << std::endl << std::endl;

			return false;
		}			
	}

	std::cout << "match succedded..." << std::endl;
	return true;
}

template <class T>
bool CompareArrays(const short int* pcm_8khz_data_ptr, const T* first_data_ptr, const T* second_data_ptr, unsigned int no_of_data)
{
	const T* org_first_data_ptr = first_data_ptr;
	const T* org_second_data_ptr = second_data_ptr;

	for(unsigned int k=0; k<no_of_data; k++)
	{
		if(*first_data_ptr++ != *second_data_ptr++)
		{
			std::cout << "match failed at index " << k << std::endl;
			std::cout << (int)pcm_8khz_data_ptr[k] << "-->" << (int) (*(--first_data_ptr)) << " - " << (int) *(--second_data_ptr) << std::endl << std::endl;

			for (unsigned int j=0; j<300; j++)
			{
				std::cout << (int)pcm_8khz_data_ptr[j] << "-->" << (int)org_first_data_ptr[j] << "-" << (int)org_second_data_ptr[j] << ((org_first_data_ptr[j] != org_second_data_ptr[j]) ? (" ... fail") : (" ")) << std::endl;
			}
			std::cout << std::endl << std::endl;
	
			return false;
		}			
	}

	std::cout << "match succedded..." << std::endl;
	return true;
}

int main(int argc, char *argv[])
{
    std::vector<int> iteration_vector;
    std::vector<int>::iterator itr;

    std::vector<float> cpu_alaw_vector;
    std::vector<float> gpu_alaw_vector;
    std::vector<float> cpu_ulaw_vector;
    std::vector<float> gpu_ulaw_vector;
    std::vector<float> cpu_g722_vector;
    std::vector<float> gpu_g722_vector;

    std::vector<float> cpu_sinusoidal_8khz_vector;
    std::vector<float> gpu_sinusoidal_8khz_vector;
    std::vector<float> cpu_sinusoidal_16khz_vector;
    std::vector<float> gpu_sinusoidal_16khz_vector;

    std::vector<float>::iterator flt_itr;

    int start_of_iterations = DEFAULT_START_OF_ITERATIONS;
    int end_of_iterations = DEFAULT_END_OF_ITERATIONS;
    int number_of_iterations = DEFAULT_NUMBER_OF_ITERATIONS;
    int number_of_averaging = DEFAULT_NUMBER_OF_AVERAGING;

    // process input arguments and generate necessary data
    if (3 == argc)
    {
        start_of_iterations = atoi(argv[1]);
        number_of_iterations = atoi(argv[2]);
    }
    else if (4==argc)
    {
        start_of_iterations = atoi(argv[1]);
        number_of_iterations = atoi(argv[2]);
        number_of_averaging = atoi(argv[3]);
    }
    else
    {
        // values are already set, print for usage for customization
        std::cout << "calculation will be done by defaut values, to customize, for example" << std::endl; 
        std::cout << "in order to calculate 13 iterations from 250 to 1024000 use..." << std::endl;
        std::cout << "--->   test 250 13" << std::endl << std::endl;
    }

    std::cout << "simulation will run with following parameters" << std::endl;
    std::cout << "start of iterations : " << start_of_iterations << "  " 
              << "end of iterations: " << end_of_iterations << "  "
              << "number of iterations : " << number_of_iterations << std::endl << std::endl;

    DetectCudaDevice();

    iteration_vector.resize(number_of_iterations);
//    float step_size = (end_of_iterations - start_of_iterations) / ((float)(number_of_iterations - 1));
    int current_iteration = start_of_iterations;

    for (itr = iteration_vector.begin(); itr != iteration_vector.end(); current_iteration = current_iteration << 1)
         *itr++ = std::floor(current_iteration + 0.5);

    std::ofstream out_file;
    out_file.open("iteration_results.txt"); 

    float avg_cpu_alaw_time = 0;
    float avg_gpu_alaw_time = 0;
    float avg_cpu_ulaw_time = 0;
    float avg_gpu_ulaw_time = 0;
    float avg_cpu_g722_time = 0;
    float avg_gpu_g722_time = 0;

    float avg_cpu_sinusoidal_8khz_time = 0;
    float avg_gpu_sinusoidal_8khz_time = 0;
    float avg_cpu_sinusoidal_16khz_time = 0;
    float avg_gpu_sinusoidal_16khz_time = 0;

    cpu_alaw_vector.resize(number_of_averaging);
    gpu_alaw_vector.resize(number_of_averaging);
    cpu_ulaw_vector.resize(number_of_averaging);
    gpu_ulaw_vector.resize(number_of_averaging);
    cpu_g722_vector.resize(number_of_averaging);
    gpu_g722_vector.resize(number_of_averaging);

    cpu_sinusoidal_8khz_vector.resize(number_of_averaging);
    gpu_sinusoidal_8khz_vector.resize(number_of_averaging);
    cpu_sinusoidal_16khz_vector.resize(number_of_averaging);
    gpu_sinusoidal_16khz_vector.resize(number_of_averaging);


    unsigned int no_of_test_packet = end_of_iterations;
    unsigned int no_of_test_data = no_of_test_packet * 160;
    short int* pcm_8khz_data_ptr = NULL;
    short int* pcm_16khz_data_ptr = NULL;

    float* parameter_8khz_data_ptr = new float[no_of_test_packet*6*MAX_SIN_8KHZ];
    float* parameter_16khz_data_ptr = new float[no_of_test_packet*6*MAX_SIN_16KHZ]; 
    short int* cpu_sinusoidal_8khz_data_ptr = new short int[no_of_test_data];
    short int* gpu_sinusoidal_8khz_data_ptr = new short int[no_of_test_data];
    short int* cpu_sinusoidal_16khz_data_ptr = new short int[no_of_test_data*2];
    short int* gpu_sinusoidal_16khz_data_ptr = new short int[no_of_test_data*2];

    unsigned char* cpu_encoded_alaw_data_ptr = NULL;
    unsigned char* cpu_encoded_ulaw_data_ptr = NULL;
    unsigned char* cpu_encoded_g722_data_ptr = NULL;	
    unsigned char* gpu_encoded_alaw_data_ptr = NULL;
    unsigned char* gpu_encoded_ulaw_data_ptr = NULL;
    unsigned char* gpu_encoded_g722_data_ptr = NULL;

    unsigned int size_of_g722_band = no_of_test_packet * 104;
    short int* cpu_band_ptr = NULL;
    short int* gpu_band_ptr = NULL;
    short int* g722_band_ptr = NULL;

    pcm_8khz_data_ptr = new short int[no_of_test_data];
    pcm_16khz_data_ptr = new short int[no_of_test_data*2];	
	
    cpu_encoded_alaw_data_ptr = new unsigned char[no_of_test_data];
    cpu_encoded_ulaw_data_ptr = new unsigned char[no_of_test_data];
    cpu_encoded_g722_data_ptr = new unsigned char[no_of_test_data];	
    gpu_encoded_alaw_data_ptr = new unsigned char[no_of_test_data];
    gpu_encoded_ulaw_data_ptr = new unsigned char[no_of_test_data];
    gpu_encoded_g722_data_ptr = new unsigned char[no_of_test_data];
	
    cpu_band_ptr = new short int[size_of_g722_band];
    gpu_band_ptr = new short int[size_of_g722_band];
    g722_band_ptr = new short int[size_of_g722_band];
    //std::memset(cpu_band_ptr, 0, size_of_g722_band * sizeof(int));
    //std::memset(gpu_band_ptr, 0, size_of_g722_band * sizeof(int));
	
    srand(time(NULL));

    // initialize gpu
    CudaGpuInitialize();


    std::cout << "iter." << '\t' 
              << "8khz C" << '\t' << "8khz G" << '\t' << "speedup" << '\t'
              << "16khz C" << '\t' << "16khz G" << '\t' << "speedup" << '\t'
              << "g711a C" << '\t' << "g711a G" << '\t' << "speedup" << '\t' 
              << "g711u C" << '\t' << "g711u G" << '\t' << "speedup" << '\t'
              << "g722 C" << '\t' << "g722 G" << '\t' << "speedup" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    out_file << "iter." << '\t' 
              << "8khz C" << '\t' << "8khz G" << '\t' << "speedup" << '\t'
              << "16khz C" << '\t' << "16khz G" << '\t' << "speedup" << '\t'
              << "g711a C" << '\t' << "g711a G" << '\t' << "speedup" << '\t' 
              << "g711u C" << '\t' << "g711u G" << '\t' << "speedup" << '\t'
              << "g722 C" << '\t' << "g722 G" << '\t' << "speedup" << std::endl;
    out_file << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;


    // loop through iteration vector
    std::cout << std::fixed << std::setprecision(5);
    for (itr = iteration_vector.begin(); itr != iteration_vector.end();)
    {	
        no_of_test_packet = *itr++;
        no_of_test_data = no_of_test_packet * 160;
        size_of_g722_band = no_of_test_packet * 104;

        avg_cpu_alaw_time = 0;
        avg_gpu_alaw_time = 0;
        avg_cpu_ulaw_time = 0;
        avg_gpu_ulaw_time = 0;
        avg_cpu_g722_time = 0;
        avg_gpu_g722_time = 0;
        avg_cpu_sinusoidal_8khz_time = 0;
        avg_gpu_sinusoidal_8khz_time = 0;
        avg_cpu_sinusoidal_16khz_time = 0;
        avg_gpu_sinusoidal_16khz_time = 0;
		
        float cpu_ulaw_time = 0;
        float cpu_alaw_time = 0;
        float cpu_g722_time = 0;
        float cpu_sinusoidal_8khz_time = 0;
        float gpu_sinusoidal_8khz_time = 0;
        float cpu_sinusoidal_16khz_time = 0;
        float gpu_sinusoidal_16khz_time = 0;

        for (int k=0; k<number_of_averaging; k++)
        { 
            // generate parameter data for sinusoidal sytnhesis
            GenerateParameterData(parameter_8khz_data_ptr, no_of_test_packet, MAX_SIN_8KHZ);
            GenerateParameterData(parameter_16khz_data_ptr, no_of_test_packet, MAX_SIN_16KHZ);

            // calculate cpu sinusoidals of case 8khz
            std::clock_t cpu_sinusoidal_8khz_start = std::clock();
            SinusoidalSynthesis(parameter_8khz_data_ptr, MAX_SIN_8KHZ, no_of_test_packet, 160, cpu_sinusoidal_8khz_data_ptr);
            std::clock_t cpu_sinusoidal_8khz_stop = std::clock(); 
            cpu_sinusoidal_8khz_time = (cpu_sinusoidal_8khz_stop - cpu_sinusoidal_8khz_start)/((float) CLOCKS_PER_SEC);

            // calculate cpu sinusoidals of case 16khz
            std::clock_t cpu_sinusoidal_16khz_start = std::clock(); 
            SinusoidalSynthesis(parameter_16khz_data_ptr, MAX_SIN_16KHZ, no_of_test_packet, 320, cpu_sinusoidal_16khz_data_ptr);
            std::clock_t cpu_sinusoidal_16khz_stop = std::clock(); 
            cpu_sinusoidal_16khz_time = (cpu_sinusoidal_16khz_stop - cpu_sinusoidal_16khz_start)/((float) CLOCKS_PER_SEC);

            // calculate gpu sinusoidals of case 8khz
            std::clock_t gpu_sinusoidal_8khz_start = std::clock(); 
            CudaSinusoidalSynthesis(parameter_8khz_data_ptr, MAX_SIN_8KHZ, no_of_test_packet, 160, gpu_sinusoidal_8khz_data_ptr);
            std::clock_t gpu_sinusoidal_8khz_stop = std::clock(); 
            gpu_sinusoidal_8khz_time = (gpu_sinusoidal_8khz_stop - gpu_sinusoidal_8khz_start)/((float) CLOCKS_PER_SEC);

            // calculate gpu sinusoidals of case 16khz
            std::clock_t gpu_sinusoidal_16khz_start = std::clock(); 
            CudaSinusoidalSynthesis(parameter_16khz_data_ptr, MAX_SIN_16KHZ, no_of_test_packet, 320, gpu_sinusoidal_16khz_data_ptr);
            std::clock_t gpu_sinusoidal_16khz_stop = std::clock(); 
            gpu_sinusoidal_16khz_time = (gpu_sinusoidal_16khz_stop - gpu_sinusoidal_16khz_start)/((float) CLOCKS_PER_SEC);



        // generate pcm test data
        GeneratePcmData(pcm_8khz_data_ptr, pcm_16khz_data_ptr, g722_band_ptr, no_of_test_data, COALESCED);
        std::memcpy(cpu_band_ptr, g722_band_ptr, size_of_g722_band * sizeof(short int));
        std::memcpy(gpu_band_ptr, g722_band_ptr, size_of_g722_band * sizeof(short int));
	
        // calculate cpu cases once since it takes too much time to complete
        //if (10  > k)
        {
        // encode ulaw using cpu
        std::clock_t cpu_ulaw_start = std::clock(); 
        PcmToG711u(pcm_8khz_data_ptr, no_of_test_data, cpu_encoded_ulaw_data_ptr);
        std::clock_t cpu_ulaw_stop = std::clock(); 
        cpu_ulaw_time = (cpu_ulaw_stop - cpu_ulaw_start)/((float) CLOCKS_PER_SEC);
        //std::cout << cpu_ulaw_time << std::endl;

        // encode alaw using cpu
        std::clock_t cpu_alaw_start = std::clock(); 
        PcmToG711a(pcm_8khz_data_ptr, no_of_test_data, cpu_encoded_alaw_data_ptr);
        std::clock_t cpu_alaw_stop = std::clock();
        cpu_alaw_time = (cpu_alaw_stop - cpu_alaw_start)/((float) CLOCKS_PER_SEC);
        //std::cout << cpu_alaw_time << std::endl;

        // encode g722 using cpu
        std::clock_t cpu_g722_start = std::clock(); 
        PcmToG722(pcm_16khz_data_ptr, cpu_band_ptr, no_of_test_data*2, cpu_encoded_g722_data_ptr);
        std::clock_t cpu_g722_stop = std::clock(); 
        cpu_g722_time = (cpu_g722_stop - cpu_g722_start)/((float) CLOCKS_PER_SEC);
        }

        // encode ulaw using gpu
        std::clock_t gpu_ulaw_start = std::clock(); 
        CudaPcmToG711u(pcm_8khz_data_ptr, no_of_test_data, gpu_encoded_ulaw_data_ptr);
        std::clock_t gpu_ulaw_stop = std::clock(); 
        float gpu_ulaw_time = (gpu_ulaw_stop - gpu_ulaw_start)/((float) CLOCKS_PER_SEC);
	
        // encode alaw using gpu
        std::clock_t gpu_alaw_start = std::clock(); 
        CudaPcmToG711a(pcm_8khz_data_ptr, no_of_test_data, gpu_encoded_alaw_data_ptr);
        std::clock_t gpu_alaw_stop = std::clock();
        float gpu_alaw_time = (gpu_alaw_stop - gpu_alaw_start)/((float) CLOCKS_PER_SEC);
	
        // decode g722 using gpu
        std::clock_t gpu_g722_start = std::clock(); 
        CudaPcmToG722(pcm_16khz_data_ptr, gpu_band_ptr, no_of_test_data*2, gpu_encoded_g722_data_ptr);
        std::clock_t gpu_g722_stop = std::clock(); 
        float gpu_g722_time = (gpu_g722_stop - gpu_g722_start)/((float) CLOCKS_PER_SEC);

        /*out_file  << no_of_test_packet << "\t" 
                  << cpu_alaw_time << "\t" << gpu_alaw_time << "\t"
                  << cpu_ulaw_time << "\t" << gpu_ulaw_time << "\t"
                  << cpu_g722_time << "\t" << gpu_g722_time << std::endl; */

        cpu_alaw_vector[k] = cpu_alaw_time;
        gpu_alaw_vector[k] = gpu_alaw_time;
        cpu_ulaw_vector[k] = cpu_ulaw_time;
        gpu_ulaw_vector[k] = gpu_ulaw_time;
        cpu_g722_vector[k] = cpu_g722_time;
        gpu_g722_vector[k] = gpu_g722_time;
            cpu_sinusoidal_8khz_vector[k] = cpu_sinusoidal_8khz_time;
            cpu_sinusoidal_16khz_vector[k] = cpu_sinusoidal_16khz_time;
            gpu_sinusoidal_8khz_vector[k] = gpu_sinusoidal_8khz_time;
            gpu_sinusoidal_16khz_vector[k] = gpu_sinusoidal_16khz_time;

        avg_cpu_alaw_time += (cpu_alaw_time/number_of_averaging);
        avg_gpu_alaw_time += (gpu_alaw_time/number_of_averaging);
        avg_cpu_ulaw_time += (cpu_ulaw_time/number_of_averaging);
        avg_gpu_ulaw_time += (gpu_ulaw_time/number_of_averaging);
        avg_cpu_g722_time += (cpu_g722_time/number_of_averaging);
        avg_gpu_g722_time += (gpu_g722_time/number_of_averaging);
            avg_cpu_sinusoidal_8khz_time += (cpu_sinusoidal_8khz_time/number_of_averaging);
            avg_cpu_sinusoidal_16khz_time += (cpu_sinusoidal_16khz_time/number_of_averaging);
            avg_gpu_sinusoidal_8khz_time += (gpu_sinusoidal_8khz_time/number_of_averaging);
            avg_gpu_sinusoidal_16khz_time += (gpu_sinusoidal_16khz_time/number_of_averaging);
        }

        // calculate std variance
        float std_cpu_alaw_time = 0;
        float std_gpu_alaw_time = 0;
        float std_cpu_ulaw_time = 0;
        float std_gpu_ulaw_time = 0;
        float std_cpu_g722_time = 0;
        float std_gpu_g722_time = 0;
        float std_cpu_sinusoidal_8khz_time = 0;
        float std_cpu_sinusoidal_16khz_time = 0;
        float std_gpu_sinusoidal_8khz_time = 0;
        float std_gpu_sinusoidal_16khz_time = 0;

        for (int k=0; k<number_of_averaging; k++)
        {
            std_cpu_alaw_time += std::pow((double)((cpu_alaw_vector[k] - avg_cpu_alaw_time)), 2.0);
            std_gpu_alaw_time += std::pow((double)((gpu_alaw_vector[k] - avg_gpu_alaw_time)), 2.0);
            std_cpu_ulaw_time += std::pow((double)((cpu_ulaw_vector[k] - avg_cpu_ulaw_time)), 2.0);
            std_gpu_ulaw_time += std::pow((double)((gpu_ulaw_vector[k] - avg_gpu_ulaw_time)), 2.0);
            std_cpu_g722_time += std::pow((double)((cpu_g722_vector[k] - avg_cpu_g722_time)), 2.0);
            std_gpu_g722_time += std::pow((double)((gpu_g722_vector[k] - avg_gpu_g722_time)), 2.0);
            std_cpu_sinusoidal_8khz_time += std::pow((double)((cpu_sinusoidal_8khz_vector[k] - avg_cpu_sinusoidal_8khz_time)), 2.0);
            std_cpu_sinusoidal_16khz_time += std::pow((double)((cpu_sinusoidal_16khz_vector[k] - avg_cpu_sinusoidal_16khz_time)), 2.0);
            std_gpu_sinusoidal_8khz_time += std::pow((double)((gpu_sinusoidal_8khz_vector[k] - avg_gpu_sinusoidal_8khz_time)), 2.0);
            std_gpu_sinusoidal_16khz_time += std::pow((double)((gpu_sinusoidal_16khz_vector[k] - avg_gpu_sinusoidal_16khz_time)), 2.0); 
        }

        // display results
        std::cout << no_of_test_packet << "\t" 
                  << avg_cpu_sinusoidal_8khz_time << "\t"
                  << avg_gpu_sinusoidal_8khz_time << "\t"
                  << std::setprecision(4) << avg_cpu_sinusoidal_8khz_time/avg_gpu_sinusoidal_8khz_time << std::setprecision(5) << "\t"
                  << avg_cpu_sinusoidal_16khz_time << "\t"
                  << avg_gpu_sinusoidal_16khz_time << "\t"
                  << std::setprecision(4) << avg_cpu_sinusoidal_16khz_time/avg_gpu_sinusoidal_16khz_time << std::setprecision(5) << "\t"
                  << avg_cpu_alaw_time << "\t" 
                  << avg_gpu_alaw_time << "\t" 
//                  << avg_cpu_alaw_time << '(' << std_cpu_alaw_time << ')' << "\t" 
//                  << avg_gpu_alaw_time << '(' << std_gpu_alaw_time << ')' << "\t" 
                  << std::setprecision(4) << avg_cpu_alaw_time/avg_gpu_alaw_time << std::setprecision(5) << "\t"
                  << avg_cpu_ulaw_time << "\t" 
                  << avg_gpu_ulaw_time << "\t" 
//                  << avg_cpu_ulaw_time << '(' << std_cpu_ulaw_time << ')' << "\t" 
//                  << avg_gpu_ulaw_time << '(' << std_gpu_ulaw_time << ')' << "\t" 
                  << std::setprecision(4) << avg_cpu_ulaw_time/avg_gpu_ulaw_time << std::setprecision(5) << "\t"
                  << avg_cpu_g722_time << "\t" 
                  << avg_gpu_g722_time << "\t"	
//                  << avg_cpu_g722_time << '(' << std_cpu_g722_time << ')' << "\t" 
//                  << avg_gpu_g722_time << '(' << std_gpu_g722_time << ')' << "\t"
                  << std::setprecision(4) << avg_cpu_g722_time/avg_gpu_g722_time << std::setprecision(5) << std::endl;
   
        out_file  << no_of_test_packet << "\t"
                  << avg_cpu_sinusoidal_8khz_time << '(' << std_cpu_sinusoidal_8khz_time << ')' << "\t"
                  << avg_gpu_sinusoidal_8khz_time << '(' << std_gpu_sinusoidal_8khz_time << ')' <<"\t"
                  << avg_cpu_sinusoidal_8khz_time/avg_gpu_sinusoidal_8khz_time << "\t"
                  << avg_cpu_sinusoidal_16khz_time << '(' << std_cpu_sinusoidal_16khz_time << ')' << "\t"
                  << avg_gpu_sinusoidal_16khz_time << '(' << std_gpu_sinusoidal_16khz_time << ')' << "\t"
                  << avg_cpu_sinusoidal_16khz_time/avg_gpu_sinusoidal_16khz_time << "\t"
                  << avg_cpu_alaw_time << '(' << std_cpu_alaw_time << ')' << "\t" 
                  << avg_gpu_alaw_time << '(' << std_gpu_alaw_time << ')' << "\t" 
		  << avg_cpu_alaw_time/avg_gpu_alaw_time << "\t"
                  << avg_cpu_ulaw_time << '(' << std_cpu_ulaw_time << ')' << "\t" 
                  << avg_gpu_ulaw_time << '(' << std_gpu_ulaw_time << ')' << "\t" 
		  << avg_cpu_ulaw_time/avg_gpu_ulaw_time << "\t"
                  << avg_cpu_g722_time << '(' << std_cpu_g722_time << ')' << "\t" 
                  << avg_gpu_g722_time << '(' << std_gpu_g722_time << ')' << "\t"
		  << avg_cpu_g722_time/avg_gpu_g722_time << std::endl; 

    }
//	for (unsigned int k=0; k<10; k++)
//		std::cout << pcm_8khz_data_ptr[k] << " ";
//	std::cout << std::endl << std::endl;

    CompareArrays<unsigned char>(cpu_encoded_alaw_data_ptr, gpu_encoded_alaw_data_ptr, no_of_test_data);
    std::cout << std::endl;

    CompareArrays<unsigned char>(cpu_encoded_ulaw_data_ptr, gpu_encoded_ulaw_data_ptr, no_of_test_data);		
    std::cout << std::endl;

    CompareArrays<unsigned char>(cpu_encoded_g722_data_ptr, gpu_encoded_g722_data_ptr, no_of_test_data*2);
    CompareArrays<short int>(cpu_band_ptr, gpu_band_ptr, size_of_g722_band);
    std::cout << std::endl;

    out_file.close();
	

    if(pcm_16khz_data_ptr)
        delete []pcm_16khz_data_ptr;
    if(pcm_8khz_data_ptr)
        delete []pcm_8khz_data_ptr;
    if(cpu_encoded_alaw_data_ptr)
        delete []cpu_encoded_alaw_data_ptr;
    if(cpu_encoded_ulaw_data_ptr)
        delete []cpu_encoded_ulaw_data_ptr;
    if(cpu_encoded_g722_data_ptr)
        delete []cpu_encoded_g722_data_ptr;
    if(cpu_band_ptr)
        delete []cpu_band_ptr;	
    if(gpu_encoded_alaw_data_ptr)
        delete []gpu_encoded_alaw_data_ptr;
    if(gpu_encoded_ulaw_data_ptr)
        delete []gpu_encoded_ulaw_data_ptr;
    if(gpu_encoded_g722_data_ptr)
        delete []gpu_encoded_g722_data_ptr;	
    if(gpu_band_ptr)
        delete []gpu_band_ptr;
    if (parameter_8khz_data_ptr)
        delete []parameter_8khz_data_ptr;
    if (parameter_16khz_data_ptr)
        delete []parameter_16khz_data_ptr;
    if (cpu_sinusoidal_8khz_data_ptr)
        delete []cpu_sinusoidal_8khz_data_ptr;
    if (gpu_sinusoidal_8khz_data_ptr)
        delete []gpu_sinusoidal_8khz_data_ptr;
    if (cpu_sinusoidal_16khz_data_ptr)
        delete []cpu_sinusoidal_16khz_data_ptr;
    if (gpu_sinusoidal_16khz_data_ptr)
        delete []gpu_sinusoidal_16khz_data_ptr;
	
    return 0;
}
