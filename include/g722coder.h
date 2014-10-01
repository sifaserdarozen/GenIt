#ifndef _G722_CODER_TYPE
#define _G722_CODER_TYPE

//const int wl[8] = {-60, -30, 58, 172, 334, 538, 1198, 3042 };
const int rl42[16] = {0, 7, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1, 0 };
const int ilb[32] = {
		2048, 2093, 2139, 2186, 2233, 2282, 2332,
		2383, 2435, 2489, 2543, 2599, 2656, 2714,
		2774, 2834, 2896, 2960, 3025, 3091, 3158,
		3228, 3298, 3371, 3444, 3520, 3597, 3676,
		3756, 3838, 3922, 4008
	};
const int qm4[16] = {
	      0, -20456, -12896,  -8968,
	  -6288,  -4240,  -2584,  -1200,
	  20456,  12896,   8968,   6288,
	   4240,   2584,   1200,      0
	};
const int qm6[64] = {
		-136,   -136,   -136,   -136,
		-24808, -21904, -19008, -16704,
		-14984, -13512, -12280, -11192,
		-10232,  -9360,  -8576,  -7856,
		-7192,  -6576,  -6000,  -5456,
		-4944,  -4464,  -4008,  -3576,
		-3168,  -2776,  -2400,  -2032,
		-1688,  -1360,  -1040,   -728,
		24808,  21904,  19008,  16704,
		14984,  13512,  12280,  11192,
		10232,   9360,   8576,   7856,
		7192,   6576,   6000,   5456,
		4944,   4464,   4008,   3576,
		3168,   2776,   2400,   2032,
		1688,   1360,   1040,    728,
		432,    136,   -432,   -136
	};
	
const int g722_consts[146] = {
	 	// wl[8]
		-60, -30, 58, 172, 334, 538, 1198, 3042,
		// r142[16]
		0, 7, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1, 0,
		// ilb[32]
		2048, 2093, 2139, 2186, 2233, 2282, 2332,
		2383, 2435, 2489, 2543, 2599, 2656, 2714,
		2774, 2834, 2896, 2960, 3025, 3091, 3158,
		3228, 3298, 3371, 3444, 3520, 3597, 3676,
		3756, 3838, 3922, 4008,
		// qm4[16]
		0, -20456, -12896,  -8968,
		-6288,  -4240,  -2584,  -1200,
		20456,  12896,   8968,   6288,
		4240,   2584,   1200,      0,
		// qm6[64]
		-136,   -136,   -136,   -136,
		-24808, -21904, -19008, -16704,
		-14984, -13512, -12280, -11192,
		-10232,  -9360,  -8576,  -7856,
		-7192,  -6576,  -6000,  -5456,
		-4944,  -4464,  -4008,  -3576,
		-3168,  -2776,  -2400,  -2032,
		-1688,  -1360,  -1040,   -728,
		24808,  21904,  19008,  16704,
		14984,  13512,  12280,  11192,
		10232,   9360,   8576,   7856,
		7192,   6576,   6000,   5456,
		4944,   4464,   4008,   3576,
		3168,   2776,   2400,   2032,
		1688,   1360,   1040,    728,
		432,    136,   -432,   -136
		};
	
struct BandType
{
	int s;
	int sp;
	int sz;
	int r[3];
	int a[3];
	int ap[3];
	int p[3];
	int d[7];
	int b[7];
	int bp[7];
	int sg[7];
	int nb;
	int det;
};

class G722DecoderType
{
private:
    BandType band;

protected:

public:
	G722DecoderType();
	~G722DecoderType();
	short int ConvertLongToShort(long int in_value);
	unsigned short Decode(unsigned char g722_data);
	void Decode(unsigned char* g722_data_ptr, int* band_data, unsigned short* pcm_data_ptr, unsigned int no_of_data);
	void Reset();
};

//**********************************************************************************

/* **** Coefficients for both transmission and reception QMF **** */
const short int coef_qmf[24] =
{
	3 * 2, -11 * 2, -11 * 2, 53 * 2, 12 * 2, -156 * 2,
	32 * 2, 362 * 2, -210 * 2, -805 * 2, 951 * 2, 3876 * 2,
	3876 * 2, 951 * 2, -805 * 2, -210 * 2, 362 * 2, 32 * 2,
	-156 * 2, 12 * 2, 53 * 2, -11 * 2, -11 * 2, 3 * 2
};

const short int misil[2][32] =
{
	{
		0x0000, 0x003F, 0x003E, 0x001F, 0x001E, 0x001D, 0x001C, 0x001B,
		0x001A, 0x0019, 0x0018, 0x0017, 0x0016, 0x0015, 0x0014, 0x0013,
		0x0012, 0x0011, 0x0010, 0x000F, 0x000E, 0x000D, 0x000C, 0x000B,
		0x000A, 0x0009, 0x0008, 0x0007, 0x0006, 0x0005, 0x0004, 0x0000
	},
	{
		0x0000, 0x003D, 0x003C, 0x003B, 0x003A, 0x0039, 0x0038, 0x0037,
		0x0036, 0x0035, 0x0034, 0x0033, 0x0032, 0x0031, 0x0030, 0x002F,
		0x002E, 0x002D, 0x002C, 0x002B, 0x002A, 0x0029, 0x0028, 0x0027,
		0x0026, 0x0025, 0x0024, 0x0023, 0x0022, 0x0021, 0x0020, 0x0000
	}
};

const short int q6[31] =
{
	0, 35, 72, 110, 150, 190, 233, 276,
	323, 370, 422, 473, 530, 587, 650, 714,
	786, 858, 940, 1023, 1121, 1219, 1339, 1458,
	1612, 1765, 1980, 2195, 2557, 2919, 3200
};

const short int ril4[16] = {0, 7, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1, 0};
const short int risil[16] = {0, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0};
const short int oq4[8] = {0, 150, 323, 530, 786, 1121, 1612, 2557};
const short int wl[8] = {-60, -30, 58, 172, 334, 538, 1198, 3042};

const short int ila[353] = {
	1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2,
	3, 3, 3, 3, 3, 3, 3, 3,
	3, 3, 3, 4, 4, 4, 4, 4,
	4, 4, 4, 5, 5, 5, 5, 5,
	5, 5, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 8, 8,
	8, 8, 8, 9, 9, 9, 9, 10,
	10, 10, 10, 11, 11, 11, 11, 12,
	12, 12, 13, 13, 13, 13, 14, 14,
	15, 15, 15, 16, 16, 16, 17, 17,
	18, 18, 18, 19, 19, 20, 20, 21,
	21, 22, 22, 23, 23, 24, 24, 25,
	25, 26, 27, 27, 28, 28, 29, 30,
	31, 31, 32, 33, 33, 34, 35, 36,
	37, 37, 38, 39, 40, 41, 42, 43,
	44, 45, 46, 47, 48, 49, 50, 51,
	52, 54, 55, 56, 57, 58, 60, 61,
	63, 64, 65, 67, 68, 70, 71, 73,
	75, 76, 78, 80, 82, 83, 85, 87,
	89, 91, 93, 95, 97, 99, 102, 104,
	106, 109, 111, 113, 116, 118, 121, 124,
	127, 129, 132, 135, 138, 141, 144, 147,
	151, 154, 157, 161, 165, 168, 172, 176,
	180, 184, 188, 192, 196, 200, 205, 209,
	214, 219, 223, 228, 233, 238, 244, 249,
	255, 260, 266, 272, 278, 284, 290, 296,
	303, 310, 316, 323, 331, 338, 345, 353,
	361, 369, 377, 385, 393, 402, 411, 420,
	429, 439, 448, 458, 468, 478, 489, 500,
	511, 522, 533, 545, 557, 569, 582, 594,
	607, 621, 634, 648, 663, 677, 692, 707,
	723, 739, 755, 771, 788, 806, 823, 841,
	860, 879, 898, 918, 938, 958, 979, 1001,
	1023, 1045, 1068, 1092, 1115, 1140, 1165, 1190,
	1216, 1243, 1270, 1298, 1327, 1356, 1386, 1416,
	1447, 1479, 1511, 1544, 1578, 1613, 1648, 1684,
	1721, 1759, 1797, 1837, 1877, 1918, 1960, 2003,
	2047, 2092, 2138, 2185, 2232, 2281, 2331, 2382,
	2434, 2488, 2542, 2598, 2655, 2713, 2773, 2833,
	2895, 2959, 3024, 3090, 3157, 3227, 3297, 3370,
	3443, 3519, 3596, 3675, 3755, 3837, 3921, 4007,
	4095};

const short int misih[2][3] = {{0, 1, 0}, {0, 3, 2}};
const short ih2[4] = {2, 1, 2, 1};
const short sih[4] = {-1, -1, 0, 0};
const short oq2[3] = {0, 202, 926};
const short wh[3] = {0, -214, 798};

const short int full_g722_consts[] = {
	// const short int coef_qmf[24]
	3 * 2, -11 * 2, -11 * 2, 53 * 2, 12 * 2, -156 * 2,
	32 * 2, 362 * 2, -210 * 2, -805 * 2, 951 * 2, 3876 * 2,
	3876 * 2, 951 * 2, -805 * 2, -210 * 2, 362 * 2, 32 * 2,
	-156 * 2, 12 * 2, 53 * 2, -11 * 2, -11 * 2, 3 * 2,

	//const short int misil[2][32]
	0x0000, 0x003F, 0x003E, 0x001F, 0x001E, 0x001D, 0x001C, 0x001B,
	0x001A, 0x0019, 0x0018, 0x0017, 0x0016, 0x0015, 0x0014, 0x0013,
	0x0012, 0x0011, 0x0010, 0x000F, 0x000E, 0x000D, 0x000C, 0x000B,
	0x000A, 0x0009, 0x0008, 0x0007, 0x0006, 0x0005, 0x0004, 0x0000,
	0x0000, 0x003D, 0x003C, 0x003B, 0x003A, 0x0039, 0x0038, 0x0037,
	0x0036, 0x0035, 0x0034, 0x0033, 0x0032, 0x0031, 0x0030, 0x002F,
	0x002E, 0x002D, 0x002C, 0x002B, 0x002A, 0x0029, 0x0028, 0x0027,
	0x0026, 0x0025, 0x0024, 0x0023, 0x0022, 0x0021, 0x0020, 0x0000,

	// const short int q6[31]
	0, 35, 72, 110, 150, 190, 233, 276,
	323, 370, 422, 473, 530, 587, 650, 714,
	786, 858, 940, 1023, 1121, 1219, 1339, 1458,
	1612, 1765, 1980, 2195, 2557, 2919, 3200,

	// const short int ril4[16]
	0, 7, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1, 0,

	// const short int risil[16]
	0, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,
	// const short int oq4[8]
	0, 150, 323, 530, 786, 1121, 1612, 2557,
	// const short int wl[8]
	-60, -30, 58, 172, 334, 538, 1198, 3042,

	// const short int ila[353]
	1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2,
	3, 3, 3, 3, 3, 3, 3, 3,
	3, 3, 3, 4, 4, 4, 4, 4,
	4, 4, 4, 5, 5, 5, 5, 5,
	5, 5, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 8, 8,
	8, 8, 8, 9, 9, 9, 9, 10,
	10, 10, 10, 11, 11, 11, 11, 12,
	12, 12, 13, 13, 13, 13, 14, 14,
	15, 15, 15, 16, 16, 16, 17, 17,
	18, 18, 18, 19, 19, 20, 20, 21,
	21, 22, 22, 23, 23, 24, 24, 25,
	25, 26, 27, 27, 28, 28, 29, 30,
	31, 31, 32, 33, 33, 34, 35, 36,
	37, 37, 38, 39, 40, 41, 42, 43,
	44, 45, 46, 47, 48, 49, 50, 51,
	52, 54, 55, 56, 57, 58, 60, 61,
	63, 64, 65, 67, 68, 70, 71, 73,
	75, 76, 78, 80, 82, 83, 85, 87,
	89, 91, 93, 95, 97, 99, 102, 104,
	106, 109, 111, 113, 116, 118, 121, 124,
	127, 129, 132, 135, 138, 141, 144, 147,
	151, 154, 157, 161, 165, 168, 172, 176,
	180, 184, 188, 192, 196, 200, 205, 209,
	214, 219, 223, 228, 233, 238, 244, 249,
	255, 260, 266, 272, 278, 284, 290, 296,
	303, 310, 316, 323, 331, 338, 345, 353,
	361, 369, 377, 385, 393, 402, 411, 420,
	429, 439, 448, 458, 468, 478, 489, 500,
	511, 522, 533, 545, 557, 569, 582, 594,
	607, 621, 634, 648, 663, 677, 692, 707,
	723, 739, 755, 771, 788, 806, 823, 841,
	860, 879, 898, 918, 938, 958, 979, 1001,
	1023, 1045, 1068, 1092, 1115, 1140, 1165, 1190,
	1216, 1243, 1270, 1298, 1327, 1356, 1386, 1416,
	1447, 1479, 1511, 1544, 1578, 1613, 1648, 1684,
	1721, 1759, 1797, 1837, 1877, 1918, 1960, 2003,
	2047, 2092, 2138, 2185, 2232, 2281, 2331, 2382,
	2434, 2488, 2542, 2598, 2655, 2713, 2773, 2833,
	2895, 2959, 3024, 3090, 3157, 3227, 3297, 3370,
	3443, 3519, 3596, 3675, 3755, 3837, 3921, 4007,
	4095,

	// const short int misih[2][3]
	0, 1, 0,
	0, 3, 2,

	// const short ih2[4]
	2, 1, 2, 1,
	//const short sih[4]
	-1, -1, 0, 0,
	// const short oq2[3]
	0, 202, 926,
	// const short wh[3]
	0, -214, 798
};

#define MIN_32 (int)(0x80000000)
#define MAX_32 (int)(0x7fffffff)
#define MAX_16 (short int)(0x7fff)
#define MIN_16 (short int)(0x8000)

struct FullBandType
{
	short int detl;
	short int deth;
	short int nbl;
	short int sl;
	short int spl;
	short int szl;
	short int nbh;
	short int sh;
	short int sph;
	short int szh;
	short int al[3];
	short int plt[3]; /* plt[0]=plt */
	short int rlt[3];
	short int ah[3];
	short int ph[3]; /* ph[0]=ph */
	short int rh[3];
	short int bl[7];
	short int dlt[7]; /* dlt[0]=dlt */
	short int bh[7];
	short int dh[7]; /* dh[0]=dh */
	short int qmf_tx_delayx[24];
	short int qmf_rx_delayx[24];
};


class G722EncoderType
{
private:
	FullBandType band;

	int SaturateAdd(int op1, int op2) const;
	int SaturateSubtract(int op1, int op2) const;
	short int SaturateSubtractShort(short int op1, short int op2) const;
	short int SaturateAddShort(short int op1, short int op2) const;
	int MultiplyAdd(int add_op, short int mul_op1, short int mul_op2) const;
	int ShiftRight(int op1, short int op2) const;
	int ShiftLeft(int op1, short int op2) const;
	short int ShiftLeftShort(short int op1, short int op2) const;
	short int ShiftRightShort(short int op1, short int op2) const;
	int Clamp15ToBits(int op) const;
	short int Saturate(int op) const;
	short int ScaledMult(short int op1, short int op2) const;
	
	short int Quantl(short int el, short int detl) const;
	short int Quanth(short int eh, short int deth) const;
	short int Invqal(short int il, short int detl) const;
	short int Invqah(short int ih, short int deth) const;
	short int Logscl(short int il, short int nbl) const;
	short int Logsch(short int ih, short int nbh) const;
	short int Scalel(short int nbpl) const;
	short int Scaleh(short int nbph) const;
	void Upzero(short int* dlt_ptr, short int* bl_ptr);
	void Uppol1(short int* al_ptr, short int* plt_ptr);
	void Uppol2(short int* al_ptr, short int* plt_ptr);
	short int Filtez(short int* dlt_ptr, short int* bl_ptr);
	short int Filtep(short int* rlt_ptr, short int* al_ptr);
	void QmfTx(short int xin0, short int xin1, short int& xl, short int& xh);
	short int LsbCod(short int xl);
	short int HsbCod(short int xh);

protected:

public:
	G722EncoderType();
	~G722EncoderType();
	unsigned short Encode(short int pcm_data);
	void Encode(const short int * pcm_data_ptr, short int* band_data, unsigned char* g722_data_ptr, unsigned int no_of_data);
	void Reset();

	void PrintBand() const;
};


#endif
