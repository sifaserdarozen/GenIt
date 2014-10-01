#include "g722coder.h"
#include <iostream>
#include <cstring>    // needed for g++ for std::memcpy

G722DecoderType::G722DecoderType()
{
    Reset();
}

G722DecoderType::~G722DecoderType()
{

}

short int G722DecoderType::ConvertLongToShort(long int in_value)
{
	if (in_value > 32767)
		return 32767;
	else if (in_value < -32768)
		return -32768;
	else
		return (short int)in_value;
}

void G722DecoderType::Decode(unsigned char* g722_data_ptr, int* band_data, unsigned short* pcm_data_ptr, unsigned int no_of_data)
{
	std::memcpy((void*)&band, (void*)band_data, 45*sizeof(int));

	int dlowt;
	int rlow;
//	int ihigh;
	int wd1;
	int wd2;
	int wd3;
	unsigned char g722_data;
	
	for (unsigned int index=0; index<no_of_data; index++)
	{
	    g722_data = g722_data_ptr[index];

	    wd1 = g722_data & 0x3F;
//	    ihigh = (g722_data >> 6) & 0x03;
	    wd2 = qm6[wd1];
	    wd1 >>= 2;

	    /********************** Block 5 *******************/
	    // INVQBL (ITU page 43), compute quantized difference signal for the decoder output in the lower sub-band
	    wd2 = (band.det * wd2) >> 15;
	    // RECONS ( ITU page 41), compute reconstructed signal for the adaptive predictor
	    rlow = band.s + wd2;
	
	    /********************** Block 6 ********************/
	    // LIMIT (ITU page 44), limit the output reconstructed signal
	    if (rlow > 16383)
		    rlow = 16383;
	    else if (rlow < -16384)
		    rlow = -16384;

	    /********************** Block 2 ***********************/	
	    // INVQAL (ITU page 37), compute the quantized differences signal for the adaptive predictor in the lower sub-band
	    wd2 = qm4[wd1];
	    dlowt = (band.det * wd2) >> 15;

	    /********************** Block 3 ************************/
	    // LOGSCL (ITU page 38), update the logarithmic quantizer scale factor in the lower sub-band
	    wd2 = rl42[wd1];
	    wd1 = (band.nb * 127) >> 7;
	    wd1 += wl[wd2];
	    if (wd1 < 0)
		    wd1 = 0;
	    else if (wd1 > 18432)
		    wd1 = 18432;
	    band.nb = wd1;

	    // SCALEL (ITU page 38), compute the quantizer scale factor in the lower sub-band 
	    wd1 = (band.nb >> 6) & 31;
	    wd2 = 8 - (band.nb >> 11);
	    wd3 = (wd2 < 0)	 ?  (ilb[wd1] << -wd2)	:  (ilb[wd1] >> wd2);
	    band.det = wd3 << 2;

	    /********************** Block 4 **************************/

	    // RECONS (ITU page 41), compute reconstructed signal for the adaptive predictor
	    band.d[0] = dlowt;
	    band.r[0] = ConvertLongToShort(band.s + dlowt);

	    // PARREC (ITU page 40), compute partially reconstructed signal
	    band.p[0] = ConvertLongToShort(band.sz + dlowt);

	    // UPPOL2 (ITU page 41), update second predictor coefficient
	    int i;  // loop variable
	    for (i = 0;	 i < 3;	 i++)
		    band.sg[i] = band.p[i] >> 15;
	    wd1 = ConvertLongToShort(band.a[1] << 2);

	    wd2 = (band.sg[0] == band.sg[1])	?  -wd1	 :  wd1;
	    if (wd2 > 32767)
		    wd2 = 32767;
	    wd3 = (band.sg[0] == band.sg[2])	?  128	:  -128;
	    wd3 += (wd2 >> 7);
	    wd3 += (band.a[2]*32512) >> 15;
	    if (wd3 > 12288)
	    	wd3 = 12288;
	    else if (wd3 < -12288)
		    wd3 = -12288;
	    band.ap[2] = wd3;

	    // UPPOL1 (ITU page 42), update first predictor coefficient
	    band.sg[0] = band.p[0] >> 15;
	    band.sg[1] = band.p[1] >> 15;
	    wd1 = (band.sg[0] == band.sg[1])	?  192	:  -192;
	    wd2 = (band.a[1]*32640) >> 15;

	    band.ap[1] = ConvertLongToShort(wd1 + wd2);
	    wd3 = ConvertLongToShort(15360 - band.ap[2]);
	    if (band.ap[1] > wd3)
		    band.ap[1] = wd3;
	    else if (band.ap[1] < -wd3)
		    band.ap[1] = -wd3;

	    // UPZERO (ITU page 41), update sixth order predictor coefficients
	    wd1 = (dlowt == 0)  ?  0  :  128;
	    band.sg[0] = dlowt >> 15;
	    for (i = 1;	 i < 7;	 i++)
	    {
		    band.sg[i] = band.d[i] >> 15;
		    wd2 = (band.sg[i] == band.sg[0])  ?  wd1  :  -wd1;
		    wd3 = (band.b[i]*32640) >> 15;
		    band.bp[i] = ConvertLongToShort(wd2 + wd3);
	    }

	    // DELAYA (ITU page 38), memory block delay 
	    for (i = 6;	 i > 0;	 i--)
	    {
		    band.d[i] = band.d[i - 1];
		    band.b[i] = band.bp[i];
	    }

	    for (i = 2;	 i > 0;	 i--)
	    {
		    band.r[i] = band.r[i - 1];
		    band.p[i] = band.p[i - 1];
		    band.a[i] = band.ap[i];
	    }

	    // FILTEP (ITU page 43), compute predictor output signal, poles
	    wd1 = ConvertLongToShort(band.r[1] + band.r[1]);
	    wd1 = (band.a[1]*wd1) >> 15;
	    wd2 = ConvertLongToShort(band.r[2] + band.r[2]);
	    wd2 = (band.a[2]*wd2) >> 15;
	    band.sp = ConvertLongToShort(wd1 + wd2);

	    // FILTEZ (ITU page 42), compute predictor output signal, zeros
	    band.sz = 0;
	    for (i = 6;	 i > 0;	 i--)
	    {
		    wd1 = ConvertLongToShort(band.d[i] + band.d[i]);
		    band.sz += (band.b[i]*wd1) >> 15;
	    }
	    band.sz = ConvertLongToShort(band.sz);

	    // PREDIC (ITU page 43), compute predictor output value
	    band.s = ConvertLongToShort(band.sp + band.sz);

	    pcm_data_ptr[index]= (short int)rlow;
	}
	std::memcpy((void*)band_data, (void*)&band, 45*sizeof(int));
}

void G722DecoderType::Reset()
{
	band.s = 0;
	band.sp = 0;
	band.sz = 0;
	
	for (int k=0; k<3; k++)
	{
		band.r[k] = 0;
		band.a[k] = 0;
		band.ap[k] = 0;
		band.p[k] = 0;
	}
	
	for (int k=0; k<7; k++)
	{
		band.d[k] = 0;
		band.b[k] = 0;
		band.bp[k] = 0;
		band.sg[k] = 0;
	}
	
	band.nb = 0;
	band.det = 32;
}



G722EncoderType::G722EncoderType()
{
	Reset();
}

G722EncoderType::~G722EncoderType()
{

}

int G722EncoderType::SaturateAdd(int op1, int op2) const
{
	int out = op1 + op2;
	if ((((op1 ^ op2) & MIN_32) == 0) && ((out ^ op1) & MIN_32))
		out = (op1 < 0) ? MIN_32 : MAX_32;
	return out;
}

int G722EncoderType::SaturateSubtract(int op1, int op2) const
{
	int out = op1 - op2;
	if ((((op1 ^ op2) & MIN_32) != 0) && ((out ^ op1) & MIN_32))
		out = (op1 < 0L) ? MIN_32 : MAX_32;
	return out;
}

int G722EncoderType::ShiftRight(int op1, short int op2) const
{
		if (op2 >= 31)
			return (op1 < 0) ? -1 : 0;
		else
			return op1 >> op2;
}

int G722EncoderType::ShiftLeft(int op1, short int op2) const
{

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


short int G722EncoderType::ShiftLeftShort(short int op1, short int op2) const
{
		int result = ((int)op1) * ((int) 1 << op2);

		if ((op2 > 15 && op1 != 0) || (result != (int) ((short int) result)))
			return (op1 > 0) ? MAX_16 : MIN_16;
		else
			return (short int)result;
}


short int G722EncoderType::ShiftRightShort(short int op1, short int op2) const
{
		if (op2 >= 15)
			return (op1 < 0) ? -1 : 0;
		else
			if (op1 < 0)
				return ~((~op1) >> op2);
			else
				return op1 >> op2;
}

int G722EncoderType::Clamp15ToBits(int op) const
{
	if (op > 16383)
		return 16383;
	else if (op < -16384)
		return -16384;
	return op;
}

int G722EncoderType::MultiplyAdd(int add_op, short int mul_op1, short int mul_op2) const
{
	return SaturateAdd(add_op, ((int)mul_op1 * (int)mul_op2));
}

short int G722EncoderType::Saturate(int op) const
{
	if (op > MAX_16)
		return MAX_16;
	else if (op < MIN_16)
		return MIN_16;
	return op;
}

short int G722EncoderType::SaturateSubtractShort(short int op1, short int op2) const
{
	return Saturate (((int)op1 - op2));
}

short int G722EncoderType::SaturateAddShort(short int op1, short int op2) const
{
	return Saturate (((int)op1 + op2));
}

short int G722EncoderType::ScaledMult(short int op1, short int op2) const
{
	int product = (((int)op1 * (int)op2) & (int)(0xffff8000)) >> 15;

	if (product & (int)0x00010000)
		product |= (int)0xffff0000;

	return (short int)product;
}

short int G722EncoderType::Quantl(short int el, short int detl) const
{
	short int sil = ShiftRightShort(el, 15);
	short int wd = SaturateSubtractShort(MAX_16,(el & MAX_16));
	short int mil = 0;

	if (sil == 0)
		wd = el;

	short int val = ScaledMult(ShiftLeftShort(q6[mil], 3), detl);
	while (SaturateSubtractShort(val,wd) <= 0)
	{
		if (SaturateSubtractShort(mil, 30) == 0)
			break;
		else
		{
			mil = SaturateAddShort(mil, 1);
			val = ScaledMult(ShiftLeftShort(q6[mil], 3), detl);
		}
	}

	sil = SaturateAddShort(sil, 1);

	return misil[sil][mil];
}

short int G722EncoderType::Quanth(short int eh, short int deth) const
{
	short int sih = ShiftRightShort(eh, 15);
	short int wd = SaturateSubtractShort(MAX_16, (eh & MAX_16));

	if (sih == 0)
		wd = eh;

	short int mih = 1;

	if (SaturateSubtractShort(wd, ScaledMult(ShiftLeftShort(564, 3), deth)) >= 0)
		mih = 2;

	sih = SaturateAddShort(sih, 1);

	return misih[sih][mih];
}

short int G722EncoderType::Invqal(short int il, short int detl) const
{
	short int ril = ShiftRightShort(il, 2);
	short int wd1 = ShiftLeftShort(oq4[ril4[ril]], 3);
	short int wd2 = -wd1;

	if (risil[ril] == 0)
		wd2 = wd1;

	return ScaledMult(detl, wd2);
}

short int G722EncoderType::Invqah(short int ih, short int deth) const
{
	short int wd1 = ShiftLeftShort(oq2[ih2[ih]], 3);
	short int wd2 = -wd1;

	if (sih[ih] == 0)
		wd2 = wd1;

	return ScaledMult(wd2, deth);
}

short int G722EncoderType::Logscl(short int il, short int nbl) const
{
	short int ril = ShiftRightShort(il, 2);
	short int wd = ScaledMult(nbl, 32512);
	short int il4 = ril4[ril];
	short int nbpl = SaturateAddShort (wd, wl[il4]);

	if (nbpl < 0)
		nbpl = 0;

	if (SaturateSubtractShort(nbpl, 18432) > 0)
		nbpl = 18432;

	return nbpl;
}

short int G722EncoderType::Logsch(short int ih, short int nbh) const
{
	short int wd = ScaledMult(nbh, 32512);
	short int nbph = SaturateAddShort(wd, wh[ih2[ih]]);

	if(nbph < 0)
		nbph = 0;

	if(SaturateSubtractShort(nbph, 22528) > 0)
		nbph = 22528;

	return nbph;
}

short int G722EncoderType::Scalel(short int nbpl) const
{
	short int wd1 = ShiftRightShort(nbpl, 6) & 511;
	short int wd2 = SaturateAddShort(wd1, 64);
	return (ShiftLeftShort(SaturateAddShort(ila[wd2], 1), 2));
}

short int G722EncoderType::Scaleh(short int nbph) const
{
	short int wd = ShiftRightShort(nbph, 6) & 511;
	return ShiftLeftShort(SaturateAddShort(ila[wd], 1), 2);
}


void G722EncoderType::Upzero(short int* dlt_ptr, short int* bl_ptr)
{
	short int wd1 = 128;

	if (dlt_ptr[0] == 0)
		 wd1 = 0;

	short int sg0 = ShiftRightShort(dlt_ptr[0], 15);

	for (short int i = 6; i > 0; i--)
	{
		short int wd2 = SaturateSubtractShort (0, wd1);
		if(sg0 == ShiftRightShort(dlt_ptr[i], 15))
			wd2 = SaturateAddShort (0, wd1);

		bl_ptr[i] = SaturateAddShort(wd2, ScaledMult(bl_ptr[i], 32640));
		dlt_ptr[i] = dlt_ptr[i - 1];
	}
}

void G722EncoderType::Uppol1(short int* al_ptr, short int* plt_ptr)
{
	short int sg0 = ShiftRightShort(plt_ptr[0], 15);
	short int sg1 = ShiftRightShort(plt_ptr[1], 15);
	short int wd1 = -192;

	if (SaturateSubtractShort(sg0, sg1) == 0)
		wd1 = 192;

	short int wd2 = ScaledMult (al_ptr[1], 32640);
	short int apl1 = SaturateAddShort(wd1, wd2);
	short int wd3 = SaturateSubtractShort(15360, al_ptr[2]);

	if (SaturateSubtractShort(apl1, wd3) > 0)
		apl1 = wd3;
	else if (SaturateAddShort(apl1, wd3) < 0)
		apl1 = -wd3;

	/* Shift of the plt signals */
	plt_ptr[2] = plt_ptr[1];
	plt_ptr[1] = plt_ptr[0];
	al_ptr[1] = apl1;
}

void G722EncoderType::Uppol2(short int* al_ptr, short int* plt_ptr)
{
	short int sg0 = ShiftRightShort(plt_ptr[0], 15);
	short int sg1 = ShiftRightShort(plt_ptr[1], 15);
	short int sg2 = ShiftRightShort(plt_ptr[2], 15);
	short int wd1 = ShiftLeftShort(al_ptr[1], 2);
	short int wd2 = SaturateAddShort(0, wd1);

	if (SaturateSubtractShort(sg0, sg1) == 0)
		wd2 = SaturateSubtractShort(0, wd1);

	wd2 = ShiftRightShort(wd2, 7);
	short int wd3 = -128;

	if (SaturateSubtractShort(sg0, sg2) == 0)
		wd3 = 128;

	short int wd4 = SaturateAddShort (wd2, wd3);
	short int wd5 = ScaledMult(al_ptr[2], 32512);
	short int apl2 = SaturateAddShort(wd4, wd5);

	if (SaturateSubtractShort(apl2, 12288) > 0)
		apl2 = 12288;

	if (SaturateSubtractShort(apl2, -12288) < 0)
		apl2 = -12288;

	al_ptr[2] = apl2;
}

short int G722EncoderType::Filtez(short int* dlt_ptr, short int* bl_ptr)
{
	short int szl = 0;

	for (short int i = 6; i > 0; i--)
	{
		short int wd = SaturateAddShort(dlt_ptr[i], dlt_ptr[i]);
		wd = ScaledMult(wd, bl_ptr[i]);
		szl = SaturateAddShort(szl, wd);
	}
	return szl;
}

short int G722EncoderType::Filtep(short int* rlt_ptr, short int* al_ptr)
{
	// shift of rlt
	rlt_ptr[2] = rlt_ptr[1];		
	rlt_ptr[1] = rlt_ptr[0];		

	short int wd1 = SaturateAddShort(rlt_ptr[1], rlt_ptr[1]);
	wd1 = ScaledMult(al_ptr[1], wd1);
	short int wd2 = SaturateAddShort(rlt_ptr[2], rlt_ptr[2]);
	wd2 = ScaledMult(al_ptr[2], wd2);
	return SaturateAddShort(wd1, wd2);
}

void G722EncoderType::QmfTx(short int xin0, short int xin1, short int& xl, short int& xh)
{
	int accuma;
	int accumb;
	int comp_low;
	int comp_high;

	const short int* pcoef = coef_qmf;
	short int* pdelayx = band.qmf_tx_delayx;

	/* Saving past samples in delay line */
	band.qmf_tx_delayx[1] = xin1;
	band.qmf_tx_delayx[0] = xin0;

	accuma = (int)*pcoef++, (int)*pdelayx++;
	accumb = (int)*pcoef++, (int)*pdelayx++;

	for(short int i = 1; i < 12; i++)
	{
		accuma = MultiplyAdd(accuma, *pcoef++, *pdelayx++);
		accumb = MultiplyAdd(accumb, *pcoef++, *pdelayx++);
	}

	/* Descaling and shift of the delay line */
	for (short int i = 0; i < 22; i++)
		band.qmf_tx_delayx[23 - i] = band.qmf_tx_delayx[21 - i];

	comp_low = SaturateAdd (accuma, accumb);
	comp_low = SaturateAdd (comp_low, comp_low);
	comp_high = SaturateSubtract (accuma, accumb);
	comp_high = SaturateAdd (comp_high, comp_high);
	xl = Clamp15ToBits (ShiftRight(comp_low, 16));
	xh = Clamp15ToBits (ShiftRight(comp_high, 16));
}

short int G722EncoderType::LsbCod(short int xl)
{
	short int il = Quantl (SaturateSubtractShort (xl, band.sl), band.detl);
	band.dlt[0] = Invqal (il, band.detl);
	short int nbpl = Logscl (il, band.nbl);
	band.nbl = nbpl;
	band.detl = Scalel (nbpl);
	band.plt[0] = SaturateAddShort (band.dlt[0], band.szl);   /* parrec */
	band.rlt[0] = SaturateAddShort (band.sl, band.dlt[0]);    /* recons */
	Upzero (band.dlt, band.bl);
	Uppol2 (band.al, band.plt);
	Uppol1 (band.al, band.plt);
	band.szl = Filtez(band.dlt, band.bl);
	band.spl = Filtep(band.rlt, band.al);
	band.sl = SaturateAddShort (band.spl, band.szl);          /* predic */

	/* Return encoded sample */
	return il;
}

short int G722EncoderType::HsbCod(short int xh)
{
	short int ih = Quanth (SaturateSubtractShort(xh, band.sh), band.deth);
	band.dh[0] = Invqah (ih, band.deth);
	short int nbph = Logsch (ih, band.nbh);
	band.nbh = nbph;
	band.deth = Scaleh (nbph);
	band.ph[0] = SaturateAddShort(band.dh[0], band.szh);   /* parrec */
	band.rh[0] = SaturateAddShort(band.sh, band.dh[0]);    /* recons */
	Upzero (band.dh, band.bh);
	Uppol2 (band.ah, band.ph);
	Uppol1 (band.ah, band.ph);
	band.szh = Filtez (band.dh, band.bh);
	band.sph = Filtep (band.rh, band.ah);
	band.sh = SaturateAddShort(band.sph, band.szh);        /* predic */

	return ih;
}

void G722EncoderType::PrintBand() const
{
/*
	short int al[3];
	short int plt[3]; 
	short int rlt[3];
	short int ah[3];
	short int ph[3];
	short int rh[3];
	short int bl[7];
	short int dlt[7]; 
	short int bh[7];
	short int dh[7];
	short int qmf_tx_delayx[24];
	short int qmf_rx_delayx[24];
*/

	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << band.detl << " " << band.deth << " " << band.nbl << " " << band.sl << " " << band. spl << " " << band.szl << " " << band.nbh << " " << band.sh << " " << band.sph << " " << band.szh << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
}

void G722EncoderType::Encode(const short int* pcm_data_ptr, short int* band_data, unsigned char* g722_data_ptr, unsigned int no_of_data)
{
	std::memcpy((void*)&band, (void*)band_data, (sizeof(short int) * 104));

	//PrintBand();

	unsigned char g722_data;
	short int xin1;
	short int xin0;
	short int xl;
	short int il;
	short int xh;
	short int ih;

	for (unsigned int index = 0; index < no_of_data; index += 2)
	{
		xin1 = pcm_data_ptr[index];
		xin0 = pcm_data_ptr[index + 1];

		// Calculation of the synthesis QMF samples 
		// qmf_tx (xin0, xin1, &xl, &xh, encoder);
		QmfTx(xin0, xin1, xl, xh);

		// Call the upper and lower band ADPCM encoders
		// il = lsbcod (xl, 0, encoder);
		il = LsbCod(xl);
		// ih = hsbcod (xh, 0, encoder);
		ih = HsbCod(xh);

		// Mount the output G722 codeword: bits 0 to 5 are the lower-band
		// portion of the encoding, and bits 6 and 7 are the upper-band
		// portion of the encoding
		// code[i] = s_and(add(shl(ih, 6), il), 0xFF);
		g722_data = (unsigned char) SaturateAddShort(ShiftLeftShort(ih, 6), il);

		g722_data_ptr[index/2] = g722_data;


//		g722_data_ptr[index/2] = (unsigned char) ((xl + xh)/2);
	}
	std::memcpy((void*)band_data, (void*)&band, (sizeof(short int) * 104));
	//PrintBand();
}

void G722EncoderType::Reset()
{
	// initialization due to lower band
	band.detl = 32;
	band.deth = 8;
	band.sl = 0;
	band.sh = 0;
	band.spl = 0;
	band.sph = 0;
	band.szl = 0;
	band.szh = 0;
	band.nbl = 0;
	band.nbh = 0;

	for (int i = 0; i < 2; i++)
	{
		band.al[i] = 0;
		band.ah[i] = 0;
		band.plt[i] = 0;
		band.rlt[i] = 0;
		band.ph[i] = 0;
		band.rh[i] = 0;
	}

	for (int i = 0; i < 7; i++)
	{
		band.bl[i] = 0;
		band.dlt[i] = 0;
		band.bh[i] = 0;
		band.dh[i] = 0;
	}

	for (int i = 0; i < 24; i++)
	{
		band.qmf_tx_delayx[i] = 0;
		band.qmf_rx_delayx[i] = 0;
	}
}
