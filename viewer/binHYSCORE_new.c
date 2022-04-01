/*
binHYSCORE(Spec, wa, wb, Amp, MaxFreq)
           binning to a 2D matrix Spec
	designed to speed up kv_hyscorefd from MatLab
	
	2020 Update (yeah, COVID year). Apparently, "they" changed the way complex arrays are stored in memory
	Old ways don't work anymore... bustards...
	If you compile this with just "mex binHYSCORE_new.c", it is not going to work...why?... because reasons...
	To make this script work, you need to compile with "mex -R2018a binHYSCORE_new.c" 
			-- makes no sense, I know ... but who am I to question gods.
	For future me ... if they change the default mex setting and nothing works, google "mxGetComplexDoubles" 
	and see what setting "du jour" they offer 
	>> alsi 2020 PennState University
*/
#include "mex.h"
#include <math.h>

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nNums;
  double *omega_a, *omega_b, dx, nPointsD, MaxFreq;
  int idx, idx1, idx2;
  int a, r, c, i, j, count, nPoints, iPeak;
  bool realOnly;
  /*-----------------------------------------------------------------
     Check number of inputs
    ----------------------------------------------------------------- */
  if (nrhs!=5)
     mexErrMsgTxt("Wrong number of input arguments. Usage: binHYSCORE(Spec, wa, wb, Amp, MaxFreq)");

	/*-----------------------------------------------------------------
   input : 
	0. Matrix to bin to
	1. omega_a
	2. omega_b
	3. Amp
	4. MaxFreq
    -----------------------------------------------------------------
    */
	
 #if MX_HAS_INTERLEAVED_COMPLEX   
	mxComplexDouble *Ampl, *Spec;
	//mxDouble *rAmpl
	a = 0;
    /* Spec: spectral storage array */
    if (mxGetNumberOfDimensions(prhs[a])!=2)
		mexErrMsgTxt("Spec must be a 2D array!");
	Spec = mxGetComplexDoubles(prhs[a]);
	nPoints = mxGetM(prhs[a]);
	nPointsD  = (double) nPoints;
	if (mxGetN(prhs[a])!=nPoints) mexErrMsgTxt("Spec must be a 2D square array!");
	if (!mxIsComplex(prhs[a])) mexErrMsgTxt("Spec must be complex!");	 
 
    a = 1;
    /* omega_a*/
    if (mxIsComplex(prhs[a])) mexErrMsgTxt("omega_a must be real!");
    omega_a = mxGetDoubles(prhs[a]);
	
	nNums = mxGetNumberOfElements(prhs[a]);

    a = 2;
    /* omega_b */
    if (mxIsComplex(prhs[a])) mexErrMsgTxt("omega_b must be real!");
    if (mxGetNumberOfElements(prhs[a])!=nNums) mexErrMsgTxt("omega_a and omega_b must have the same number of elements!");
	omega_b = mxGetDoubles(prhs[a]);
	
    a = 3;
    /* Amplitudes: peak amplitudes */
    if (mxGetNumberOfElements(prhs[a])!=nNums) mexErrMsgTxt("Amplitudes has wrong size!");
	
	if (!mxIsComplex(prhs[a])) 
	{	//mexErrMsgTxt("Amplitudes must be complex!");
		Ampl = mxGetDoubles(prhs[a]);
		realOnly = true;
	}
	else
	{	Ampl = mxGetComplexDoubles(prhs[a]);
		realOnly = false;
	}
	
    a = 4;
    /* Max Frequency*/
    if (mxIsComplex(prhs[a])) mexErrMsgTxt("MaxFreq must be real!");
	MaxFreq = mxGetScalar(prhs[a]);
	
   iPeak = 0;
   dx = (nPointsD-1.0)/2.0/ MaxFreq;

    
  for (r=0; r<nNums; r++) {
    idx1 = floor((omega_a[r]+MaxFreq )*dx);
	
	if (idx1<0) continue;
	if (idx1>nPoints-1) continue;
	idx2 = floor((omega_b[r]+MaxFreq )*dx);
	if (idx2<0) continue;
	if (idx2>nPoints-1) continue;
		
    idx = idx1 + idx2*nPoints;
	
	//mexPrintf("--%d %f  ", idx, Ampl[r]);
	if (realOnly)
	{	Spec[idx].real += Ampl[r].real;
		//  for symmetric HYSCORE 
		idx = idx2 + idx1*nPoints;
		Spec[idx].real += Ampl[r].real;
	}
	else
    {
		Spec[idx].real += Ampl[r].real;
		Spec[idx].imag += Ampl[r].imag;
			
		//mexPrintf("-- %f  \n", Spec[idx]);
		
		//  for symmetric HYSCORE 
		idx = idx2 + idx1*nPoints;
		Spec[idx].real += Ampl[r].real;
		Spec[idx].imag += Ampl[r].imag;
	}
	
	
  }
#else
	mexPrintf(" This routine was not compiled properly. Need ""mex -R2018a binHYSCORE_new.c"" because ... reasons...\n");
#endif

}
