/*
binHYSCORE(Spec, wa, wb, Amp, MaxFreq)
           binning to a 2D matrix Spec
	desgned to speed up kv_hyscorefd from MAtLab
	
	to complile in MatLab, chose "lcc" compliler from "mex -setup" 
	and type "mex binHSYCORE.c " from the directory of the Kazan Viewer
*/
#include "mex.h"
#include <math.h>

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nNums;
  double *rAmpl, *iAmpl, *rSpec, *iSpec, MaxFreq;
  double *omega_a, *omega_b, dx, nPointsD;
  int ridx, idx, idx1, idx2;

  int a, r, c, nPoints, iPeak;
  mxArray *imagZeros, *outR, *outL ;

  /*-----------------------------------------------------------------
     Check number of input and output arguments
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

  a = 0;
  /* Spec: spectral storage array */
  if (mxGetNumberOfDimensions(prhs[a])!=2)
    mexErrMsgTxt("Spec must be a 2D array!");

  rSpec = mxGetPr(prhs[a]);
    nPoints = mxGetM(prhs[a]);
  nPointsD  = (double) nPoints;
  if (mxGetN(prhs[a])!=nPoints)
    mexErrMsgTxt("Spec must be a 2D square array!");
  if (!mxIsComplex(prhs[a]))
    mexErrMsgTxt("Spec must be complex!");	
  iSpec = mxGetPi(prhs[a]);

	
  a = 1;
  /* omega_a*/
  omega_a = mxGetPr(prhs[a]);
  if (mxIsComplex(prhs[a]))
    mexErrMsgTxt("omega_a must be real!");
  nNums = mxGetNumberOfElements(prhs[a]);

  a = 2;
  /* omega_b */
  omega_b = mxGetPr(prhs[a]);
  if (mxIsComplex(prhs[a]))
    mexErrMsgTxt("omega_b must be real!");
  if (mxGetNumberOfElements(prhs[a])!=nNums)
    mexErrMsgTxt("omega_a and omega_b must have the same number of elements!");

  a = 3;
  /* Amplitudes: peak amplitudes */
  if (mxGetNumberOfElements(prhs[a])!=nNums)
    mexErrMsgTxt("Amplitudes has wrong size!");
  rAmpl = mxGetPr(prhs[a]);
  imagZeros = mxCreateDoubleMatrix(nNums,1,mxREAL);
  iAmpl = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);

    a = 4;
  /* Max Frequency*/
  MaxFreq = mxGetScalar(prhs[a]);
  if (mxIsComplex(prhs[a]))
    mexErrMsgTxt("MaxFreq must be real!");

	
  iPeak = 0;
  dx = (nPointsD-1.0)/2.0/ MaxFreq;

  for (r=0; r<nNums; r++) {
     //ind1 = 1+floor( (Freqs(:, 1)+Exp.MaxFreq)/(abs(2*Exp.MaxFreq/(Exp.nPoints-1))) )
    idx1 = floor((omega_a[r]+MaxFreq )*dx);
	
	if (idx1<0) continue;
	if (idx1>nPoints-1) continue;
	idx2 = floor((omega_b[r]+MaxFreq )*dx);
	if (idx2<0) continue;
	if (idx2>nPoints-1) continue;
		
    idx = idx1 + (idx2)*nPoints;
    rSpec[idx] += rAmpl[r];
    iSpec[idx] += iAmpl[r];
	
	/*  for symmetric HYSCORE */
	idx = idx2 + (idx1)*nPoints;
    rSpec[idx] += rAmpl[r];
    iSpec[idx] += iAmpl[r];
	
	//iPeak++;
  }
 // rSpec[0] = iPeak;
  mxDestroyArray(imagZeros);

  
  //mxSetPr(plhs[0], rSpec);
  //mxSetPi(plhs[0], iSpec);
}
