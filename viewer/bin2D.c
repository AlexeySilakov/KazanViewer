/*
bin2D(Spec, ind1, ind2, Amp)
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
  double *rAmpl, *iAmpl, *rSpec, *iSpec, MaxFreq, *rOut, *iOut;
  double *omega_a, *omega_b, dx, nPointsD;
  int ridx, idx;
  double *idx1, *idx2;
  bool iflag;
  int a, r, c, nPoints, iPeak;
  mxArray *imagZeros, *outR, *outL ;

  /*-----------------------------------------------------------------
     Check number of input and output arguments
    ----------------------------------------------------------------- */
  if (nrhs!=4)
     mexErrMsgTxt("Wrong number of input arguments. Usage: bin2D(Spec, ind1, ind2, Amp)");

  /*-----------------------------------------------------------------
   input : 
	0. Matrix to bin to
	1. ind1
	2. ind2
	3. Amp
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
	iflag = false;
  else //   mexErrMsgTxt("Spec must be complex!");	
  { 
	iflag = true;
	iSpec = mxGetPi(prhs[a]);
  }

	
  a = 1;
  /* idx1*/
  idx1 = mxGetPr(prhs[a]);
   /* if (mxIsComplex(prhs[a]))
    mexErrMsgTxt("idx1 must be real!");*/
  nNums = mxGetNumberOfElements(prhs[a]);

  a = 2;
  /* idx2 */
  idx2 = mxGetPr(prhs[a]);
  /*if (mxIsComplex(prhs[a]))
    mexErrMsgTxt("idx2 must be real!");*/
  if (mxGetNumberOfElements(prhs[a])!=nNums)
    mexErrMsgTxt("idx1 and idx2 must have the same number of elements!");

  a = 3;
  /* Amplitudes: peak amplitudes */
  if (mxGetNumberOfElements(prhs[a])!=nNums)
    mexErrMsgTxt("Amplitudes has wrong size!");
  rAmpl = mxGetPr(prhs[a]);
  imagZeros = mxCreateDoubleMatrix(nNums,1,mxREAL);
  if (iflag)
	iAmpl = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  
  if (iflag)
  {
	for (r=0; r<nNums; r++) {

		
    idx = (idx1[r]-1 + (idx2[r]-1)*nPoints); // 
	// mexPrintf("\n idx1 = %d",  idx);
	if (idx<0) continue; 
	if (idx>nPoints*nPoints -1) continue;
    rSpec[idx] += rAmpl[r];
    iSpec[idx] += iAmpl[r];

	}
	}
  else
	for (r=0; r<nNums; r++) {
	
    idx = (idx1[r]-1 + (idx2[r]-1)*nPoints); // (int) (*(idx1 + r) + (*(idx2 + r))*nPoints);
	//mexPrintf("\n idx1 = %d", idx);
	if (idx<0) continue; 
	if (idx>nPoints*nPoints -1) continue;
    rSpec[idx] += rAmpl[r];

	}  
	
 if (nlhs==1)
	plhs[0] = mxDuplicateArray(prhs[0]);

  mxDestroyArray(imagZeros);

}
