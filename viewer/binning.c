#include "mex.h"

/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *spec, *specout, *absciss, *babsciss, *bamp;
  double step, first, cc;

  int i,j, count, outrows, outcols, elem, rindex;
  int xrows, yrows, xcols, ycols;
  int binxrows, binyrows, binxcols, binycols;
  
  char *errmsg;
  
  /*  Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using
     mexErrMsgTxt within an if statement. It will never
     get to the else statement if mexErrMsgTxt is executed.
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */
  if(nrhs != 4) 
    mexErrMsgTxt("Four inputs required: x=binning(y, x, binx, biny).");
  if(nlhs != 1) 
    mexErrMsgTxt("One output required.");
  
  /* Get the dimensions of the ordinate array. */
  spec = mxGetPr(prhs[0]); /* Get a pointer to the ordinate */
  yrows = mxGetM(prhs[0]);
  ycols = mxGetN(prhs[0]); /* For future. */

  /* Get the dimensions of the abscissa array. */
  absciss = mxGetPr(prhs[1]); /* Get a pointer to the absciss*/
  xrows = mxGetM(prhs[1]);
  xcols = mxGetN(prhs[1]);
  step     = *(absciss+1) - *absciss;
  first    = *absciss;

  /* Get the dimensions of the abscissa array. */
  babsciss = mxGetPr(prhs[2]); /* Get a pointer to the binned absciss */
  binxrows = mxGetM(prhs[2]);
  binxcols = mxGetN(prhs[2]);

  /* Get the dimensions of the ordinate array. */
  bamp      = mxGetPr(prhs[3]); /* Get a pointer to the binned ordinate */
  binyrows  = mxGetM(prhs[3]);
  binycols  = mxGetN(prhs[3]); 

  /*We are operating along rows*/
  if(xrows != yrows) mexErrMsgTxt("length(x)!=size(Spec,2)");

  if(binxrows != binyrows || binxcols != binycols)
    mexErrMsgTxt("size(binx)~=size(biny)");


  /* Get the dimensions of the spectrum array. */
  outrows = mxGetM(prhs[0]);
  outcols = mxGetN(prhs[0]);
  
  /* Prepare output array */
  plhs[0] = mxCreateDoubleMatrix(outrows,outcols, mxREAL);
  specout = mxGetPr(plhs[0]);
  elem = outrows*outcols; 

  /* Copy input array to output */
  count = 0;
  for(i = 0; i < outrows; i++)
    for(j = 0; j < outcols; j++) {
      *(specout+count) = *(spec+count);
      count++;
      };

  /* Debug */
  /*  errmsg = "*********************************"; */
  
  /* Do binning */
  count = 0;
  for(i = 0; i < binxrows; i++)
    for(j = 0; j < binxcols; j++) 
      {
       rindex = (int)((*(babsciss+count) - first) / step);

       cc = 1.;
       if (rindex < elem && rindex >= 0 ) 
         {
          cc = (*(babsciss+count) - *(absciss+rindex)) / step;
          *(specout + rindex) = *(specout + rindex) + *(bamp+count)*(1.-cc);
         }
       if (rindex+1 < elem && rindex+1 >= 0 ) 
          *(specout + rindex + 1) = *(specout + rindex + 1) + *(bamp+count)*cc;
       count++;
      };
  
}

