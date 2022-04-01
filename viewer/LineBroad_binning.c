#include "mex.h"

/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *ydata, *specout, *x, *dim, *bounds;
  double step1, first1, step2, first2, cc;

  int i,j, count, outrows, outcols, nbounds, nfwhh;
  int nxrows, nxcols, ny, dimrws, make2d;
  int binxrows, binyrows, binxcols, binycols;
  
  char *errmsg;
  
  /*  Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using
     mexErrMsgTxt within an if statement. It will never
     get to the else statement if mexErrMsgTxt is executed.
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */
  if(nrhs != 5) 
    mexErrMsgTxt("Four inputs required: yout=binning(y, x, dim, bounds, fwhh).");
  if(nlhs != 1) 
    mexErrMsgTxt("One output required.");
  
  /* Get the dimensions of the ordinate array. */
  ydata = mxGetPr(prhs[0]); /* Get a pointer to the ordinate */
  ny = mxGetM(prhs[0]);
  //ycols = mxGetN(prhs[0]); /* For future. */

  /* Get the dimensions of the abscissa array. */
  x  = mxGetPr(prhs[1]); /* Get a pointer to the absciss*/
  nxrows = mxGetM(prhs[1]);
  nxcols = mxGetN(prhs[1]);
    
  if (nxrows !=ny),
	mexErrMsgTxt("Lenght of the input y's is not equal to input x's");
	
  dim = mxGetPr(prhs[2]); /*dimentions of the requered data*/
  dimrws = mxGetM(prhs[2]);
  
  if (dimrws ==2),
	{ if ( (*(dim+1) > 1) && (nxcols<2) ),
		mexErrMsgTxt("Length of the second dimension is not equal to "1" but X does not have second column");
    }
  bounds = mxGetPr(prhs[3]); 
  nbounds = mxGetM(prhs[3]);		
  
  if (nbounds !=2) || (nbounds !=4) 
	mexErrMsgTxt("Bounds should include at either to numbers [xmin, xmax] or four [x1min, x1max, x2min, x2max]");
  
  make2d = 0;
    
  step1     = ( bounds[1]-bounds[0] )/dim[0] ;
  first1    = bounds[0];
  
  step2 = step1;
  first2 = first2;
  
  if (nbounds ==2) && (dimrws==2)
	{ if (dim[1]>1)
	 {	make2d = 1;
		  step2     = ( bounds[1]-bounds[0] )/dim[1] ;
          first2    = bounds[0];
	 }
	}
 else if (nbounds ==4)
    {	step2     = ( bounds[3]-bounds[2] )/dim[1] ;
	    first2    = bounds[2];	
    }
	
  /* Get the dimensions of the spectrum array. */
  outrows = dim[0];
  if (dimrws ==2)
	outcols = dim[1];
  else
	outcols = 1;
  
  fwhh = mxGetPr(prhs[4]); /*dimentions of the requered data*/
  nfwhh = mxGetM(prhs[4]);
  
  fwhh1 = fwhh[0];
  fwhh2 = fwhh[0];
  if (nfwhh==2)
	fwhh2 = fwhh[1];
	
  /* Prepare output array */
  plhs[0] = mxCreateDoubleMatrix(outrows,outcols, mxREAL);
  specout = mxGetPr(plhs[0]);
  
  elem = outrows*outcols; 

  /* Debug */
  /*  errmsg = "*********************************"; */
  
  /* Do stuff */
  
  if (outcols ==1)
  {
  for (i= 0; i < nxrows; i++)
	{
		1Dgauss(specout, ydata[i], (x[i]-first1)/step1, fwhh1/step1, outcols)
	}
  }
  else
  {
  count = 0;
  for (i=0; i<nxrows; i++)
	for (j=0; j < nxcols; j++)
	{
		2Dgauss(specout, ydata[i], (x[i]-first1)/step1, (x[i+nxrows]-first2)/step2, fwhh1/step1, fwhh2/step2, outcols, outrows)
	}
  }
  
}

void 1Dgauss(double spec[], double amp, double x, double fwhh, int npoints)
/* makes gaussian line with x [points] as a center and fwhh [points] */
{
for (i=0; i<npoints; i++)
{
	spec[i] = spec[i] + exp(-2*(i-x)**2*2*log(2)/fwhh**2);
};
}

void 2Dgauss(double spec[], double y, double x1, double x2, double fwhh1, double fwhh2, int npoints1, int npoints2)
{
int count, i, j;
count = 0;

for (j=0; j<npoints2; j++)
 for (i=0; i<npoints1; i++)
 {
 	spec[count] = spec[count] + y*exp(-2*(i-x1)**2*2*log(2)/fwhh1**2) *exp(-2*(j-x2)**2*2*log(2)/fwhh2**2);
	count++;
 };
}



