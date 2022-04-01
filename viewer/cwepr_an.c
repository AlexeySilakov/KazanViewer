/* CW EPR first order simulation program
Alexey Silakov MPI-BAC 2011*/

#include "mex.h"
#include "math.h"

/*double ellip(double k)
{
    double a0, b0, c0, an, bn, cn, KK;
    int cnt;
    double pi = 3.141592653589793;
    
    a0 = 1.0;
    b0 = sqrt(1.0-k);
    c0 = sqrt(k);
    an = a0;
    bn = b0;
    cn = c0;
    cnt = 0;
    KK=  0.0;
    while (cn>1e-15)
    {
        an = (a0+b0)/2.0;
        bn = sqrt(a0*b0);
        cn = (a0-b0)/2.0;
        a0 = an;
        b0 = bn;
        c0 = cn;
        cnt++;
        if (cnt>1000) break; // just a protection
    }
    KK = pi/2.0/an;
    //mexPrintf("1: %f\n", KK );
    return KK;
}
*/
double amps(double BB, double tB1, double tB2, double tB3)
{
    double a0, b0, c0, an, bn, cn, KK, S1, k;
    int cnt;
    double pi = 3.141592653589793;
    
    if (BB<=tB2)
    { S1 = 2.0/pi*tB1*tB2*tB3 /BB/BB / sqrt( (tB1*tB1 - tB2*tB2)*(BB*BB - tB3*tB3) );
    k = (tB1*tB1-BB*BB)*(tB2*tB2-tB3*tB3)/(tB1*tB1-tB2*tB2)/(BB*BB - tB3*tB3);
    }
    else if (BB>tB2)
    {
        S1 = 2/pi*tB1*tB2*tB3 / BB/BB / sqrt( (tB2*tB2-tB3*tB3)*(tB1*tB1 - BB*BB) );
        k = (tB1*tB1-tB2*tB2)*(BB*BB - tB3*tB3)/(tB1*tB1-BB*BB)/(tB2*tB2-tB3*tB3);
    }
    
        a0 = 1.0;
    b0 = sqrt(1.0-k);
    c0 = sqrt(k);
    an = a0;
    bn = b0;
    cn = c0;
    cnt = 0;
    KK=  0.0;
    while (cn>1e-15)
    {
        an = (a0+b0)/2.0;
        bn = sqrt(a0*b0);
        cn = (a0-b0)/2.0;
        a0 = an;
        b0 = bn;
        c0 = cn;
        cnt++;
        if (cnt>1000) break; // just a protection
    }
    KK = pi/2.0/an;
    
///    KK = ellip(k1_2);
    return S1*KK;
    
}
/* Routine for calculating EPR spectra based on analytical solution */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	double B1, B2, B3, tB1, tB2, tB3, pi, nPointsD, nHStrainD, acc, didx1;
	double dx, SAmp, BB, S1, S2, k1_2, k2_2, AA, AA1, dii, ang, sth2, dnphi, djj;
	double *rSpec, *X, *Bizz, *dB, *HStrain, *Xax, *npp, dzz;
	int a, ii, jj, nphi, bc, nHStrain, idx1, idx2, idx3, nPoints, nNums, cnt;
    int m, m2, cend, ch, der, nel, zz;
    double *ss, lw, gamma, gamma1, dch, *func,lorgau, nsmooth;
    //mwSize buflen;
    
    mxArray *ffff, *ffff1;
    
	pi = 3.141592653589793;
	  a = 0;	
	  /* X axis */
	  Xax = mxGetPr(prhs[a]);
	  nPoints = mxGetNumberOfElements(prhs[a]);

	  nPointsD  = (double) nPoints;
	  
	 /* a = 1;
	  // Spec: spectral storage array 
	  rSpec = mxGetPr(prhs[a]);
	  nPoints = mxGetNumberOfElements(prhs[a]);
	 */ 
	 plhs[0] = mxCreateDoubleMatrix(nPoints, 1, mxREAL); 	
     rSpec = mxGetPr(plhs[0]);
     
     for (ii=0; ii<nPoints; ii++)
	 rSpec[ii] = 0.0;
     
	  a++;
	  // [B1, B2, B3]
	  Bizz = mxGetPr(prhs[a]);
	  nNums = mxGetNumberOfElements(prhs[a]);
		B1 = *Bizz;
		B2 = *(Bizz+1);
		B3 = *(Bizz+2);
        
	  a++;
	  /* dB */
	  dB = mxGetPr(prhs[a]);

	  a++;
	  /* HStrain */
	  HStrain = mxGetPr(prhs[a]);
		
      // defaut value
      lw = 0.0;
	  if (nrhs>(a-1))
      {    a++;
          lw =  *mxGetPr(prhs[a]);
      }
      
      // defaut values
      nHStrain = 30; // nGrid for HStrain
      nHStrainD = 30.0;
      der = 0; // 0 - absorption, 1 - first derivative, 2 - second derivative
      lorgau = 0.0; // for pseudo voight lor/gau ratio
      if (nrhs>(a-1))
      {
        a++;
        nel = mxGetNumberOfElements(prhs[a]);
        if (nel!=3) return;
	  
        npp =  mxGetPr(prhs[a]);
      
        /* HStrain numbers*/
        nHStrainD  = (double) npp[0];
        nHStrain  = (int) npp[0];
        // prepare spectrum
          	/*mexPrintf("%d\n", nrhs);
            return;	  */

        der = (int) npp[1]; 
        lorgau = npp[2];
       // mexPrintf(" %f\n", nHStrainD  );
       // mexPrintf(" %d\n", der  );
       // mexPrintf(" %f\n", lorgau  );
       // return;
      }
     

		
	//nphi = 201;
	
	for (ii=0; ii<nHStrain; ii++)
	{
        dii = (double) ii;
		dx = (dii - (nHStrainD-0.5)/2.0)/ nHStrainD;
		tB1 = B1 + 4.0*HStrain[0]*dx;
		tB2 = B2 + 4.0*HStrain[1]*dx;
		tB3 = B3 + 4.0*HStrain[2]*dx;
        if (tB1>tB2)
        {        
            tB1 = B2 + 4.0*HStrain[1]*dx;
            tB2 = B1 + 4.0*HStrain[0]*dx;
        }       
        if (tB3<tB2)
        {
            tB3 = tB2;
            tB2 = B3 + 4.0*HStrain[2]*dx;
        }
        if (tB1==tB2) tB1 = tB1-0.0000001;
        if (tB2==tB3) tB3 = tB3+0.0000001;
       
		idx1 = ceil((tB1 - Xax[0])/dB[0]);
		idx2 = floor((tB2 - Xax[0])/dB[0]);
		idx3 = floor((tB3- Xax[0])/dB[0]);
        gamma = 0.25;
		SAmp = sqrt(2.0/pi)/gamma*exp(-2.0* dx*dx/ gamma/gamma);

        
        if (idx1<=0) idx1=1;
        if (idx3>=nPoints) idx3 = nPoints-1;
        
        
		for (bc=idx1; bc<=idx3; bc++)
		{
            BB = (double) bc;
			BB = BB* dB[0]+Xax[0];
            
            
			AA = amps(BB, tB1, tB2, tB3);
           // mexPrintf(" 2: %f\n", KK );
            //mexPrintf(" 3: %f\n", ellip(k1_2) );
            //return;
           // return;
            ///
			
            // for smooth edges
            if ((bc==idx1)&(bc!=0))
            {
                didx1 = ((double) idx1 - (tB1 - Xax[0])/dB[0] );
                rSpec[bc-1] = rSpec[bc-1] +AA*SAmp*didx1;
            }
            if ((bc==idx3)&(bc<nPoints))
            {
                didx1 = ((tB3 - Xax[0])/dB[0] -(double) idx3);
                rSpec[bc+1] = rSpec[bc+1] +AA*SAmp*didx1;
            }
            if ((bc==idx2)|(bc==idx2+1)) // to smooth the "inf" problem
            {
                nsmooth = 101.0;
                AA = 0.0;
                for (zz = 0; zz<nsmooth-1; zz++)
                {
                    BB = (double) bc;
                    dzz = (double) zz;

                    BB = BB*dB[0]+Xax[0] + dB[0]*(dzz/(nsmooth-1) - 0.5); // half a point down to up
                    AA = AA  + amps(BB, tB1, tB2, tB3);
                }
                AA = AA/(nsmooth-1);
            }
                
            rSpec[bc] = rSpec[bc] +AA*SAmp;
			//return;
		}	
	}
          
    if (lw >0.0)
    {

         
          m = 2*floor(4.0*lw/ *dB)+1;
          m2= (m-1)/2;
          
          ffff = mxCreateDoubleMatrix((nPoints+2*m2), 1, mxREAL);
          ffff1 = mxCreateDoubleMatrix(m, 1, mxREAL);
          ss = mxGetPr(ffff);
          func =  mxGetPr(ffff1);
          //nn = nPoints+2*m2;
          
          //for (ii=0, ii<m, ii++) hh[ii]=h[m-ii+1];
          
          //y = zeros(1,nn);
          cend = 0;
          
          for (ii=0; ii<m2; ii++)
          {
              ss[cend] = rSpec[0];
              
              cend++;
          }
            
          for(ii=0; ii<nPoints; ii++) 
          {
              ss[cend] = rSpec[ii];
              cend++;
          }
          
          for (ii=0; ii<m2; ii++)
          {
              ss[cend] = rSpec[nPoints-1];
              cend++;
          }
          
          gamma = lw/ *dB / 1.177410022515475; // 1/sqrt(2*ln(2))
          gamma1 = lw/ *dB;
          if (der==0)
              for (ii=0; ii<m; ii++)
              {    dch= (double) ii-m2;
                  func[ii] = (1-lorgau)*exp( -2*dch*dch / gamma/gamma)*sqrt(2.0/pi)/gamma;
                  func[ii] = func[ii] + lorgau*gamma1/2/pi/(dch*dch + gamma1*gamma1);
                }
          else if (der==1)
              for (ii=0; ii<m; ii++)
              {    dch= (double) ii-m2;
                   func[ii]  = 4.0*dch*exp( -2*dch*dch / gamma/gamma)*sqrt(2.0/pi)/gamma/gamma/gamma;
                   func[ii] = func[ii] + lorgau*gamma1*dch/pi/(dch*dch + gamma1*gamma1)/(dch*dch + gamma1*gamma1);
                  
              }
          else
              for (ii=0; ii<m; ii++)
              {    dch= (double) ii-m2;
                 func[ii]  = -4.0*exp( -2*dch*dch / gamma/gamma)*sqrt(2.0/pi)/gamma/gamma/gamma;
                 func[ii]  = func[ii] +16.0*dch*dch*exp( -2*dch*dch / gamma/gamma)*sqrt(2.0/pi)/gamma/gamma/gamma/gamma/gamma;

                 func[ii] = func[ii] + lorgau*4*gamma1*dch*dch/pi/(dch*dch + gamma1*gamma1)/(dch*dch + gamma1*gamma1)/(dch*dch + gamma1*gamma1);
                 func[ii] = func[ii] - lorgau*gamma1/2/pi/(dch*dch + gamma1*gamma1)/(dch*dch + gamma1*gamma1);
              }
          
          /*mexPrintf("aaa %f\n", *(func+ii)  );
          return; */        
          for (ii = 0; ii<(nPoints); ii++)
          {     rSpec[ii]=0.0;
                for (ch=0; ch<m; ch++) 
                {
                    dch= (double) ch-m2;
                    rSpec[ii] = rSpec[ii]+ ss[ii+ch]*func[ch];
                }
          }
          mxDestroyArray(ffff);
          mxDestroyArray(ffff1);
    }
     
//y = y(m2+1:nn-m2);
return;
}
