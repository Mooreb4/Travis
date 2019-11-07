#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "numrec.h"

#define TSUN 4.92569043916e-6                                           // mass to seconds conversion
#define PI 3.1415926535897931159979634685442            // Pi  
#define PI2 9.869604401089357992304940125905      // Pi^2 
#define TPI 6.2831853071795862319959269370884     // 2 Pi
#define RT2PI 2.5066282746310005024               // sqrt(2 Pi)
#define RTPI 1.772453850905516                    // sqrt(Pi)

// gcc -o qtransform qtransform.c -lm

double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX);
void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n);
double f_nwip(double *a, double *b, int n);
void Transform(double *a, double *times, double **tfmap, double Q, double fmin, double dfx, int n, int m);


int main(int argc, char **argv)
{
  int i, j, k, M, N;
  double SNR;
  double junk, Tobs, fix, f, t;
	double fmax, fmin, dfx, Q;
  double *times, *data;
  double **tfmap, **tfavH1, **tfavL1;
	
  char filename[1024];

  int n;

  FILE *in;
  FILE *out;
  FILE *samp;
	
  // Number of samples in time series
  N = 32768*2;
	
  // Choose the range of the spectrogram
  fmax = 100;
  fmin = 0.1;

  // Set the frequency resolution
  dfx = 0.1;
    
  // Set the Q of the transform
  Q = atof(argv[1]);
    
  M = (int)((fmax-fmin)/dfx);
	
  tfmap = dmatrix(0,M-1,0,N-1);
  tfavH1 = dmatrix(0,M-1,0,N-1);
  tfavL1 = dmatrix(0,M-1,0,N-1);
	
  times = dvector(0,N-1);
  data = dvector(0,N-1);
    
  // In this example we read in a data file with the format t, h(t)
        
  in = fopen("H1_inject.dat","r");
  for(i=0; i< N; i++)
  {
	  fscanf(in,"%lf%lf",&times[i], &data[i]);
  }
	
  // FFT
  Tobs = times[N-1]-times[0];
  fix = sqrt(Tobs)/((double)N);
  drealft(data-1,N,1);
	
  for(i = 0; i < N; i++) data[i] *= fix;
	
  out = fopen("wavefH1.dat","w");
  for(i = 1; i < N/2; i++)
  {
	  f = (double)(i)/Tobs;  
	  fprintf(out ,"%e %e\n", f, sqrt(data[2*i+1]*data[2*i+1]+data[2*i]*data[2*i]));
  }
  fclose(out);
    
  Transform(data, times, tfmap, Q, fmin, dfx, N, M);
	
  out = fopen("transformH1.dat","w");
  samp = fopen("freqsamp.dat","w");
	for(j = 0; j < M; j++)
	{
		f = fmin + dfx*(double)(j);
        fprintf(samp,"%e\n",f);
		
		for(i = 0; i < N; i++)
		{
			t = times[i];
			
			///fprintf(out,"%e %e %e\n", t, f, tfmap[j][i]);
            
            // change how it's printed to get along with python plotting
            fprintf(out,"%e ",tfmap[j][i]);
		}
		
		fprintf(out,"\n");
	}
	fclose(out);
    fclose(samp);
    
    samp = fopen("timesamp.dat","w");
    for(i = 0; i < N; i++){
        fprintf(samp,"%e\n",times[i]);
    }
    fclose(samp);
    
	
	

}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v=0;
	
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) fprintf(stderr,"allocation failure in dvector()");
    return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

/* ********************************************************************************** */
/*																					  */
/*                                 Fourier Routines                                   */
/*																					  */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign)
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi, swap;
	
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

void drealft(double data[], unsigned long n, int isign)
{
    void dfour1(double data[], unsigned long nn, int isign);
    unsigned long i,i1,i2,i3,i4,np3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;
	
    theta=3.141592653589793/(double) (n>>1);
    if (isign == 1) {
        c2 = -0.5;
        dfour1(data,n>>1,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        dfour1(data,n>>1,-1);
    }
}

void Transform(double *a, double *times, double **tfmap, double Q, double fmin, double dfx, int n, int m)
{
    double max, min;
    int i, j, k, p;
    int index;
    double tmp, mx;
	double f, t, dt, df, x;
    double *AC, *AF;
    double *corr, *b;
	double *params;
    double psh;
	double bmag;
	double Tobs;
	
	Tobs = times[n-1]-times[0];
	
	// [0] t0 [1] f0 [2] Q [3] Amp [4] phi
    
    params= dvector(0,4);
	
	params[0] = 0.0;
	params[2] = Q;
	params[3] = 1.0;
	params[4] = 0.0;
	
	dt = Tobs/(double)n;
    df = 1.0/Tobs;
	
	AC=dvector(0,n-1);  AF=dvector(0,n-1);
    corr = dvector(0,n-1);
	
	b = dvector(0,n-1);

	
	for(j = 0; j < m; j++)
	{
		
		f = fmin + dfx*(double)(j);
		
		// printf("%f\n", f);
		
		params[1] = f;
		
		SineGaussianF(b, params, Tobs, n);
		
		bmag = sqrt(f_nwip(b, b, n));
	    
		for(i = 0; i < n; i++) corr[i] = 0.0;
		
		phase_blind_time_shift(AC, AF, a, b, n);
		
		for(i = 0; i < n; i++) corr[i] += sqrt(AC[i]*AC[i]+AF[i]*AF[i])/bmag;
		
		for(i = 0; i < n; i++)
		{
		 tfmap[j][i] = corr[i];
		}

		
    }
	
	
	
	
    
    free_dvector(AC,0,n-1);  free_dvector(AF,0,n-1);
    free_dvector(corr, 0,n-1); 
	free_dvector(b, 0,n-1); 
    
}






void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
	int nb2, i, l, k, j;
	int imax, imin;
	
	nb2 = n / 2;
	
	corr[0] = 0.0;
	corr[1] = 0.0;
	corrf[0] = 0.0;
	corrf[1] = 0.0;
	
	for (i=1; i < nb2; i++)
    {
		l=2*i;
       k=l+1;

		corr[l]	= ( data1[l]*data2[l] + data1[k]*data2[k]);
		corr[k]	= ( data1[k]*data2[l] - data1[l]*data2[k]);
		corrf[l] = ( data1[l]*data2[k] - data1[k]*data2[l]);
		corrf[k] = ( data1[k]*data2[k] + data1[l]*data2[l]);

    }

  drealft(corr-1, n, -1);
  drealft(corrf-1, n, -1);

}


void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmax, fmin, fac;
    double phi, f, t, x, y, z, zold;
    double tau;
	int imin, imax;
    
    int i, id, N;
	
    t0 = sigpar[0];
    f0 = sigpar[1];
    Q = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
	
    tau = Q/(TPI*f0);
    
    fmax = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmin = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    fac = sqrt(Tobs);
    
	imin = (int)(fmin*Tobs);
	imax = (int)(fmax*Tobs);
	
	if(imin < 0) imin = 0;
	if(imax > NMAX/2) imax = NMAX/2;
    
    for(i = 0; i < NMAX/2; i++)
    {
		hs[2*i] = 0.0;
		hs[2*i+1] = 0.0;
		
		if(i > imin && i < imax)
		{
			f = (double)(i)/Tobs;
			sf = (Amp/fac)*RTPI/2.0*tau*exp(-PI*PI*tau*tau*(f-f0)*(f-f0));
			sx = exp(-Q*Q*f/f0);
			hs[2*i] = sf*(cos(TPI*f*t0-phi)+sx*cos(TPI*f*t0+phi));
			hs[2*i+1] = sf*(sin(TPI*f*t0-phi)+sx*sin(TPI*f*t0+phi));
		}
		
		// printf("%d %e %e\n", i, hs[2*i], hs[2*i+1]);
    }
	
    
}



double f_nwip(double *a, double *b, int n)
{
	int i, j, k;
	double arg, product;
	double test;
	double ReA, ReB, ImA, ImB;
	
	arg = 0.0;
	for(i=1; i<n/2; i++) 
    {
		j = i * 2;
		k = j + 1;
		ReA = a[j]; ImA = a[k];
		ReB = b[j]; ImB = b[k];
		product = ReA*ReB + ImA*ImB; 
		arg += product;										
    }
	
	return(4.0*arg);
	
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v;
	
    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) fprintf(stderr,"allocation failure in ivector()");
    return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;
	
    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
	
    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if (!m[nrl]) fprintf(stderr,"allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
	
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
    /* return pointer to array of pointers to rows */
    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    double ***t;
	
    /* allocate pointers to pointers to rows */
    t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
    if (!t) fprintf(stderr,"allocation failure 1 in d3tensor()");
    t += NR_END;
    t -= nrl;
	
    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
    if (!t[nrl]) fprintf(stderr,"allocation failure 2 in d3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
	
    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
    if (!t[nrl][ncl]) fprintf(stderr,"allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
	
    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }
	
    /* return pointer to array of pointers to rows */
    return t;
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float d3tensor allocated by d3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}
