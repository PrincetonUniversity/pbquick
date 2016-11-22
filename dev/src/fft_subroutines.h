#ifndef __FFT_SUBROUTINES_H_INCLUDED__
#define __FFT_SUBROUTINES_H_INCLUDED__

static void fftmx(double *a, double *b, int ntot, int n, int nspan, int isn,
                  int m, int kt, double *at, double *ck, double *bt, double *sk,
                  int *np, int *nfac);

void fft_factor(int n, int *pmaxf, int *pmaxp);

Rboolean fft_work(double *a, double *b, int nseg, int n, int nspn, int isn,
                  double *work, int *iwork);

#endif
