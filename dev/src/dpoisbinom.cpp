#include <Rcpp.h>
#include <fftw3.h>

// [[Rcpp::export]]
Rcpp::NumericVector dpoisbinom(Rcpp::IntegerVector& nvec,
			       Rcpp::NumericVector& pp)
{

  using namespace Rcpp;

  fftw_complex *in, *out;
  fftw_plan p;
  
  int n = pp.size();
  int nn = nvec.size();
  int m = n + 1;

  
  NumericVector res(nn);
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);        
  std::complex<double> C(0.0,2);
  C = exp(C * 3.1415926535897 / ((double)m));
  double C_real = C.real();
  double C_imag = C.imag();


  //build input vector
  double tmp_real;
  double tmp_imag;
  std::complex<double> temp;
  std::complex<double> f(1.,0.0);
  in[0][0] = 1.0;
  in[0][1] = 0.0;

  
  int halfn = ceil((n) / 2) + 1; 
  for (std::size_t i = 1; i <= halfn; ++i){
    temp = 1.;
    tmp_real = f.real();
    tmp_imag = f.imag();
    f.real() = tmp_real * C_real - tmp_imag * C_imag;
    f.imag() = tmp_imag * C_real + tmp_real * C_imag;
    for(std::size_t j = 0; j < n; ++j){
      temp *= (1. + (f - 1.) * pp[j]);
    }
    
    in[i][0] = temp.real();
    in[i][1] = temp.imag();
    in[m-i][0] = temp.real();
    in[m-i][1] = - temp.imag();
  }


  //dft  
  p = fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  //form return object
  int kk; 
  for(std::size_t k = 0; k < nn; ++k)
  {
    kk = nvec[k];
    res[k] = out[kk][0] / m;
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return(res); 
}     
     
