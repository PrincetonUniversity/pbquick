#include <Rcpp.h>
#include <R_ext/Utils.h> // for findInterval();
#include <fftw3.h>


void dft_pmf(fftw_complex* out, int m,  Rcpp::NumericVector& pp);

// [[Rcpp::export]]
Rcpp::NumericVector dpoisbinom(Rcpp::IntegerVector& invec,
			       Rcpp::NumericVector& pp,
			       bool log_d = false)
{
  int m = pp.size() + 1;
  int nn = invec.size();
  fftw_complex* out;  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);

  //dft
  dft_pmf(out, m, pp);

 
  //form return object
  Rcpp::NumericVector res(nn);
  double scale = 1.0 / m;
  int kk;
  for(std::size_t k = 0; k < nn; ++k)
    {
      kk = invec[k];
      res[k] = out[kk][0] * scale;
    }
  
  //Destroy dft object
  fftw_free(out);
  
  if(log_d)
    return(log(res));
  else
    return(res);
}

// [[Rcpp::export]]
Rcpp::NumericVector ppoisbinom(Rcpp::IntegerVector& invec,
			       Rcpp::NumericVector& pp,
			       bool lower_tail = true,
			       bool log_p = false)
{
  int m = pp.size() + 1;
  int nn = invec.size();
  int max_q = max(invec) + 1; // maximum of quantiles plus one
  fftw_complex* out;  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);

  //dft
  dft_pmf(out, m, pp);

  // form cumulative probabilities
  Rcpp::NumericVector csum(max_q);
  double scale = 1.0 / m;
  csum[0] = out[0][0] * scale;
  int k = 1;    
  do
    {
      csum[k] = out[k][0] * scale + csum[k - 1];
      ++k;
    }
  while(k < max_q);

  //Destroy dft object
  fftw_free(out);
  
  //form return object
  Rcpp::NumericVector res(nn);
  int kk;
  for(std::size_t k = 0; k < nn; ++k)
    {
      kk = invec[k];
      res[k] = csum[kk];
    }

  if(!lower_tail)
    res = 1.0 - res;
  
  if(log_p)
    return(log(res));
  else
    return(res);
}


// [[Rcpp::export]]
Rcpp::IntegerVector qpoisbinom(Rcpp::NumericVector& invec,
			       Rcpp::NumericVector& pp,
			       bool lower_tail = true,
			       bool log_p = false)
{
  if (log_p) invec = exp(invec);
  
  int nn = invec.size();
  int m = pp.size() + 1;
  fftw_complex* out;  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);

  //dft
  dft_pmf(out, m, pp);

  //form cumulative probabilities
  Rcpp::NumericVector csum(m);
  double scale = 1.0 / m;
  csum[0] = out[0][0] * scale;
  int k = 1;    
  do
    {
      csum[k] = out[k][0] * scale + csum[k - 1];
      ++k;
    }
  while(k < m);

  //Destroy dft object
  fftw_free(out);
  
  //sort keeping track of original order
  Rcpp::NumericVector s_invec = Rcpp::clone(invec).sort();
  Rcpp::IntegerVector order = Rcpp::match(s_invec, invec);

  //find interval on sorted vector, and form return object
  Rcpp::IntegerVector res(nn);
  int kk;
  int flag;
  int t_res;
  for(std::size_t k = 0; k < nn; ++k)
    {
      t_res = findInterval(&csum[0], csum.size(), s_invec[k], FALSE, FALSE, t_res, &flag);
      kk = order[k];
      res[kk-1] = t_res;
    }
  
  return(res);
    
}

void dft_pmf(fftw_complex* out, int m,  Rcpp::NumericVector& pp)
{
  int n = m - 1;
  fftw_complex* in;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);
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
    f.real(tmp_real * C_real - tmp_imag * C_imag);
    f.imag(tmp_imag * C_real + tmp_real * C_imag);
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
  fftw_destroy_plan(p);
  fftw_free(in);
}

