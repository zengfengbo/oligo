#include <Rcpp.h>
#include "oligotm.h"
#include <string.h>
using namespace Rcpp;
using namespace std;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}


double vacc3a(double age, bool female, bool ily=true){
  double p = 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily;
  p = p * (female ? 1.25 : 0.75);
  p = std::max(p, 0.0);
  p = std::min(p, 1.0);
  return p;
}

// [[Rcpp::export]]
NumericVector calcTm(CharacterVector seq,
               double dna_conc = 50,
               double salt_conc = 50,
               double divalent_conc = 0,
               double dntp_conc = 0,
               int tm_method = 0,
               int salt_corrections = 0
){
  int n = seq.size();
  vector< string > seqs = as< vector< string > >(seq);

  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    string seq_i = seqs[i];
    char * cstr = new char[seq_i.size() + 1];
    strcpy(cstr, seq_i.c_str());
    //tm = oligotm(seq, d, mv, dv, n, (tm_method_type) tm_santalucia, (salt_correction_type) salt_corrections);
    out[i] = oligotm(cstr, dna_conc, salt_conc, divalent_conc, dntp_conc, tm_method, salt_corrections);
    delete [] cstr;
  }

  return out;
}

NumericVector di_to_mo(NumericVector divalent, double dntp) {
  int n = divalent.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = divalent_to_monovalent(divalent[i], dntp);
  }

  return out;
}

NumericVector vacc3(NumericVector age, LogicalVector female,
                    LogicalVector ily) {
  int n = age.size();
  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = vacc3a(age[i], female[i], ily[i]);
  }

  return out;
}
