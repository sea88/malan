#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.hpp"


//' @export
// [[Rcpp::export]]
void pedigree_populate_father_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, NumericVector mutation_rates) {  
  Pedigree* p = ped;
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);

  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
    
  ped->populate_father_haplotypes(loci, mut_rates);
}

//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_father_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, NumericVector mutation_rates, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_father_haplotypes(loci, mut_rates);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' @export
// [[Rcpp::export]]
List pedigree_get_father_haplotypes_pids(Rcpp::XPtr<Population> population, IntegerVector pids) {  
 
  size_t N = pids.size();
  List haps(N);
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = population->get_individual(pids[i]);
    haps(i) = indv->get_father_haplotype();
  }
  
  return haps;
}


