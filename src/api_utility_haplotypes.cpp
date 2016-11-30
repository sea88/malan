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
std::vector<int> pedigree_get_father_haplotype(Rcpp::XPtr<Individual> individual) {
  return individual->get_father_haplotype();
}


//' @export
// [[Rcpp::export]]
int count_father_haplotype_occurrences_individuals(const List individuals, const Rcpp::IntegerVector haplotype) {
  int n = individuals.size();
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<int> h = Rcpp::as< std::vector<int> >(haplotype);
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    std::vector<int> indv_h = indv->get_father_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      count += 1;
    }
  }
  
  return count;
}


//' @export
// [[Rcpp::export]]
int count_father_haplotype_occurrences_pedigree(Rcpp::XPtr<Pedigree> pedigree, const Rcpp::IntegerVector haplotype, int generation_upper_bound_in_result = -1) {
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<int> h = Rcpp::as< std::vector<int> >(haplotype);

  std::vector<Individual*>* family = pedigree->get_all_individuals();
  
  for (auto dest : *family) {    
    int generation = dest->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    std::vector<int> dest_h = dest->get_father_haplotype();
    
    if (dest_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (dest_h == h) {
      count += 1;
    }    
  }
  
  return count;
}

