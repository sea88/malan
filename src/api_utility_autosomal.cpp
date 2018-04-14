/**
 api_utility_haplotypes.cpp
 Purpose: Logic related to haplotypes.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.h"
#include "api_utility_individual.h"

//' Calculate genotype probabilities with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> calc_autosomal_genotype_probs(Rcpp::NumericVector allele_dist,
                                                  double theta) {
  
  if (any(allele_dist < 0).is_true() || any(allele_dist > 1).is_true()) {
    Rcpp::stop("allele_dist's elements must be between 0 and 1, both included");
  }
  
  if (theta < 0 || theta > 1) {
    Rcpp::stop("theta must be between 0 and 1, both included");
  }
  
  std::vector<double> ps = Rcpp::as< std::vector<double> >(allele_dist);
  double ps_sum = std::accumulate(ps.begin(), ps.end(), 0.0);
  const int alleles_count = ps.size();          
  
  // Normalisation
  for (int i = 0; i < alleles_count; ++i) {
    ps[i] = ps[i] / ps_sum;
  }
  
  std::vector<double> allele_dist_theta(alleles_count * (alleles_count + 1) / 2);
  int k = 0;
                                    
  for (int i = 0; i < alleles_count; ++i) {
    for (int j = 0; j <= i; ++j) {   
      if (i == j) { // homozyg
        allele_dist_theta[k] = theta*ps[i] + (1.0-theta)*ps[i]*ps[i];        
      } else { // hetegozyg
        allele_dist_theta[k] = (1.0-theta)*2.0*ps[i]*ps[j];
      }

      k++;
    }
  }
  
  return allele_dist_theta;                                   
}
             
//' Sample genotype with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//'
//' @export
// [[Rcpp::export]]
std::vector<int> sample_autosomal_genotype(Rcpp::NumericVector allele_dist,
                                           double theta) {
                                           
  const int alleles_count = allele_dist.size();
  const std::vector<double> allele_dist_theta = calc_autosomal_genotype_probs(allele_dist, theta);
  
  std::vector<double> allele_cumdist_theta(allele_dist_theta.size());
  std::partial_sum(allele_dist_theta.begin(), allele_dist_theta.end(), allele_cumdist_theta.begin(), std::plus<double>());
  
  std::vector<int> geno = draw_autosomal_genotype(allele_cumdist_theta, alleles_count);
  
  return geno;
}
                                      
//' Populate 1-locus autosomal DNA profile in pedigrees.
//' 
//' Populate 1-locus autosomal DNA profile from founder and down in all pedigrees.
//' Note, that only alleles from ladder is assigned and 
//' that all founders draw type randomly.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate haplotypes
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//' @param mutation_rate Mutation rate between 0 and 1 (both included)
//' @param progress Show progress
//'
//' @seealso [pedigrees_all_populate_haplotypes_custom_founders()] and 
//' [pedigrees_all_populate_haplotypes_ladder_bounded()].
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_autosomal(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                      Rcpp::NumericVector allele_dist,
                                      double theta,
                                      double mutation_rate,
                                      bool progress = true) {  
  std::vector<Pedigree*> peds = (*pedigrees);

  const int alleles_count = allele_dist.size();
  const std::vector<double> allele_dist_theta = calc_autosomal_genotype_probs(allele_dist, theta);

  std::vector<double> allele_cumdist_theta(allele_dist_theta.size());
  std::partial_sum(allele_dist_theta.begin(), allele_dist_theta.end(), allele_cumdist_theta.begin(), std::plus<double>());

  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_autosomal(allele_cumdist_theta, alleles_count, mutation_rate);
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

