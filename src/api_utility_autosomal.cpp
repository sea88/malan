/**
 api_utility_haplotypes.cpp
 Purpose: Logic related to haplotypes.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>
#include <unordered_map>
#include <unordered_set>

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

//' Calculate conditional genotype cumulative probabilities with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//' 
//' @return Matrix: row i: conditional cumulative distribution of alleles given allele i
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_autosomal_genotype_conditional_cumdist(
    Rcpp::NumericVector allele_dist,
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
  
  Rcpp::NumericMatrix dists(alleles_count, alleles_count);

  for (int i = 0; i < alleles_count; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (i == j) { // homozyg
        double p = theta*ps[i] + (1.0-theta)*ps[i]*ps[i];
        dists(i, i) = p;
      } else { // hetegozyg
        double p = (1.0-theta)*ps[i]*ps[j];
        dists(i, j) = p;
        dists(j, i) = p;
      }
    }
  }
  
  // Multiple passes, but easier to follow:
  
  // Get row i to sum to 1; ps[i] = sum(dists(i, Rcpp::_))
  for (int i = 0; i < alleles_count; ++i) {
    Rcpp::NumericVector row = dists(i, Rcpp::_) / ps[i];
    
    // cumsum, for some reason Rcpp::cumsum doesn't work...
    Rcpp::NumericVector res(row.size());
    std::partial_sum(row.begin(), row.end(), res.begin(), std::plus<double>());
    dists(i, Rcpp::_) = res;
  }
  
  return dists;                                   
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

  // For drawing founder types ->
  const int alleles_count = allele_dist.size();
  const std::vector<double> allele_dist_theta = calc_autosomal_genotype_probs(allele_dist, theta);
  std::vector<double> allele_cumdist_theta(allele_dist_theta.size());
  std::partial_sum(allele_dist_theta.begin(), allele_dist_theta.end(), allele_cumdist_theta.begin(), std::plus<double>());
  // <- founder
  
  // For children's ->
  Rcpp::NumericMatrix cumdist_mat = calc_autosomal_genotype_conditional_cumdist(allele_dist, theta);
  
  if (cumdist_mat.nrow() != alleles_count) {
    Rcpp::stop("Unexpected error");
  }
  std::vector< std::vector<double> > cumdists(alleles_count);
  
  for (int i = 0; i < alleles_count; ++i) {
    Rcpp::NumericVector row_vec = cumdist_mat(i, Rcpp::_);
    std::vector<double> row = Rcpp::as< std::vector<double> >(row_vec);
    cumdists[i] = row;
  }
  // <- 
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_autosomal(cumdists, allele_cumdist_theta, alleles_count, mutation_rate);
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}



// boost::hash_combine
// https://stackoverflow.com/a/27952689/3446913
size_t hash_combine(size_t lhs, size_t rhs) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}

// https://stackoverflow.com/a/20602159/3446913
struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return hash_combine(x.first, x.second);
  }
};

// [[Rcpp::export]]
std::unordered_map<int, int> hash_colisions(int p) {
  std::unordered_map<int, int> tab;
  
  for (int i = 0; i < (p-1); ++i) {
    for (int j = (i+1); j < p; ++j) {
      int hash = hash_combine(i, j);
      tab[hash] += 1;
    }
  }
  
  return tab;
}




Rcpp::List estimate_theta_1subpop(const std::unordered_map<int, double>& allele_p,
                                  const std::unordered_map<std::pair<int, int>, double, pairhash>& genotype_p,
                                  const std::unordered_set<std::pair<int, int>, pairhash>& genotypes_unique) {
  Rcpp::List theta;
  
  // Loop over unique genotypes
  std::unordered_set<std::pair<int, int>, pairhash>::const_iterator it;
  int K = genotypes_unique.size();
  
  if (K == 1) {
    theta["estimate"] = NA_REAL;
    theta["error"] = true;
    theta["details"] = "Only one genotype observed";
    return theta;
  }
  
  int k = 0;
  arma::mat X(K, 1, arma::fill::none);
  arma::vec y(K, arma::fill::none);
  
  for (it = genotypes_unique.begin(); it != genotypes_unique.end(); ++it) {
    std::pair<int, int> geno = *it;
    int a1 = geno.first;
    int a2 = geno.second;
    
    // homozyg
    if (a1 == a2) {
      double p_i = allele_p.at(a1);
      double p_ii = genotype_p.at(geno);
      double p_i2 = p_i*p_i;
      X(k, 0) = p_i - p_i2;
      y(k) = p_ii - p_i2;
    } else {
      // heterozyg
      double p_i = allele_p.at(a1);
      double p_j = allele_p.at(a2);
      double p_ij = genotype_p.at(geno);
      double tmp = -2.0*p_i*p_j;
      X(k, 0) = tmp;
      y(k) = p_ij + tmp;
    }
    
    ++k;
  }
  
  // minimisze (Xb - y)^2 for b
  arma::mat Q, R;
  bool status = arma::qr_econ(Q, R, X);
  
  if (!status) {
    theta["estimate"] = NA_REAL;
    theta["error"] = true;
    theta["details"] = "Could not make QR decomposition";
  } else {
    arma::vec coef = arma::solve(R, Q.t() * y, arma::solve_opts::no_approx);
    
    if (coef[0] >= 0 && coef[0] <= 1) {
      theta["estimate"] = coef[0];
      theta["error"] = false;
      theta["details"] = "OK";
    } else {
      theta["estimate"] = coef[0];
      theta["error"] = true;
      theta["details"] = "Estimate outside range of (0, 1)";
    }
  }
  
  return theta;
}


void estimate_theta_1subpop_fill_containers(int a1,
                                            int a2,
                                            const double one_over_n,
                                            const double one_over_2n,
                                            std::unordered_map<int, double>& allele_p,
                                            std::unordered_map<std::pair<int, int>, double, pairhash>& genotype_p,
                                            std::unordered_set<std::pair<int, int>, pairhash>& genotypes_unique) {
  
  if (a2 < a1) {
    int tmp = a1;
    a1 = a2;
    a2 = tmp;
  }
  
  std::pair<int, int> geno = std::make_pair(a1, a2);
  genotypes_unique.insert(geno);
  
  genotype_p[geno] += one_over_n;
  
  if (a1 == a2) {
    allele_p[a1] += one_over_n; // 2*one_over_2n = one_over_n
  } else {
    allele_p[a1] += one_over_2n;
    allele_p[a2] += one_over_2n;
  }
}

//' Estimate theta from genetypes
//' 
//' Estimate theta for one subpopulation given a sample of genotypes.
//' 
//' @param x Matrix of genotypes: two columns (allele1 and allele2) and a row per individual
//' 
//' @return List:
//' * estimate: Vector of length 1 containing estimate of theta or NA if it could not be estimated
//' * error: true if an error happened, false otherwise
//' * details: contains description if an error happened
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_1subpop_sample(Rcpp::IntegerMatrix x) {
  int n = x.nrow();
  
  if (n <= 0) {
    Rcpp::stop("x cannot be empty");
  }
  
  if (x.ncol() != 2) {
    Rcpp::stop("x must have exactly two columns");
  }
  
  // Build count tables
  std::unordered_map<int, double> allele_p;
  std::unordered_map<std::pair<int, int>, double, pairhash> genotype_p;
  std::unordered_set<std::pair<int, int>, pairhash> genotypes_unique;
  
  double one_over_n = 1.0 / (double)n;
  double one_over_2n = 1.0 / (2.0 * (double)n);
  
  for (int i = 0; i < n; ++i) {
    int a1 = x(i, 0);
    int a2 = x(i, 1);
    
    estimate_theta_1subpop_fill_containers(a1, a2, one_over_n, one_over_2n, 
                                           allele_p, genotype_p, genotypes_unique);
  }
  
  Rcpp::List theta = estimate_theta_1subpop(allele_p, genotype_p, genotypes_unique);
  return theta;
}


//' Estimate theta from individuals
//' 
//' Estimate theta for one subpopulation given a sample of genotypes.
//' 
//' @param individuals Individuals to get haplotypes for.
//' 
//' @return List:
//' * estimate: Vector of length 1 containing estimate of theta or NA if it could not be estimated
//' * error: true if an error happened, false otherwise
//' * details: contains description if an error happened
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_1subpop_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {
  int n = individuals.size();
  
  if (n <= 0) {
    Rcpp::stop("No individuals given");
  }
  
  if (!(individuals[0]->is_haplotype_set())) {
    Rcpp::stop("Haplotypes not yet set");
  }
  
  int loci = individuals[0]->get_haplotype().size();
  
  if (loci != 2) {
    Rcpp::stop("Expected exactly 2 autosomal loci");
  }
  
  // Build count tables
  std::unordered_map<int, double> allele_p;
  std::unordered_map<std::pair<int, int>, double, pairhash> genotype_p;
  std::unordered_set<std::pair<int, int>, pairhash> genotypes_unique;
  
  double one_over_n = 1.0 / (double)n;
  double one_over_2n = 1.0 / (2.0 * (double)n);

  for (int i = 0; i < n; ++i) {
    Individual* individual = individuals[i];
    std::vector<int> hap = individual->get_haplotype();
    
    estimate_theta_1subpop_fill_containers(hap[0], hap[1], one_over_n, one_over_2n, 
                                           allele_p, genotype_p, genotypes_unique);
  }
  
  Rcpp::List theta = estimate_theta_1subpop(allele_p, genotype_p, genotypes_unique);
  
  return theta;
}


