#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List indices_in_mixture(Rcpp::IntegerMatrix haplotypes, Rcpp::IntegerVector H1, Rcpp::IntegerVector H2) { 
  size_t N = haplotypes.nrow();
  
  Rcpp::List res;
  
  if (N == 0) {
    return res;
  }

  // mainly count wanted, but indices are good for debuggin
  Rcpp::IntegerVector res_in_mixture;
  Rcpp::IntegerVector res_H1;
  Rcpp::IntegerVector res_H2;
  
  size_t loci = haplotypes.ncol();

  for (size_t i = 0; i < N; ++i) {
    Rcpp::IntegerVector h = haplotypes(i, Rcpp::_);
    
    bool in_mixture = true;
    bool match_H1 = true; // faster than Rcpp equal/all sugar 
    bool match_H2 = true;
    
    for (size_t locus = 0; locus < loci; ++locus) {
      if (in_mixture && (h[locus] != H1[locus]) && (h[locus] != H2[locus])) {
        in_mixture = false;
      }
      
      if (match_H1 && (h[locus] != H1[locus])) {
        match_H1 = false;
      }
      
      if (match_H2 && (h[locus] != H2[locus])) {
        match_H2 = false;
      }
      
      // if neither have a chance, just stop
      if (!in_mixture && !match_H1 && !match_H2) {
        break;
      }
    }
    
    if (in_mixture) {
      res_in_mixture.push_back(i + 1); // R indexing
    }
    
    if (match_H1) {
      res_H1.push_back(i + 1); // R indexing
    }
    
    if (match_H2) {
      res_H2.push_back(i + 1); // R indexing
    }
  }
  
  res["in_mixture"] = res_in_mixture;
  res["match_H1"] = res_H1;
  res["match_H2"] = res_H2;

  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::List pedigree_get_haplotypes_pids(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids) {  
  size_t N = pids.size();
  Rcpp::List haps(N);

  for (size_t i = 0; i < N; ++i) {
    Individual* indv = population->get_individual(pids[i]);
    haps(i) = indv->get_haplotype();
  }

  return haps;
}
 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix individuals_get_haplotypes(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals) {   
  size_t n = individuals.size();
 
  if (n <= 0) {
    Rcpp::IntegerMatrix empty_haps(0, 0);
    return empty_haps;
  }
 
  size_t loci = individuals[0]->get_haplotype().size();

  if (loci <= 0) {
    Rcpp::stop("Expected > 0 loci");
    Rcpp::IntegerMatrix empty_haps(0, 0);
    return empty_haps;
  }

  Rcpp::IntegerMatrix haps(n, loci);

  for (size_t i = 0; i < n; ++i) {
    std::vector<int> hap = individuals[i]->get_haplotype();

    if (hap.size() != loci) {
      Rcpp::stop("Expected > 0 loci for all haplotypes");
      Rcpp::IntegerMatrix empty_haps(0, 0);
      return empty_haps;
    }
    
    Rcpp::IntegerVector h = Rcpp::wrap(hap);
    haps(i, Rcpp::_) = h;
  }

  return haps;
}

//' @export
// [[Rcpp::export]]
void pedigree_populate_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, Rcpp::NumericVector mutation_rates) {  
  Pedigree* p = ped;
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);

  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
    
  ped->populate_haplotypes(loci, mut_rates);
}

//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, Rcpp::NumericVector mutation_rates, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (loci != mut_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_haplotypes(loci, mut_rates);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' @export
// [[Rcpp::export]]
std::vector<int> get_haplotype(Rcpp::XPtr<Individual> individual) {
  return individual->get_haplotype();
}


//' @export
// [[Rcpp::export]]
int count_haplotype_occurrences_individuals(const Rcpp::List individuals, const Rcpp::IntegerVector haplotype) {
  int n = individuals.size();
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<int> h = Rcpp::as< std::vector<int> >(haplotype);
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    std::vector<int> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      count += 1;
    }
  }
  
  return count;
}

// MIKL: DEPRECATED FOR pedigree_matches_with_mutations?
//
// individuals typically only generation 0 to 2...
// 
// meiosis distance only well defined for pedigrees
//
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector meiosis_dist_haplotype_matches_individuals(const Rcpp::XPtr<Individual> suspect, const Rcpp::List individuals) {
  std::vector<int> h = suspect->get_haplotype();
  
  int n = individuals.size();
  int loci = h.size();
  Rcpp::IntegerVector meiosis_dists;
  
  for (int i = 0; i < n; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    
    if (indv->get_pedigree_id() != suspect->get_pedigree_id()) {
      continue; // no meiosis distance for individuals in different pedigrees
    }
    
    std::vector<int> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      int dist = suspect->meiosis_dist_tree(indv);
      
      if (dist == -1) {
        meiosis_dists.push_back(R_PosInf);        
      } else {    
        meiosis_dists.push_back(dist);
      }
    }
  }
  
  return meiosis_dists;
}

// There are count_haplotype_occurrences_pedigree matches. 
// This gives details on meiotic distance and the max L1 distance of the haplotypes on the path between the suspect and the matching individual in the pedigree. Some of these matches may have (back)mutations between in between them.
//
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(const Rcpp::XPtr<Individual> suspect, int generation_upper_bound_in_result = -1) {
  const std::vector<int> h = suspect->get_haplotype();

  const Pedigree* pedigree = suspect->get_pedigree();
  const int suspect_pedigree_id = suspect->get_pedigree_id();
  const std::vector<Individual*>* family = pedigree->get_all_individuals();
  
  std::vector<int> meiosis_dists;
  std::vector<int> max_L1_dists;
  
  // includes suspect by purpose
  for (auto dest : *family) { 
    int generation = dest->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    // only considering within pedigree matches
    if (dest->get_pedigree_id() != suspect_pedigree_id) {
      continue;
    }
    
    std::vector<int> dest_h = dest->get_haplotype();
    
    if (dest_h.size() != h.size()) {
      Rcpp::stop("haplotype and dest_h did not have same number of loci");
    }
    
    if (dest_h == h) {
      std::vector<Individual*> path = suspect->calculate_path_to(dest);  
      int meiosis_dist = suspect->meiosis_dist_tree(dest);
      
      int meiosis_dist_from_path = path.size() - 1; // n vertices means n-1 edges (tree)
      //Rcpp::Rcout << ">> path from " << suspect->get_pid() << " to " << dest->get_pid() << " has length = " << meiosis_dist_from_path << " and meioses = " << meiosis_dist << (meiosis_dist_from_path == meiosis_dist ? " ok" : " ERROR") << ": " << std::endl;
      
      int max_L1 = 0;
      
      //Rcpp::Rcout << "  ";
      
      for (auto intermediate_node : path) { 
        //Rcpp::Rcout << intermediate_node->get_pid();
        
        int d = suspect->get_haplotype_L1(intermediate_node);
        
        if (d > max_L1) {
          max_L1 = d;
          //Rcpp::Rcout << "!";
        }
        
        //Rcpp::Rcout << " ";
      }
      
      //Rcpp::Rcout << std::endl;      
      
      if (meiosis_dist == -1) {
        Rcpp::stop("Cannot occur in pedigree!");
      }
      
      meiosis_dists.push_back(meiosis_dist);
      max_L1_dists.push_back(max_L1);
    }
  }
  
  size_t n = meiosis_dists.size();
  
  Rcpp::IntegerMatrix matches(n, 2);
  colnames(matches) = Rcpp::CharacterVector::create("meioses", "max_L1");
  
  for (size_t i = 0; i < n; ++i) {
    matches(i, 0) = meiosis_dists[i];
    matches(i, 1) = max_L1_dists[i];
  }
  
  return matches;
}

//' @export
// [[Rcpp::export]]
int count_haplotype_occurrences_pedigree(Rcpp::XPtr<Pedigree> pedigree, const Rcpp::IntegerVector haplotype, int generation_upper_bound_in_result = -1) {
  int loci = haplotype.size();
  int count = 0;
  
  std::vector<int> h = Rcpp::as< std::vector<int> >(haplotype);

  std::vector<Individual*>* family = pedigree->get_all_individuals();

  for (auto dest : *family) {    
    int generation = dest->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    std::vector<int> dest_h = dest->get_haplotype();
    
    if (dest_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (dest_h == h) {
      count += 1;
    }    
  }
  
  return count;
}

