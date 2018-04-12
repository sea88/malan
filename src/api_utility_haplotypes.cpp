#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.hpp"
#include "api_utility_individual.hpp"

//' Mixture information about 2 persons' mixture of donor1 and donor2.
//' 
//' @param individuals Individuals to consider as possible contributors and thereby get information from.
//' @param donor1 Contributor1/donor 1
//' @param donor2 Contributor2/donor 2
//' @return A list with mixture information about the mixture \code{donor1}+\code{donor2}+\code{donor3} from \code{individuals}
//' 
//' @seealso \code{\link{mixture_info_by_individuals_3pers}}
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List mixture_info_by_individuals(const Rcpp::List individuals, Rcpp::XPtr<Individual>& donor1, Rcpp::XPtr<Individual>& donor2) { 
  size_t N = individuals.size();
  
  Rcpp::List res;
  
  if (N == 0) {
    return res;
  }

  // mainly count wanted, but indices are good for debugging
  Rcpp::IntegerVector res_comp_with_mixture;
  Rcpp::List res_comp_with_mixture_dists;
  Rcpp::IntegerVector res_match_donor1;
  Rcpp::IntegerVector res_match_donor2;
  Rcpp::IntegerVector res_others_included;
  
  std::vector<int> H1 = donor1->get_haplotype();
  std::vector<int> H2 = donor2->get_haplotype();
  
  size_t loci = H1.size();
  
  if (H2.size() != loci) {
    Rcpp::stop("H2.size() != H1.size()");
  }
  
  size_t loci_not_matching = 0;
  
  for (size_t locus = 0; locus < loci; ++locus) {
    if (H1[locus] != H2[locus]) {
      loci_not_matching += 1;
    }
  }

  for (size_t i = 0; i < N; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    std::vector<int> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("indv_h.size() != H1.size()");
    }
    
    bool in_mixture = true;
    bool match_H1 = true; // faster than Rcpp equal/all sugar 
    bool match_H2 = true;
    
    for (size_t locus = 0; locus < loci; ++locus) {
      if (in_mixture && (indv_h[locus] != H1[locus]) && (indv_h[locus] != H2[locus])) {
        in_mixture = false;
      }
      
      if (match_H1 && (indv_h[locus] != H1[locus])) {
        match_H1 = false;
      }
      
      if (match_H2 && (indv_h[locus] != H2[locus])) {
        match_H2 = false;
      }
      
      // if neither have a chance, just stop
      if (!in_mixture && !match_H1 && !match_H2) {
        break;
      }
    }
    
    int pid = indv->get_pid();
    
    if (in_mixture) {
      res_comp_with_mixture.push_back(pid); // R indexing
      
      //int dist_donor1 = donor1->calculate_path_to(indv);
      //int dist_donor2 = donor2->calculate_path_to(indv);      
      int dist_donor1 = donor1->meiosis_dist_tree(indv);
      int dist_donor2 = donor2->meiosis_dist_tree(indv);
      
      Rcpp::List r = Rcpp::List::create(
        Rcpp::Named("indv_pid") = pid,
        //Rcpp::Named("pid_donor1") = donor1->get_pid(),
        //Rcpp::Named("pid_donor2") = donor2->get_pid(),
        Rcpp::Named("dist_donor1") = dist_donor1,
        Rcpp::Named("dist_donor2") = dist_donor2);   
        
      res_comp_with_mixture_dists.push_back(r);
      
      if (match_H1) {
        res_match_donor1.push_back(pid);
      }
      
      if (match_H2) {
        res_match_donor2.push_back(pid);
      }
      
      if (!match_H1 && !match_H2) {
        res_others_included.push_back(pid);
      }
    }    
  }
  
  res["pids_included_in_mixture"] = res_comp_with_mixture;
  res["pids_included_in_mixture_info"] = res_comp_with_mixture_dists;
  res["pids_matching_donor1"] = res_match_donor1;
  res["pids_matching_donor2"] = res_match_donor2;
  res["pids_others_included"] = res_others_included;
  res["pids_donor12_meiotic_dist"] = donor1->meiosis_dist_tree(donor2);
  res["donor1_family_info"] = get_family_info(donor1);
  res["donor2_family_info"] = get_family_info(donor2);
  res["donor1_profile"] = H1;
  res["donor2_profile"] = H2;
  res["donor1_pid"] = donor1->get_pid();
  res["donor2_pid"] = donor2->get_pid();
  res["loci_not_matching"] = loci_not_matching;

  return res;
}













//' Mixture information about 3 persons' mixture of donor1, donor2 and donor3.
//' 
//' @inherit mixture_info_by_individuals
//' @param donor3 Contributor2/donor 3
//' 
//' @seealso \code{\link{mixture_info_by_individuals}}
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List mixture_info_by_individuals_3pers(const Rcpp::List individuals, 
    Rcpp::XPtr<Individual>& donor1, 
    Rcpp::XPtr<Individual>& donor2, 
    Rcpp::XPtr<Individual>& donor3) { 
    
  size_t N = individuals.size();
  
  Rcpp::List res;
  
  if (N == 0) {
    return res;
  }

  // mainly count wanted, but indices are good for debugging
  Rcpp::IntegerVector res_comp_with_mixture;
  Rcpp::IntegerVector res_match_donor1;
  Rcpp::IntegerVector res_match_donor2;
  Rcpp::IntegerVector res_match_donor3;
  Rcpp::IntegerVector res_others_included;
  
  std::vector<int> H1 = donor1->get_haplotype();
  std::vector<int> H2 = donor2->get_haplotype();
  std::vector<int> H3 = donor3->get_haplotype();
  
  size_t loci = H1.size();
  
  if (H2.size() != loci) {
    Rcpp::stop("H2.size() != H1.size()");
  }

  if (H3.size() != loci) {
    Rcpp::stop("H3.size() != H1.size()");
  }

  for (size_t i = 0; i < N; ++i) {
    Rcpp::XPtr<Individual> indv = individuals[i];
    std::vector<int> indv_h = indv->get_haplotype();
    
    if (indv_h.size() != loci) {
      Rcpp::stop("indv_h.size() != H1.size()");
    }
    
    bool in_mixture = true;
    bool match_H1 = true; // faster than Rcpp equal/all sugar 
    bool match_H2 = true;
    bool match_H3 = true;
    
    for (size_t locus = 0; locus < loci; ++locus) {
      if (in_mixture && 
            (indv_h[locus] != H1[locus]) && 
            (indv_h[locus] != H2[locus]) && 
            (indv_h[locus] != H3[locus])) {
        in_mixture = false;
      }
      
      if (match_H1 && (indv_h[locus] != H1[locus])) {
        match_H1 = false;
      }
      
      if (match_H2 && (indv_h[locus] != H2[locus])) {
        match_H2 = false;
      }
      
      if (match_H3 && (indv_h[locus] != H3[locus])) {
        match_H3 = false;
      }
      
      // if neither have a chance, just stop
      if (!in_mixture && !match_H1 && !match_H2 && !match_H3) {
        break;
      }
    }
    
    int pid = indv->get_pid();
    
    if (in_mixture) {
      res_comp_with_mixture.push_back(pid); // R indexing
            
      if (match_H1) {
        res_match_donor1.push_back(pid);
      }
      
      if (match_H2) {
        res_match_donor2.push_back(pid);
      }
      
      if (match_H3) {
        res_match_donor3.push_back(pid);
      }
      
      if (!match_H1 && !match_H2 && !match_H3) {
        res_others_included.push_back(pid);
      }
    }    
  }
  
  res["pids_included_in_mixture"] = res_comp_with_mixture;
  res["pids_matching_donor1"] = res_match_donor1;
  res["pids_matching_donor2"] = res_match_donor2;
  res["pids_matching_donor3"] = res_match_donor3;
  res["pids_others_included"] = res_others_included;
  res["donor1_profile"] = H1;
  res["donor2_profile"] = H2;
  res["donor3_profile"] = H3;
  res["donor1_pid"] = donor1->get_pid();
  res["donor2_pid"] = donor2->get_pid();
  res["donor3_pid"] = donor3->get_pid();

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
 
//' Get haplotype matrix from list of individuals
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

//' Unbounded, founders all zero
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                       int loci, 
                                       Rcpp::NumericVector mutation_rates, 
                                       bool progress = true) {
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

//' Unbounded, custom founders
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes_custom_founders(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                       Rcpp::NumericVector mutation_rates,
                                       Rcpp::Nullable<Rcpp::Function> get_founder_haplotype = R_NilValue,
                                       bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  
  if (get_founder_haplotype.isNull()) {
    Rcpp::stop("get_founder_haplotype must not be NULL");
  }  
  
  Rcpp::Function g_founder_hap = Rcpp::as<Rcpp::Function>(get_founder_haplotype);

  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_haplotypes_custom_founders(mut_rates, g_founder_hap);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}

//' Bounded, custom founders
//' 
//' Populate haplotypes such that they are all on-ladder
//' 
//' @param get_founder_haplotype has no default as it is not know in advance how many loci there are and what the ladder is; see \code{\link{generate_get_founder_haplotype}}
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_haplotypes_ladder_bounded(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                                      Rcpp::NumericVector mutation_rates, 
                                                      Rcpp::IntegerVector ladder_min,
                                                      Rcpp::IntegerVector ladder_max,
                                                      Rcpp::Nullable<Rcpp::Function> get_founder_haplotype = R_NilValue,
                                                      bool progress = true) {
  //https://stackoverflow.com/questions/36992627/can-rcppfunction-be-null
  
  std::vector<Pedigree*> peds = (*pedigrees);

  std::vector<double> mut_rates = Rcpp::as< std::vector<double> >(mutation_rates);
  std::vector<int> lad_min = Rcpp::as< std::vector<int> >(ladder_min);
  std::vector<int> lad_max = Rcpp::as< std::vector<int> >(ladder_max);
  

  if (mutation_rates.size() != lad_min.size()) {
    Rcpp::stop("mutation_rates and ladder_min must have same length");
  }

  if (mutation_rates.size() != lad_max.size()) {
    Rcpp::stop("mutation_rates and ladder_max must have same length");
  }

  if (get_founder_haplotype.isNull()) {
    Rcpp::stop("get_founder_haplotype must not be NULL");
  }  
  
  Rcpp::Function g_founder_hap = Rcpp::as<Rcpp::Function>(get_founder_haplotype);
    
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_haplotypes_ladder_bounded(mut_rates, lad_min, lad_max, g_founder_hap);
    
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
  std::vector<int> pids;
  
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
      pids.push_back(dest->get_pid());
    }
  }
  
  size_t n = meiosis_dists.size();
  
  Rcpp::IntegerMatrix matches(n, 3);
  colnames(matches) = Rcpp::CharacterVector::create("meioses", "max_L1", "pid");
  
  for (size_t i = 0; i < n; ++i) {
    matches(i, 0) = meiosis_dists[i];
    matches(i, 1) = max_L1_dists[i];
    matches(i, 2) = pids[i];
  }
  
  return matches;
}
  
//' @export
// [[Rcpp::export]]
int meiotic_dist(Rcpp::XPtr<Individual> ind1, Rcpp::XPtr<Individual> ind2) {
  return ind1->meiosis_dist_tree(ind2);
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

