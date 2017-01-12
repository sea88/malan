/*
Previous functions, unused for now
*/

/*

//' Get number of children for each individual in the population
//' 
//' @export
// [[Rcpp::export]]
IntegerMatrix get_number_of_children(Rcpp::XPtr<Population> population, bool progress = true) {
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  size_t N = pop.size();
  size_t i = 0;
  IntegerMatrix children_count(N, 2);
  colnames(children_count) = CharacterVector::create("pid", "boys_n");

  Progress p(N, progress);
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    Individual* ind = it->second;
    children_count(i, 0) = ind->get_pid();
    
    std::vector<Individual*>* children = ind->get_children();
    
    for (auto &child : (*children)) {
      children_count(i, 1) += 1;
    }
    
    ++i;
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      return children_count;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  if (i != N) {
    stop("Expected N indviduals...");
  }
  
  return children_count;
}
*/


/*
//' @export
// [[Rcpp::export]]
std::map<int, int> meioses_distribution(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  Pedigree* ped = i->get_pedigree();
  std::vector<Individual*>* family = ped->get_all_individuals();
  std::map<int, int> tab;
  
  for (auto dest : *family) {    
    int dist = i->meiosis_dist_tree(dest);
    tab[dist] += 1;
  } 
  
  return tab;
}
*/


/*
//' @export
// [[Rcpp::export]]
std::map<int, std::map<int, int> > meioses_generation_distribution_OLD(Rcpp::XPtr<Individual> individual, int generation_upper_bound = -1) {  
  Individual* i = individual;
  
  Pedigree* ped = i->get_pedigree();
  std::vector<Individual*>* family = ped->get_all_individuals();
  std::map<int, std::map<int, int> > tab;
  
  for (auto dest : *family) {    
    int generation = dest->get_generation();
    
    if (generation_upper_bound != -1 && generation > generation_upper_bound) {
      continue;
    }
    
    int dist = i->meiosis_dist_tree(dest);

    (tab[generation])[dist] += 1;    
  } 
  
  return tab;
}
*/

/*

//' get pids in pedigree with certain criteria
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pids_in_pedigree_criteria(Rcpp::XPtr<Pedigree> ped, bool must_be_alive, bool use_birth_year, int birth_year_min, int birth_year_max) {    
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  IntegerVector res;

  for (auto ind : *inds) {   
    bool skip = false;
    
    if (must_be_alive && ind->get_alive_status() == false) {
      skip = true;
    }
    
    if (!skip && use_birth_year) {
      int birth_year = ind->get_birth_year();
      
      if (birth_year < birth_year_min || birth_year > birth_year_max) {
        skip = true;
      }
    }
    
    if (skip) {
      continue;
    }
  
    res.push_back(ind->get_pid());
  } 
  
  return res;
}
*/

/*


//[[Rcpp::export]]
int meiosis_dist_tree(Rcpp::XPtr<Individual> src, Rcpp::XPtr<Individual> dest) {
  Individual* i1 = src;
  Individual* i2 = dest;
  int dist = i1->meiosis_dist_tree(i2);
  
  if (dist == -1) {
    Rcpp::stop("Individuals were not in the same pedigree");
  }
  
  return dist;
}

//' @export
//[[Rcpp::export]]
IntegerMatrix meiosis_dist_tree_matrix(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  size_t n = inds->size();
  
  CharacterVector nms(n);
  IntegerMatrix res(n, n);
  int i = 0;
  
  for (size_t i = 0; i < (n-1); ++i) {
    Individual* i1 = inds->at(i);
    nms[i] = std::to_string(i1->get_pid());
    res(i,i) = 0;
    
    for (size_t j = i+1; j < n; ++j) {
      Individual* i2 = inds->at(j);
      
      int dist = i1->meiosis_dist_tree(i2);
      res(i,j) = dist;
      res(j,i) = dist;
    }
  }
  
  nms[n-1] = std::to_string(inds->at(n-1)->get_pid());
  
  rownames(res) = nms;
  colnames(res) = nms;
  
  return res;
}

*/



/*
 

 
 */













