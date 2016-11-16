#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.hpp"

//[[Rcpp::export]]
void malan_test() {
  Rcout << "mikl was here 1324" << std::endl;
}


//[[Rcpp::export]]
int pop_size(Rcpp::XPtr<Population> population) {
  return population->get_population_size();
}


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


//' Get number of pedigrees
//' 
//' @export
// [[Rcpp::export]]
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  return pedigrees->size();
}

//' Get pedigree size
//' 
//' @export
// [[Rcpp::export]]
int pedigree_size(Rcpp::XPtr<Pedigree> ped) {  
  return ped->get_all_individuals()->size();
}

//' Get distribution of pedigree sizes
//' 
//' @export
//[[Rcpp::export]]
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  std::vector<Pedigree*>* peds = pedigrees;
  std::unordered_map<int, int> tab;
  
  for (auto it = peds->begin(); it != peds->end(); ++it) {
    tab[(*it)->get_all_individuals()->size()] += 1;
  }
  
  return tab;
}

//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index) {  
  std::vector<Pedigree*>* peds = pedigrees;
  Pedigree* p = peds->at(index);
  
  //Rcpp::XPtr<Pedigree> res(p, true);
  Rcpp::XPtr<Pedigree> res(p, false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = CharacterVector::create("malan_pedigree", "externalptr");
  
  return res;
}

//[[Rcpp::export]]
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid) {  
  Population* pop = population;
  
  Individual* ind = population->get_individual(pid);
  //Rcpp::XPtr<Individual> res(ind, true);
  Rcpp::XPtr<Individual> res(ind, false); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  res.attr("class") = CharacterVector::create("malan_individual", "externalptr");
  
  return res;
}

//' @export
// [[Rcpp::export]]
void print_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  int pid_f = (i->get_father() != nullptr) ? i->get_father()->get_pid() : -1;
  std::vector<Individual*>* children = i->get_children();
  
  Rcpp::Rcout << "  " << i->get_pid() << " with father " << pid_f << " and";
  
  if (children->size() == 0) {
    Rcpp::Rcout << " no children" << std::endl;
  } else {
    Rcpp::Rcout << " children (n = " << children->size() << "): " << std::endl;

    for (auto child : *children) {    
      std::vector<Individual*>* child_children = child->get_children();
      
      Rcpp::Rcout << "  " << child->get_pid() << " with father " << pid_f << " and " <<  child_children->size() << " children" << std::endl;
    }
  }
}


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


//[[Rcpp::export]]
void print_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::Rcout << "Pedigree with " << p->get_all_individuals()->size() << " individuals:" << std::endl;
  
  for (auto i : *inds) {    
    int pid_f = (i->get_father() != NULL) ? i->get_father()->get_pid() : -1;
    
    Rcpp::Rcout << "  " << i->get_pid() << " with father " << pid_f << std::endl;
  } 
}

//' get pids in pedigree
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  IntegerVector res(inds->size());
  int i = 0;
  for (auto ind : *inds) {   
    res(i) = ind->get_pid();
    ++i;
  } 
  
  return res;
}

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

//[[Rcpp::export]]
CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  CharacterMatrix edgelist(rels->size(), 2);
  int i = 0;
  
  for (auto pair: *rels) {
    edgelist(i, 0) = std::to_string(pair->first->get_pid());
    edgelist(i, 1) = std::to_string(pair->second->get_pid());
    ++i;
  }
  
  return edgelist;
}


//' Get pedigree information as graph (mainly intended for plotting)
//' 
//' @export
// [[Rcpp::export]]
List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  CharacterVector nodes(inds->size());
  
  int i = 0;
  for (auto individual : *inds) {
    nodes(i) = std::to_string(individual->get_pid());   
    ++i;
  }
  
  List ret;
  ret["nodes"] = nodes;
  ret["edgelist"] = get_pedigree_edgelist(ped);
  
  return ret;
}

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

//' Get pedigree from individual
//' 
//' @export
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;  
  Rcpp::XPtr<Pedigree> res(i->get_pedigree(), false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = CharacterVector::create("malan_pedigree", "externalptr");
  
  return res;
}

//' Get pedigree id from pid
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, IntegerVector pids) {  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  
  int N = pids.size();
  IntegerVector pedigree_ids(N);
  
  for (int i = 0; i < N; ++i) {
    Individual* ind = population->get_individual(pids[i]);
    pedigree_ids[i] = ind->get_pedigree_id();
  }
  
  return pedigree_ids;
}







//' @export
// [[Rcpp::export]]
void pedigree_populate_father_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, double mutation_rate) {  
  Pedigree* p = ped;
  ped->populate_father_haplotypes(loci, mutation_rate);
}

//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_father_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, double mutation_rate, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_father_haplotypes(loci, mutation_rate);
    
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


