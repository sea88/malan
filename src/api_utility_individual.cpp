#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.hpp"


//[[Rcpp::export]]
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid) {  
  Population* pop = population;
  
  Individual* ind = population->get_individual(pid);
  //Rcpp::XPtr<Individual> res(ind, true);
  Rcpp::XPtr<Individual> res(ind, false); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  res.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");
  
  return res;
}


//' @export
// [[Rcpp::export]]
int get_pid(Rcpp::XPtr<Individual> individual) {  
  return individual->get_pid();
}

//' @export
// [[Rcpp::export]]
void print_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  int pid_f = (i->get_father() != nullptr) ? i->get_father()->get_pid() : -1;
  std::vector<Individual*>* children = i->get_children();
  
  Rcpp::Rcout << "  pid = " << i->get_pid() << " with father pid = " << pid_f << " and";
  
  if (children->size() == 0) {
    Rcpp::Rcout << " no children" << std::endl;
  } else {
    Rcpp::Rcout << " children (n = " << children->size() << "): " << std::endl;

    for (auto child : *children) {    
      std::vector<Individual*>* child_children = child->get_children();
      
      Rcpp::Rcout << "    pid = " << child->get_pid() << " with father pid = " << pid_f << " and " <<  child_children->size() << " children" << std::endl;
    }
  }
}

//' Get pedigree from individual
//' 
//' @export
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;  
  Rcpp::XPtr<Pedigree> res(i->get_pedigree(), false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = Rcpp::CharacterVector::create("malan_pedigree", "externalptr");
  
  return res;
}

//' Get pedigree id from pid
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids) {  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  
  int N = pids.size();
  Rcpp::IntegerVector pedigree_ids(N);
  
  for (int i = 0; i < N; ++i) {
    Individual* ind = population->get_individual(pids[i]);
    pedigree_ids[i] = ind->get_pedigree_id();
  }
  
  return pedigree_ids;
}


