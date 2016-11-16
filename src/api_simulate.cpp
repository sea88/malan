#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "malan_types.hpp"

using namespace Rcpp;



// Number in {0, 1, ..., population_size - 1}: perfect for 0-indexed
int sample_person(size_t population_size) {
  return R::runif(0, 1)*((double)population_size);
}


//' @export
// [[Rcpp::export]]
List sample_geneology(size_t population_size, size_t generations, bool progress = true, bool verbose_result = false) {
  if (population_size <= 1) {
    Rcpp::stop("Please specify population_size > 1");
  }
  if (generations <= 1) {
    Rcpp::stop("Please specify generations > 1");
  }
  
  Progress progress_bar(generations, progress);
  
  
  
  IntegerMatrix individual_pids;
  IntegerMatrix father_pids;
  IntegerMatrix father_indices;
  
  
  if (verbose_result) {
    individual_pids = IntegerMatrix(population_size, generations);
    father_pids = IntegerMatrix(population_size, generations);
    father_indices = IntegerMatrix(population_size, generations);
    
    std::fill(individual_pids.begin(), individual_pids.end(), NA_INTEGER);
    std::fill(father_pids.begin(), father_pids.end(), NA_INTEGER);
    std::fill(father_indices.begin(), father_indices.end(), NA_INTEGER);
  
  }
  
  
  std::unordered_map<int, Individual*>* population_map = new std::unordered_map<int, Individual*>(); // pid's are garanteed to be unique
  Population* population = new Population(population_map);
  
  int individual_id = 1;
  std::vector<Individual*> end_generation(population_size);

  // Current generation: set-up
  for (size_t i = 0; i < population_size; ++i) {
    Individual* indv = new Individual(individual_id++);
    end_generation[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    individual_pids(i, 0) = indv->get_pid();
  }
  progress_bar.increment();
  
  // Next generation  
  //std::vector<Individual*>* children_generation = &end_generation;
  std::vector<Individual*> children_generation(population_size);
  for (size_t i = 0; i < population_size; ++i) children_generation[i] = end_generation[i];
  std::vector<Individual*> fathers_generation(population_size);
  
  // now, find out who the fathers to the children are
  for (size_t generation = 1; generation < generations; ++generation) {
    // clear
    for (size_t i = 0; i < population_size; ++i) { // necessary?
      fathers_generation[i] = nullptr;
    }
    
    // now, run through children to pick each child's father
    for (size_t i = 0; i < population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      // child [i] in [generation-1]/children_generation has father [father_i] in [generation]/fathers_generation
      int father_i = sample_person(population_size);
      
      // if this is the father's first child, create the father
      if (fathers_generation[father_i] == nullptr) {
        Individual* father = new Individual(individual_id++);
        fathers_generation[father_i] = father;
        (*population_map)[father->get_pid()] = father;      
        
        //individual_pids(i, generation) = father->get_pid();
        individual_pids(father_i, generation) = father->get_pid();
      }
      
      father_pids(i, generation-1) = fathers_generation[father_i]->get_pid();
      father_indices(i, generation-1) = father_i + 1; // 1 to get R's 1-indexed
      
      children_generation[i]->set_father(fathers_generation[father_i]);
      fathers_generation[father_i]->add_child(children_generation[i]);
    }
    
    //children_generation = &fathers_generation;
    for (size_t i = 0; i < population_size; ++i) children_generation[i] = fathers_generation[i];
    
    if (Progress::check_abort()) {
      stop("Aborted");
    }
    
    if (progress) {
      progress_bar.increment();
    }
  }
  
  Rcpp::XPtr<Population> population_xptr(population, true);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  List res;
  res["population"] = population_xptr;
  
  if (verbose_result) {
    res["individual_pids"] = individual_pids;
    res["father_pids"] = father_pids;
    res["father_indices"] = father_indices;
  }
  
  return res;
}



