#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "malan_types.hpp"

using namespace Rcpp;


// based on sample_geneology
// @param generations -1 for simulate to 1 founder, else simulate this number of generations
//' @export
// [[Rcpp::export]]
List sample_geneology(size_t population_size, int generations, 
  double gamma_parameter_shape = 7, double gamma_parameter_scale = 7, 
  bool enable_gamma_variance_extension = false,
  bool progress = true, int individuals_generations_return = 2, bool verbose_result = false) {
  
  if (population_size <= 1) {
    Rcpp::stop("Please specify population_size > 1");
  }
  if (generations < -1 || generations == 0) {
    Rcpp::stop("Please specify generations as -1 (for simulation to 1 founder) or > 0");
  }

  if (enable_gamma_variance_extension) {
    if (gamma_parameter_shape <= 0.0) {
      Rcpp::stop("gamma_parameter_shape must be > 0.0");
    }
    if (gamma_parameter_scale <= 0.0) {
      Rcpp::stop("gamma_parameter_scale must be > 0.0");
    }
  }

  WFRandomFather wf_random_father(population_size);
  GammaVarianceRandomFather gamma_variance_father(population_size, gamma_parameter_shape, gamma_parameter_scale);  
  SimulateChooseFather* choose_father = &wf_random_father;
  
  if (enable_gamma_variance_extension) {
    choose_father = &gamma_variance_father;
  }
  
  bool simulate_fixed_number_generations = (generations == -1) ? false : true;
  
  //Rcpp::Rcout << simulate_fixed_number_generations << std::endl;
  
  Progress progress_bar((simulate_fixed_number_generations) ? generations : 1000, progress);
  
  IntegerMatrix individual_pids;
  IntegerMatrix father_pids;
  IntegerMatrix father_indices;
  
  std::vector<IntegerVector> individual_pids_tmp;
  std::vector<IntegerVector> father_pids_tmp;
  std::vector<IntegerVector> father_indices_tmp;
  IntegerVector individual_pids_tmp_vec;
  IntegerVector father_pids_tmp_vec;
  IntegerVector father_indices_tmp_vec;
  
  
  if (verbose_result) {
    if (simulate_fixed_number_generations) {
      individual_pids = IntegerMatrix(population_size, generations);
      father_pids = IntegerMatrix(population_size, generations);
      father_indices = IntegerMatrix(population_size, generations);
      
      std::fill(individual_pids.begin(), individual_pids.end(), NA_INTEGER);
      std::fill(father_pids.begin(), father_pids.end(), NA_INTEGER);
      std::fill(father_indices.begin(), father_indices.end(), NA_INTEGER);
    }
  }
  
  std::unordered_map<int, Individual*>* population_map = new std::unordered_map<int, Individual*>(); // pid's are garanteed to be unique
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  
  int individual_id = 1;
  std::vector<Individual*> end_generation(population_size);
  List end_generation_individuals(population_size);
  List last_k_generations_individuals;

  // Current generation: set-up
  if (verbose_result) {
    if (simulate_fixed_number_generations == false) {
      individual_pids_tmp_vec = IntegerVector(population_size);
      std::fill(individual_pids_tmp_vec.begin(), individual_pids_tmp_vec.end(), NA_INTEGER);
    }
  }
  
  for (size_t i = 0; i < population_size; ++i) {
    Individual* indv = new Individual(individual_id++, 0);
    end_generation[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    if (verbose_result) {
      if (simulate_fixed_number_generations) {
        individual_pids(i, 0) = indv->get_pid();
      } else {
        individual_pids_tmp_vec[i] = indv->get_pid();
      }
    }
    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_individuals[i] = indv_xptr;
    
    if (individuals_generations_return >= 0) {
      last_k_generations_individuals.push_back(indv_xptr);
    }
  }
  if (verbose_result) {
    if (simulate_fixed_number_generations == false) {
      individual_pids_tmp.push_back(individual_pids_tmp_vec);
    }
  }
  progress_bar.increment();
  
  // Next generation  
  //std::vector<Individual*>* children_generation = &end_generation;
  std::vector<Individual*> children_generation(population_size);
  for (size_t i = 0; i < population_size; ++i) children_generation[i] = end_generation[i];
  std::vector<Individual*> fathers_generation(population_size);
  
  int founders_left = population_size;
  
  // now, find out who the fathers to the children are
  //for (size_t generation = 1; generation < generations; ++generation) {
  size_t generation = 1;
  while ((simulate_fixed_number_generations == true && generation < generations) || (simulate_fixed_number_generations == false && founders_left > 1)) {
    int new_founders_left = 0;
    //Rcpp::Rcerr << "Generation " << generation << std::endl;
    
    // clear
    for (size_t i = 0; i < population_size; ++i) { // necessary?
      fathers_generation[i] = nullptr;
    }
    
    // for verbose result
    if (verbose_result) {
      individual_pids_tmp_vec = IntegerVector(population_size);
      father_pids_tmp_vec = IntegerVector(population_size);
      father_indices_tmp_vec = IntegerVector(population_size);
      
      std::fill(individual_pids_tmp_vec.begin(), individual_pids_tmp_vec.end(), NA_INTEGER);
      std::fill(father_pids_tmp_vec.begin(), father_pids_tmp_vec.end(), NA_INTEGER);
      std::fill(father_indices_tmp_vec.begin(), father_indices_tmp_vec.end(), NA_INTEGER);
    }
    
    choose_father->update_state_new_generation();
    
    // now, run through children to pick each child's father
    for (size_t i = 0; i < population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      // child [i] in [generation-1]/children_generation has father [father_i] in [generation]/fathers_generation
      //int father_i = sample_person_weighted(population_size, fathers_prob, fathers_prob_perm);
      int father_i = choose_father->get_father_i();
      
      // if this is the father's first child, create the father
      if (fathers_generation[father_i] == nullptr) {
        Individual* father = new Individual(individual_id++, generation);
        fathers_generation[father_i] = father;
        (*population_map)[father->get_pid()] = father;      
        
        if (verbose_result) {
          individual_pids_tmp_vec[father_i] = father->get_pid();
        }
        
        new_founders_left++;

        if (generation <= individuals_generations_return) {
          Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG);
          last_k_generations_individuals.push_back(father_xptr);
        }
      }
      
      if (verbose_result) {
        father_pids_tmp_vec[i] = fathers_generation[father_i]->get_pid();
        father_indices_tmp_vec[i] = father_i + 1; // 1 to get R's 1-indexed
      }      
      
      children_generation[i]->set_father(fathers_generation[father_i]);
      fathers_generation[father_i]->add_child(children_generation[i]);
    }
    
    // verbose result
    if (verbose_result) {
      if (simulate_fixed_number_generations) {
        individual_pids(Rcpp::_, generation) = individual_pids_tmp_vec;
        father_pids(Rcpp::_, generation - 1) = father_pids_tmp_vec;
        father_indices(Rcpp::_, generation - 1) = father_indices_tmp_vec;      
      } else {
        individual_pids_tmp.push_back(individual_pids_tmp_vec);
        father_pids_tmp.push_back(father_pids_tmp_vec);
        father_indices_tmp.push_back(father_indices_tmp_vec);  
      }
    }
        
    // children_generation = &fathers_generation;
    for (size_t i = 0; i < population_size; ++i) children_generation[i] = fathers_generation[i];
    
    if (Progress::check_abort()) {
      stop("Aborted");
    }
    
    if (progress) {
      progress_bar.increment();
    }
    
    founders_left = new_founders_left;
    generation += 1;
  }
  
  if (verbose_result) {
    if (simulate_fixed_number_generations == false) {
      // Fill in last NA column
      father_pids_tmp_vec = IntegerVector(population_size);
      father_indices_tmp_vec = IntegerVector(population_size);
      std::fill(father_pids_tmp_vec.begin(), father_pids_tmp_vec.end(), NA_INTEGER);
      std::fill(father_indices_tmp_vec.begin(), father_indices_tmp_vec.end(), NA_INTEGER);
      father_pids_tmp.push_back(father_pids_tmp_vec);
      father_indices_tmp.push_back(father_indices_tmp_vec); 


      int generations_final = generation;
      
      individual_pids = IntegerMatrix(population_size, generations_final);
      father_pids = IntegerMatrix(population_size, generations_final);
      father_indices = IntegerMatrix(population_size, generations_final);
      
      std::fill(individual_pids.begin(), individual_pids.end(), NA_INTEGER);
      std::fill(father_pids.begin(), father_pids.end(), NA_INTEGER);
      std::fill(father_indices.begin(), father_indices.end(), NA_INTEGER);
      
      for (int g = 0; g < individual_pids_tmp.size(); ++g) {
        for (size_t i = 0; i < population_size; ++i) {
          individual_pids(i, g) = individual_pids_tmp[g][i];
          father_pids(i, g) = father_pids_tmp[g][i];
          father_indices(i, g) = father_indices_tmp[g][i];
        }      
      }      
    }
  }
  
  List res;
  res["population"] = population_xptr;
  res["generations"] = generation;
  res["founders"] = founders_left;
  res["type"] = (enable_gamma_variance_extension) ? "GammaVariation" : "SimpleWF";  
  res["end_generation_individuals"] = end_generation_individuals;
  res["individuals_generations"] = last_k_generations_individuals;

  if (verbose_result) {
    res["individual_pids"] = individual_pids;
    res["father_pids"] = father_pids;
    res["father_indices"] = father_indices;
  }
  
  return res;
}

