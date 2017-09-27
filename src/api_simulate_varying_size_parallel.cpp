// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#if RCPP_PARALLEL_USE_TBB
#include "tbb/concurrent_vector.h"
#endif

#include <progress.hpp>

#include "malan_types.hpp"
#include "api_simulate.hpp"

#ifdef MALAN_PROFILE
#include <gperftools/profiler.h>
#endif

using namespace Rcpp;

#ifndef RCPP_PARALLEL_USE_TBB

List sample_geneology_varying_size_parallel(
  IntegerVector population_sizes,
  int extra_generations_full = 0,  
  double gamma_parameter_shape = 7, double gamma_parameter_scale = 7, 
  bool enable_gamma_variance_extension = false,
  bool progress = true, 
  int individuals_generations_return = 2,
  int threads = 1) {
  
  List res;
  stop("Intel TBB support required, but it was not found. Please use non-parallel version.");
  return res;
}

#else

struct InitialisePopulation : public RcppParallel::Worker {
  tbb::concurrent_vector<Individual*>* m_init_indvs;
  tbb::concurrent_vector<Individual*>* m_children_gen;
  tbb::concurrent_vector<Individual*>* m_end_generation;
  std::unordered_map<int, Individual*>* m_population_map;
  tbb::mutex* m_population_map_mutex;

  InitialisePopulation(tbb::concurrent_vector<Individual*>* init_indvs, 
                       tbb::concurrent_vector<Individual*>* children_gen, 
                       tbb::concurrent_vector<Individual*>* end_generation,
                       std::unordered_map<int, Individual*>* population_map,
                       tbb::mutex* population_map_mutex) {
    m_init_indvs = init_indvs;
    m_children_gen = children_gen;
    m_end_generation = end_generation;
    m_population_map = population_map;
    m_population_map_mutex = population_map_mutex;
  }

  void operator()(std::size_t begin, std::size_t end) {
    //Rcpp::Rcout << "i = " << begin << "; i < " << end << std::endl;
    
    for (size_t i = begin; i < end; ++i) {
      const int individual_id = i + 1;
      
      Individual* indv = new Individual(individual_id, 0);
      m_init_indvs->at(i) = indv;
      m_children_gen->at(i) = indv;
      m_end_generation->at(i) = indv;
      
      m_population_map_mutex->lock();
      m_population_map->insert({ individual_id, indv });
      m_population_map_mutex->unlock();
    }      
  }
};

/*
struct EvolveGeneration : public RcppParallel::Worker {
  int m_generation;
  tbb::concurrent_vector<Individual*>* m_children_gen;
  tbb::concurrent_vector<int>* m_child_i_father_i_map;
  tbb::concurrent_vector<Individual*>* m_fathers_generation;
  std::unordered_map<int, Individual*>* m_population_map;
  tbb::mutex* m_population_map_mutex;
  int* m_individual_id;
  int* m_new_founders_left;

  EvolveGeneration(int generation,
                   tbb::concurrent_vector<Individual*>* children_gen,
                   tbb::concurrent_vector<int>* child_i_father_i_map,
                   tbb::concurrent_vector<Individual*>* fathers_generation,
                   std::unordered_map<int, Individual*>* population_map,
                   tbb::mutex* population_map_mutex,
                   int* individual_id,
                   int* new_founders_left) {
    m_generation = generation;
    m_children_gen = children_gen;
    m_child_i_father_i_map = child_i_father_i_map;
    m_fathers_generation = fathers_generation;
    m_population_map = population_map;
    m_population_map_mutex = population_map_mutex;
    m_individual_id = individual_id;
    m_new_founders_left = new_founders_left;
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (m_children_gen->at(i) == nullptr) {
        continue;
      }
      
      // FIXME: Consider locking...
      m_population_map_mutex->lock();
      
      // child [i] in [generation-1]/m_children_gen has father [father_i] in [generation]/fathers_generation
      int father_i = m_child_i_father_i_map->at(i);
      
      // if this is the father's first child, create the father
      Individual* father = m_fathers_generation->at(father_i);
      
      if (father == nullptr) {
        // m_fathers_generation[father_i] is now made an Individual* and not nullptr any more
        
        //m_population_map_mutex->lock();
        int ind_id = *m_individual_id;
        father = new Individual(ind_id, m_generation);
        m_population_map->insert({ ind_id, father });
        *m_individual_id = ind_id + 1;
        *m_new_founders_left = (*m_new_founders_left) + 1;
        //m_population_map_mutex->unlock();
        
        m_fathers_generation->at(father_i) = father;
      }
      
      m_population_map_mutex->unlock();
      
      Individual* child = m_children_gen->at(i);
      child->set_father(father);
      father->add_child(child);
    }      
  }
};
*/

struct GenerateFathers : public RcppParallel::Worker {
  int m_generation;
  tbb::concurrent_vector<Individual*>* m_fathers_generation;
  std::vector<int>* m_fathers_to_create;
  std::vector<int>* m_fathers_pids;
  std::unordered_map<int, Individual*>* m_population_map;
  tbb::mutex* m_population_map_mutex;

  GenerateFathers(int generation, 
                  tbb::concurrent_vector<Individual*>* fathers_generation,
                  std::vector<int>* fathers_to_create,
                  std::vector<int>* fathers_pids,
                  std::unordered_map<int, Individual*>* population_map,
                  tbb::mutex* population_map_mutex) {

    m_generation = generation;
    m_fathers_generation = fathers_generation;
    m_fathers_to_create = fathers_to_create;
    m_fathers_pids = fathers_pids;
    m_population_map = population_map;
    m_population_map_mutex = population_map_mutex;
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (m_fathers_to_create->at(i) == 0) {
        continue;
      }
      
      int father_pid = m_fathers_pids->at(i);
      
      Individual* father = new Individual(father_pid, m_generation);
      m_fathers_generation->at(i) = father;
      
      m_population_map_mutex->lock();
      m_population_map->insert({ father_pid, father });
      m_population_map_mutex->unlock();
    }
  }
};


//' Simulate a geneology with varying population size.
//' 
//' This function simulates a geneology with varying population size specified
//' by a vector of population sizes, one for each generation. 
//' 
//' By the backwards simulating process of the Wright-Fisher model, 
//' individuals with no descendants in the end population are not simulated 
//' If for some reason additional full generations should be simulated, 
//' the number can be specified via the \code{extra_generations_full} parameter.
//' This can for example be useful if one wants to simulate the 
//' final 3 generations although some of these may not get (male) children.
//' 
//' Let \eqn{\alpha} be the parameter of a symmetric Dirichlet distribution 
//' specifying each man's probability to be the father of an arbitrary 
//' male in the next generation. When \eqn{\alpha = 5}, a man's relative probability 
//' to be the father has 95\% probability to lie between 0.32 and 2.05, compared with a 
//' constant 1 under the standard Wright-Fisher model and the standard deviation in 
//' the number of male offspring per man is 1.10 (standard Wright-Fisher = 1).
//' 
//' This symmetric Dirichlet distribution is implemented by drawing 
//' father (unscaled) probabilities from a Gamma distribution with 
//' parameters \code{gamma_parameter_shape} and \code{gamma_parameter_scale} 
//' that are then normalised to sum to 1. 
//' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
//' the following must be used:
//' \eqn{\code{gamma_parameter_shape} = \alpha}
//' and 
//' \eqn{\code{gamma_parameter_scale} = 1/\alpha}.
//' 
//' @param population_sizes The size of the population at each generation, g. 
//'        population_sizes[g] is the population size at generation g.
//'        The length of population_sizes is the number of generations being simulated.
//' @param extra_generations_full Additional full generations to be simulated.
//' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
//' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
//' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
//' @param progress Show progress.
//' @param individuals_generations_return How many generations back to return (pointers to) individuals for.
//' 
//' @return A malan_simulation / list with the following entries:
//' \itemize{
//'   \item \code{population}. An external pointer to the population.
//'   \item \code{generations}. Generations actually simulated, mostly useful when parameter \code{generations = -1}.
//'   \item \code{founders}. Number of founders after the simulated \code{generations}.
//'   \item \code{growth_type}. Growth type model.
//'   \item \code{sdo_type}. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on \code{enable_gamma_variance_extension}.
//'   \item \code{end_generation_individuals}. Pointers to individuals in end generation.
//'   \item \code{individuals_generations}. Pointers to individuals in end generation in addition to the previous \code{individuals_generations_return}.
//' }
//' 
//' @import Rcpp
//' @import RcppProgress
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
List sample_geneology_varying_size_parallel(
  IntegerVector population_sizes,
  int extra_generations_full = 0,  
  double gamma_parameter_shape = 7, double gamma_parameter_scale = 7, 
  bool enable_gamma_variance_extension = false,
  bool progress = true, 
  int individuals_generations_return = 2,
  int threads = 0,
  int do_parallel_when_individuals_above = 0) {
  
  #ifdef MALAN_PROFILE
  ProfilerStart("/home/mikl/work-aau/research-projects/2016-Balding-shared-haplotypes/model-31-malan-speed-up/prof.log");
  #endif

  if (threads < 0) {
    Rcpp::stop("threads must be >= 0 (0 for automatic)");
  }
  
  if (threads > 0) {
    static tbb::task_scheduler_init* s_pTaskScheduler = NULL;
    try
    {
      if (!s_pTaskScheduler) {
        s_pTaskScheduler = new tbb::task_scheduler_init(threads, 0);
      } else {
        s_pTaskScheduler->terminate();
        s_pTaskScheduler->initialize(threads, 0); 
      }
    } catch (...) {
      stop("Error loading TBB: (Unknown error)");
    }
  }

  // boolean chosen like this to obey NA's
  bool all_gt_1 = is_true(all(population_sizes >= 1));
  if (!all_gt_1) {
    Rcpp::stop("Please specify only population_sizes >= 1");
  }
  
  int generations = population_sizes.length();
  
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
  
  //Rcpp::Rcout << "vary 1" << std::endl;

  
  Progress progress_bar(generations, progress);
  
  // pid's are garanteed to be unique
  std::unordered_map<int, Individual*>* population_map = new std::unordered_map<int, Individual*>();
  Population* population = new Population(population_map);
  tbb::mutex population_map_mutex;  
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  
  //std::vector<Individual*> end_generation(population_sizes[generations-1]);
  tbb::concurrent_vector<Individual*> end_generation(population_sizes[generations-1]);
  List end_generation_individuals(population_sizes[generations-1]);
  List last_k_generations_individuals;
  
  // FIXME: Make last_k_generations_individuals big enough already here? The size
  //        can easily be calculated.
  
  if (individuals_generations_return >= 0) {
    last_k_generations_individuals = List(population_sizes[generations-1]);
  }  
  
  int initial_popsize = population_sizes[generations-1];
  tbb::concurrent_vector<Individual*> initial_individuals(initial_popsize);
  tbb::concurrent_vector<Individual*> children_generation(initial_popsize);
  
  InitialisePopulation init_pop(&initial_individuals, 
                                &children_generation, 
                                &end_generation, 
                                population_map, 
                                &population_map_mutex);
                                
  if (initial_popsize >= do_parallel_when_individuals_above) {
    parallelFor(0, initial_popsize, init_pop); // does i = 0; i < initial_popsize
  } else {
    init_pop(0, initial_popsize);
  }
  

  // Cannot be parallised as R memory is called (Rcpp::XPtr, List, ...) 
  // and that must only be done in the main thread.
  for (size_t i = 0; i < population_sizes[generations-1]; ++i) {
    const int individual_id = i + 1;
    Individual* indv = initial_individuals[i];

    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_individuals[i] = indv_xptr;
    
    if (individuals_generations_return >= 0) {
      last_k_generations_individuals[i] = indv_xptr;
    }
  }
  
  // Be ready with id's for next generations:
  // In the initialisation, we created population_sizes[generations-1] individuals
  int individual_id = population_sizes[generations - 1] + 1; // +1: next individual must be the next
  
  if (progress) {
    progress_bar.increment();
  }
  
  // Next generation ready from children_generation prepared by InitialisePopulation
  

  //std::vector<Individual*> fathers_generation;
  tbb::concurrent_vector<Individual*> fathers_generation;
  
  int founders_left = population_sizes[generations-1];
  
  // Now, find out who the fathers to the children are
  for (size_t generation = 1; generation < generations; ++generation) {
    int population_size = population_sizes[generations-(generation+1)];    
    int children_population_size = population_sizes[generations-generation];
    
    // This is not possible to parallelise as it uses 
    // R/Rcpp API: R::runif / Rcpp::rgamma
    WFRandomFather wf_random_father(population_size);
    GammaVarianceRandomFather gamma_variance_father(population_size, gamma_parameter_shape, gamma_parameter_scale);  
    SimulateChooseFather* choose_father = &wf_random_father;
    if (enable_gamma_variance_extension) {
      choose_father = &gamma_variance_father;
    }
    
    //Rcpp::Rcout << "vary 5-" << generation << std::endl;
    
    fathers_generation.clear();
    fathers_generation.resize(population_size);
    
    // <- Init
    
    int new_founders_left = 0;
    //Rcpp::Rcerr << "Generation " << generation << std::endl;
    
    // clear
    for (size_t i = 0; i < population_size; ++i) { // necessary?
      fathers_generation[i] = nullptr;
    }
    
    choose_father->update_state_new_generation();
    
    /*
    FIXME: New, more explicit logic? Maybe copy create_father_update_simulation_state_varying_size into here?
    */
    
    /*
    // now, run through children to pick each child's father
    // FIXME: PARALLELISE?? BLOCK ON
    //     int* individual_id, 
    //     std::vector<Individual*>& fathers_generation, 
    //     std::unordered_map<int, Individual*>* population_map, 
    //     int* new_founders_left,
    //     List& last_k_generations_individuals
    for (size_t i = 0; i < children_population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      // child [i] in [generation-1]/children_generation has father [father_i] in [generation]/fathers_generation
      //int father_i = sample_person_weighted(population_size, fathers_prob, fathers_prob_perm);
      int father_i = choose_father->get_father_i();
      
      // if this is the father's first child, create the father
      if (fathers_generation[father_i] == nullptr) {        
        // fathers_generation[father_i] is now an Individual* and not nullptr any more
        create_father_update_simulation_state_varying_size_parallel(father_i, &individual_id, generation, 
              individuals_generations_return, fathers_generation, population_map, 
              &new_founders_left, last_k_generations_individuals);
      }

      children_generation[i]->set_father(fathers_generation[father_i]);
      fathers_generation[father_i]->add_child(children_generation[i]);
    }
    */
    
    // cannot be parallised due to choose_father
    // is only going to be read afterwards, so stl vector would suffice
    std::vector<int> child_i_father_i_map(children_population_size);
    std::vector<int> fathers_to_create(population_size);
    std::vector<int> fathers_pids(population_size);
    
    for (size_t i = 0; i < children_population_size; ++i) {
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      int father_id = choose_father->get_father_i();
      child_i_father_i_map[i] = father_id;
      
      if (fathers_to_create[father_id] == 0) {
        fathers_pids[father_id] = individual_id;
        individual_id += 1;
        new_founders_left += 1;
      }
      
      fathers_to_create[father_id] += 1;
    }
    
    GenerateFathers fathers_gen(generation,                                
                                &fathers_generation,
                                &fathers_to_create,
                                &fathers_pids, 
                                population_map, 
                                &population_map_mutex);
    
    // Only in parallel if it makes sense
    if (new_founders_left >= do_parallel_when_individuals_above) {
      parallelFor(0, population_size, fathers_gen); // does i = 0; i < children_population_size
    } else {
      fathers_gen(0, population_size);
    }
    
    
    for (size_t i = 0; i < children_population_size; ++i) {
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      int father_id = child_i_father_i_map[i];
      Individual* child = children_generation[i];
      Individual* father = fathers_generation[father_id];
      
      child->set_father(father);
      father->add_child(child);
      
      if (generation <= individuals_generations_return) {
        Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG);
        last_k_generations_individuals.push_back(father_xptr);
      }
    }
    
    
    
    //Rcpp::Rcout << "vary 8-" << generation << std::endl;
    
    //Rcpp::Rcout << "generation = " << generation << "; extra_generations_full = " << extra_generations_full << std::endl;
    
    // create additional fathers (without children) if needed:
    if (generation <= extra_generations_full) {
      for (size_t father_i = 0; father_i < population_size; ++father_i) {
        //Rcpp::Rcout << "vary 8-dong-" << father_i << std::endl;
        
        if (fathers_generation[father_i] != nullptr) {
          continue;
        }        
        
        // create father, no children etc.
        // create_father_update_simulation_state_varying_size ->
          Individual* father = new Individual(individual_id, generation);
          fathers_generation[father_i] = father;
          population_map->insert({ individual_id, father });
          
          individual_id = individual_id + 1;
          new_founders_left = new_founders_left + 1;

          if (generation <= individuals_generations_return) {
            //Rcpp::Rcout << "create_father_update_simulation_state_varying_size_parallel: generation = " << generation << "; individuals_generations_return = " << individuals_generations_return << std::endl;
            
            Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG);
            last_k_generations_individuals.push_back(father_xptr);
          }  
        // <- create_father_update_simulation_state_varying_size
      }      
    }
    
    //Rcpp::Rcout << "vary 9-" << generation << std::endl;
    
    // children_generation = &fathers_generation;
    // FIXME mikl 2017-06-26 19:09
    //children_generation = fathers_generation;
    //vs:
    children_generation.clear();
    children_generation.resize(population_size);
    for (size_t i = 0; i < population_size; ++i) {
      children_generation[i] = fathers_generation[i];
    }
    //<-
    
    if (Progress::check_abort()) {
      stop("Aborted");
    }
    
    if (progress) {
      progress_bar.increment();
    }
    
    founders_left = new_founders_left;
  }
  
  List res;
  res["population"] = population_xptr;
  res["generations"] = generations;
  res["founders"] = founders_left;
  res["growth_type"] = "VaryingPopulationSize";
  res["sdo_type"] = (enable_gamma_variance_extension) ? "GammaVariation" : "StandardWF";
  res["end_generation_individuals"] = end_generation_individuals;
  res["individuals_generations"] = last_k_generations_individuals;
  res["threads"] = threads;

  res.attr("class") = CharacterVector::create("malan_simulation", "list");
  
  #ifdef MALAN_PROFILE
  ProfilerStop();
  #endif
  
  return res;
}

#endif // #ifdef RCPP_PARALLEL_USE_TBB


