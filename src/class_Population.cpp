#include <RcppArmadillo.h>
#include "malan_types.hpp"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/*
==========================================
Individual
==========================================
*/
Population::Population(std::unordered_map<int, Individual*>* population) {
  m_population = population;  
}

Population::~Population() {
  std::unordered_map<int, Individual*> pop = *m_population;
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    if (it->second == nullptr) {
      continue;
    }
    
    delete (it->second);
  }
  
  delete m_population;
}

std::unordered_map<int, Individual*>* Population::get_population() const {
  return m_population;
}

Individual* Population::get_individual(int pid) const {
  std::unordered_map<int, Individual*>::const_iterator got = m_population->find(pid);
  
  if (got == m_population->end()) {
    Rcpp::Rcerr << "Individual with pid = " << pid << " not found!" << std::endl;
    Rcpp::stop("Individual not found");
  }
  
  return got->second;
}

int Population::get_population_size() const {
  return m_population->size();
}
