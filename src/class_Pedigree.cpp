#include "malan_types.hpp"

#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?

/*
==========================================
Pedigree
==========================================
*/

Pedigree::Pedigree(int id) {
  //Rcpp::Rcout << "Pedigree with id = " << id << " created" << std::endl;
  
  m_pedigree_id = id;
  
  m_all_individuals = new std::vector<Individual*>();
  m_relations = new std::vector< std::pair<Individual*, Individual*>* >();
}

Pedigree::~Pedigree() {
  delete m_all_individuals;
  
  for (auto it = m_relations->begin(); it != m_relations->end(); ++it) {
    delete *it;
  }
  
  delete m_relations;
}

int Pedigree::get_id() const {
  return m_pedigree_id;
}

void Pedigree::add_member(Individual* i) {
  m_all_individuals->push_back(i);
}

void Pedigree::add_relation(Individual* lhs, Individual* rhs) {
  //std::pair<Individual*, Individual*>* pair = new std::pair(lhs, rhs);
  std::pair<Individual*, Individual*>* pair = new std::pair<Individual*, Individual*>(lhs, rhs);
  m_relations->push_back(pair);
}

std::vector<Individual*>* Pedigree::get_all_individuals() const {
  return m_all_individuals;
}

std::vector< std::pair<Individual*, Individual*>* >* Pedigree::get_relations() const {
  return m_relations;
}



Individual* Pedigree::get_root() {
  if (m_root != nullptr) {
    return m_root;
  }
  
  /* FIXME: Exploits tree */
  bool root_set = false;
  
  for (auto &individual : (*m_all_individuals)) {
    if (individual->get_father() == nullptr) {
      if (root_set) {
        Rcpp::stop("Only expected one root in male pedigree!");
      } else {
        m_root = individual;
        root_set = true;
      }
      
      break;
    }
  }

  if (!root_set || m_root == nullptr) {
    Rcpp::stop("Expected a root in male pedigree!");
  }

  return m_root;
}


void Pedigree::populate_haplotypes(int loci, std::vector<double>& mutation_rates) {
  /* FIXME: Exploits tree */
  Individual* root = this->get_root();
  
  std::vector<int> h(loci); // initialises to 0, 0, ..., 0
  
  root->set_haplotype(h);
  root->pass_haplotype_to_children(true, mutation_rates);
}




