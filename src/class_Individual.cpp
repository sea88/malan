#include "malan_types.hpp"

#include <stdexcept>

#include <Rcpp.h> // FIXME: Avoid Rcpp here? Only in api_* files?

/*
==========================================
Individual
==========================================
*/
Individual::Individual(int pid, int generation) {
  m_pid = pid;
  m_generation = generation;
  
  m_children = new std::vector<Individual*>();
}

Individual::~Individual() {
  delete m_children;
}

int Individual::get_pid() const {
  return m_pid;
}

int Individual::get_generation() const {
  return m_generation;
}

void Individual::add_child(Individual* child) {
  m_children->push_back(child);
}


void Individual::set_father(Individual* i) {
  // FIXME: Check sex of i?
  m_father = i;
}

Individual* Individual::get_father() const {
  return m_father;
}

std::vector<Individual*>* Individual::get_children() const {
  return m_children;
}

int Individual::get_children_count() const {
  return m_children->size();
}

bool Individual::pedigree_is_set() const {
  return (m_pedigree_id != 0);
}

int Individual::get_pedigree_id() const {
  return m_pedigree_id;
}

Pedigree* Individual::get_pedigree() const {
  return m_pedigree;
}

void Individual::set_pedigree_id(int id, Pedigree* ped, int* pedigree_size) {
  if (this->pedigree_is_set()) {
    return;
  }
  
  m_pedigree = ped;
  m_pedigree_id = id;
  *pedigree_size += 1;
  ped->add_member(this);
  
  if (m_father != nullptr) {  
    m_father->set_pedigree_id(id, ped, pedigree_size);
  }
  
  for (auto &child : (*m_children)) {
    ped->add_relation(this, child);
    child->set_pedigree_id(id, ped, pedigree_size);
  }
}

void Individual::set_alive_status(bool is_alive) {
  m_is_alive = is_alive;
}

bool Individual::get_alive_status() const {
  return m_is_alive;
}

void Individual::set_birth_year(int birth_year) {
  m_birth_year = birth_year;
}

int Individual::get_birth_year() const {
  return m_birth_year;
}

void Individual::dijkstra_reset() {
  m_dijkstra_visited = false;
  m_dijkstra_distance = 0;
}

void Individual::dijkstra_tick_distance(int step) {
  m_dijkstra_distance += step;
}

void Individual::dijkstra_set_distance_if_less(int dist) {
  if (m_dijkstra_distance < dist) {
    m_dijkstra_distance = dist;
  }
}

void Individual::dijkstra_mark_visited() {
  m_dijkstra_visited = true;
}

int Individual::dijkstra_get_distance() const {
  return m_dijkstra_distance; 
}

bool Individual::dijkstra_was_visited() const {
  return m_dijkstra_visited; 
}

// ASSUMES TREE!
//FIXME: Heavily relies on it being a tree, hence there is only one path connecting every pair of nodes
void Individual::meiosis_dist_tree_internal(Individual* dest, int* dist) const {
  if (this->get_pid() == dest->get_pid()) {
    //FIXME: Heavily relies on it being a tree, hence there is only one path connecting every pair of nodes
    *dist = dest->dijkstra_get_distance();
    return;
  }
  
  if (dest->dijkstra_was_visited()) {
    return;
  }
  
  dest->dijkstra_mark_visited();
  dest->dijkstra_tick_distance(1);
  int m = dest->dijkstra_get_distance();
  
  // FIXME: If not tree, then distance must be somehow checked if shorter and then adjusted
  
  Individual* father = dest->get_father();
  if (father != nullptr) {  
    //tree: ok
    father->dijkstra_tick_distance(m);

    // general? FIXME Correct?
    //father->dijkstra_set_distance_if_less(m);
    
    this->meiosis_dist_tree_internal(father, dist); 
  }
  
  std::vector<Individual*>* children = dest->get_children();
  for (auto child : *children) {
    //tree: ok
    child->dijkstra_tick_distance(m);

    // general? FIXME Correct?
    //child->dijkstra_set_distance_if_less(m);
    
    this->meiosis_dist_tree_internal(child, dist);
  }
}

// ASSUMES TREE!
int Individual::meiosis_dist_tree(Individual* dest) const {
  if (!(this->pedigree_is_set())) {
    throw std::invalid_argument("!(this->pedigree_is_set())");
  }
  
  if (dest == nullptr) {
    throw std::invalid_argument("dest is NULL");
  }
  
  if (!(dest->pedigree_is_set())) {
    throw std::invalid_argument("!(dest->pedigree_is_set())");
  }
  
  if (this->get_pedigree_id() != dest->get_pedigree_id()) {
    return -1;
  }
  
  std::vector<Individual*>* inds = this->get_pedigree()->get_all_individuals();
  for (auto child : *inds) {
    child->dijkstra_reset();
  }
  
  // At this point, the individuals this and dest belong to same pedigree
  int dist = 0;
  this->meiosis_dist_tree_internal(dest, &dist);
  return dist;
}



/*
Father haplotype
FIXME mutation_model?
*/
void Individual::father_haplotype_mutate(std::vector<double>& mutation_rates) {
  if (!m_father_haplotype_set) {
    Rcpp::stop("Father haplotype not set yet, so cannot mutate");
  }
  if (m_father_haplotype.size() != mutation_rates.size()) {
    Rcpp::stop("Number of loci specified in haplotype must equal number of mutation rates specified");
  }
  if (m_father_haplotype_mutated) {
    Rcpp::stop("Father haplotype already set and mutated");
  }
  
  
  for (int loc = 0; loc < m_father_haplotype.size(); ++loc) {
    if (R::runif(0.0, 1.0) < mutation_rates[loc]) {
      if (R::runif(0.0, 1.0) < 0.5) {
        m_father_haplotype[loc] = m_father_haplotype[loc] - 1;
      } else {
        m_father_haplotype[loc] = m_father_haplotype[loc] + 1;
      }
    }
  }
}

bool Individual::is_father_haplotype_set() const {
  return m_father_haplotype_set; 
}

void Individual::set_father_haplotype(std::vector<int> h) {
  m_father_haplotype = h;
  m_father_haplotype_set = true;
}

std::vector<int> Individual::get_father_haplotype() const {
  return m_father_haplotype;
}

void Individual::pass_haplotype_to_children(bool recursive, std::vector<double>& mutation_rates) {
  for (auto &child : (*m_children)) {
    child->set_father_haplotype(m_father_haplotype);
    child->father_haplotype_mutate(mutation_rates);
    
    if (recursive) {
      child->pass_haplotype_to_children(recursive, mutation_rates);
    }
  }
}

