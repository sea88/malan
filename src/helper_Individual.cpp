#include <RcppArmadillo.h>
#include "malan_types.hpp"


bool find_path_from_root_to_dest(Individual* root, std::vector<Individual*>& path, const Individual* dest) {
  if (root == NULL || root == nullptr) {
    return false;
  }
  
  int dest_pid = dest->get_pid();

  path.push_back(root);

  if (root->get_pid() == dest_pid) {
    return true;
  }

  std::vector<Individual*>* children = root->get_children();
  for (auto child : *children) {
    if (find_path_from_root_to_dest(child, path, dest)) {
      return true;
    }
  }

  path.pop_back();
  return false;
}

