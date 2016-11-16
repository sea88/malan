#include "malan_types.hpp"

#include <vector>
#include <unordered_map>

class Population {
private:
  std::unordered_map<int, Individual*>* m_population = NULL;

public:
  Population(std::unordered_map<int, Individual*>* population);
  ~Population();
  std::unordered_map<int, Individual*>* get_population() const;
  int get_population_size() const;
  Individual* get_individual(int pid) const;
};

