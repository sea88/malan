#ifndef MALAN_TYPES_H
#define MALAN_TYPES_H

#define CHECK_ABORT_EVERY 100000


//#define RCPP_XPTR_2ND_ARG true // ensures that finaliser is called
#define RCPP_XPTR_2ND_ARG false // do not call finalisers, I guess we live with some memory leaks for now...!!!

class Individual;
class Pedigree;
class Population;

#include "helper_Individual.hpp"

#include "class_Individual.hpp"
#include "class_Pedigree.hpp"
#include "class_Population.hpp"


#endif
