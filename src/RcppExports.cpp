// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "malan_types.hpp"
#include <Rcpp.h>

using namespace Rcpp;

// wipe_pedigrees
void wipe_pedigrees(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP malan_wipe_pedigrees(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    wipe_pedigrees(pedigrees);
    return R_NilValue;
END_RCPP
}
// build_pedigrees
Rcpp::XPtr< std::vector<Pedigree*> > build_pedigrees(Rcpp::XPtr<Population> population, bool progress);
RcppExport SEXP malan_build_pedigrees(SEXP populationSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(build_pedigrees(population, progress));
    return rcpp_result_gen;
END_RCPP
}
// sample_geneology
List sample_geneology(size_t population_size, int generations, bool progress, int individuals_generations_return, bool verbose_result);
RcppExport SEXP malan_sample_geneology(SEXP population_sizeSEXP, SEXP generationsSEXP, SEXP progressSEXP, SEXP individuals_generations_returnSEXP, SEXP verbose_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type generations(generationsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< int >::type individuals_generations_return(individuals_generations_returnSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_result(verbose_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_geneology(population_size, generations, progress, individuals_generations_return, verbose_result));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_populate_father_haplotypes
void pedigree_populate_father_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, NumericVector mutation_rates);
RcppExport SEXP malan_pedigree_populate_father_haplotypes(SEXP pedSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mutation_rates(mutation_ratesSEXP);
    pedigree_populate_father_haplotypes(ped, loci, mutation_rates);
    return R_NilValue;
END_RCPP
}
// pedigrees_all_populate_father_haplotypes
void pedigrees_all_populate_father_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, NumericVector mutation_rates, bool progress);
RcppExport SEXP malan_pedigrees_all_populate_father_haplotypes(SEXP pedigreesSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mutation_rates(mutation_ratesSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    pedigrees_all_populate_father_haplotypes(pedigrees, loci, mutation_rates, progress);
    return R_NilValue;
END_RCPP
}
// pedigree_get_father_haplotype
std::vector<int> pedigree_get_father_haplotype(Rcpp::XPtr<Individual> individual);
RcppExport SEXP malan_pedigree_get_father_haplotype(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_get_father_haplotype(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_individual
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid);
RcppExport SEXP malan_get_individual(SEXP populationSEXP, SEXP pidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type pid(pidSEXP);
    rcpp_result_gen = Rcpp::wrap(get_individual(population, pid));
    return rcpp_result_gen;
END_RCPP
}
// get_pid
int get_pid(Rcpp::XPtr<Individual> individual);
RcppExport SEXP malan_get_pid(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pid(individual));
    return rcpp_result_gen;
END_RCPP
}
// print_individual
void print_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP malan_print_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    print_individual(individual);
    return R_NilValue;
END_RCPP
}
// get_pedigree_from_individual
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP malan_get_pedigree_from_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_from_individual(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_id_from_pid
IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, IntegerVector pids);
RcppExport SEXP malan_get_pedigree_id_from_pid(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pids(pidsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_id_from_pid(population, pids));
    return rcpp_result_gen;
END_RCPP
}
// malan_test
void malan_test();
RcppExport SEXP malan_malan_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    malan_test();
    return R_NilValue;
END_RCPP
}
// pop_size
int pop_size(Rcpp::XPtr<Population> population);
RcppExport SEXP malan_pop_size(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    rcpp_result_gen = Rcpp::wrap(pop_size(population));
    return rcpp_result_gen;
END_RCPP
}
// meioses_generation_distribution
IntegerMatrix meioses_generation_distribution(Rcpp::XPtr<Individual> individual, int generation_upper_bound_in_result);
RcppExport SEXP malan_meioses_generation_distribution(SEXP individualSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(meioses_generation_distribution(individual, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// population_size_generation
int population_size_generation(Rcpp::XPtr<Population> population, int generation_upper_bound_in_result);
RcppExport SEXP malan_population_size_generation(SEXP populationSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(population_size_generation(population, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_size_generation
int pedigree_size_generation(Rcpp::XPtr<Pedigree> pedigree, int generation_upper_bound_in_result);
RcppExport SEXP malan_pedigree_size_generation(SEXP pedigreeSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_size_generation(pedigree, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// pedigrees_count
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP malan_pedigrees_count(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigrees_count(pedigrees));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_size
int pedigree_size(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP malan_pedigree_size(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_size(ped));
    return rcpp_result_gen;
END_RCPP
}
// pedigrees_table
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP malan_pedigrees_table(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigrees_table(pedigrees));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index);
RcppExport SEXP malan_get_pedigree(SEXP pedigreesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree(pedigrees, index));
    return rcpp_result_gen;
END_RCPP
}
// print_pedigree
void print_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP malan_print_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    print_pedigree(ped);
    return R_NilValue;
END_RCPP
}
// get_pids_in_pedigree
IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP malan_get_pids_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pids_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_edgelist
CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP malan_get_pedigree_edgelist(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_edgelist(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_as_graph
List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP malan_get_pedigree_as_graph(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_as_graph(ped));
    return rcpp_result_gen;
END_RCPP
}
