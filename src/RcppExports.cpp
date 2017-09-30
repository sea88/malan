// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "malan_types.hpp"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// wipe_pedigrees
void wipe_pedigrees(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _malan_wipe_pedigrees(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    wipe_pedigrees(pedigrees);
    return R_NilValue;
END_RCPP
}
// build_pedigrees
Rcpp::XPtr< std::vector<Pedigree*> > build_pedigrees(Rcpp::XPtr<Population> population, bool progress);
RcppExport SEXP _malan_build_pedigrees(SEXP populationSEXP, SEXP progressSEXP) {
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
List sample_geneology(size_t population_size, int generations, int extra_generations_full, double gamma_parameter_shape, double gamma_parameter_scale, bool enable_gamma_variance_extension, bool progress, int individuals_generations_return, bool verbose_result);
RcppExport SEXP _malan_sample_geneology(SEXP population_sizeSEXP, SEXP generationsSEXP, SEXP extra_generations_fullSEXP, SEXP gamma_parameter_shapeSEXP, SEXP gamma_parameter_scaleSEXP, SEXP enable_gamma_variance_extensionSEXP, SEXP progressSEXP, SEXP individuals_generations_returnSEXP, SEXP verbose_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type generations(generationsSEXP);
    Rcpp::traits::input_parameter< int >::type extra_generations_full(extra_generations_fullSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_shape(gamma_parameter_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_scale(gamma_parameter_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type enable_gamma_variance_extension(enable_gamma_variance_extensionSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< int >::type individuals_generations_return(individuals_generations_returnSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_result(verbose_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_geneology(population_size, generations, extra_generations_full, gamma_parameter_shape, gamma_parameter_scale, enable_gamma_variance_extension, progress, individuals_generations_return, verbose_result));
    return rcpp_result_gen;
END_RCPP
}
// sample_geneology_varying_size
List sample_geneology_varying_size(IntegerVector population_sizes, int extra_generations_full, double gamma_parameter_shape, double gamma_parameter_scale, bool enable_gamma_variance_extension, bool progress, int individuals_generations_return);
RcppExport SEXP _malan_sample_geneology_varying_size(SEXP population_sizesSEXP, SEXP extra_generations_fullSEXP, SEXP gamma_parameter_shapeSEXP, SEXP gamma_parameter_scaleSEXP, SEXP enable_gamma_variance_extensionSEXP, SEXP progressSEXP, SEXP individuals_generations_returnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type population_sizes(population_sizesSEXP);
    Rcpp::traits::input_parameter< int >::type extra_generations_full(extra_generations_fullSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_shape(gamma_parameter_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_parameter_scale(gamma_parameter_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type enable_gamma_variance_extension(enable_gamma_variance_extensionSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< int >::type individuals_generations_return(individuals_generations_returnSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_geneology_varying_size(population_sizes, extra_generations_full, gamma_parameter_shape, gamma_parameter_scale, enable_gamma_variance_extension, progress, individuals_generations_return));
    return rcpp_result_gen;
END_RCPP
}
// indices_in_mixture
Rcpp::List indices_in_mixture(Rcpp::IntegerMatrix haplotypes, Rcpp::IntegerVector H1, Rcpp::IntegerVector H2);
RcppExport SEXP _malan_indices_in_mixture(SEXP haplotypesSEXP, SEXP H1SEXP, SEXP H2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type haplotypes(haplotypesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type H1(H1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type H2(H2SEXP);
    rcpp_result_gen = Rcpp::wrap(indices_in_mixture(haplotypes, H1, H2));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_get_haplotypes_pids
Rcpp::List pedigree_get_haplotypes_pids(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids);
RcppExport SEXP _malan_pedigree_get_haplotypes_pids(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pids(pidsSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_get_haplotypes_pids(population, pids));
    return rcpp_result_gen;
END_RCPP
}
// individuals_get_haplotypes
Rcpp::IntegerMatrix individuals_get_haplotypes(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals);
RcppExport SEXP _malan_individuals_get_haplotypes(SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf< Rcpp::XPtr<Individual> > >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(individuals_get_haplotypes(individuals));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_populate_haplotypes
void pedigree_populate_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, Rcpp::NumericVector mutation_rates);
RcppExport SEXP _malan_pedigree_populate_haplotypes(SEXP pedSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    pedigree_populate_haplotypes(ped, loci, mutation_rates);
    return R_NilValue;
END_RCPP
}
// pedigrees_all_populate_haplotypes
void pedigrees_all_populate_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, Rcpp::NumericVector mutation_rates, bool progress);
RcppExport SEXP _malan_pedigrees_all_populate_haplotypes(SEXP pedigreesSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    pedigrees_all_populate_haplotypes(pedigrees, loci, mutation_rates, progress);
    return R_NilValue;
END_RCPP
}
// pedigrees_all_populate_haplotypes_ladder_bounded
void pedigrees_all_populate_haplotypes_ladder_bounded(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, Rcpp::NumericVector mutation_rates, Rcpp::IntegerVector ladder_max_dist_0, bool progress);
RcppExport SEXP _malan_pedigrees_all_populate_haplotypes_ladder_bounded(SEXP pedigreesSEXP, SEXP lociSEXP, SEXP mutation_ratesSEXP, SEXP ladder_max_dist_0SEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mutation_rates(mutation_ratesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ladder_max_dist_0(ladder_max_dist_0SEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    pedigrees_all_populate_haplotypes_ladder_bounded(pedigrees, loci, mutation_rates, ladder_max_dist_0, progress);
    return R_NilValue;
END_RCPP
}
// get_haplotype
std::vector<int> get_haplotype(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_get_haplotype(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype(individual));
    return rcpp_result_gen;
END_RCPP
}
// count_haplotype_occurrences_individuals
int count_haplotype_occurrences_individuals(const Rcpp::List individuals, const Rcpp::IntegerVector haplotype);
RcppExport SEXP _malan_count_haplotype_occurrences_individuals(SEXP individualsSEXP, SEXP haplotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type haplotype(haplotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(count_haplotype_occurrences_individuals(individuals, haplotype));
    return rcpp_result_gen;
END_RCPP
}
// meiosis_dist_haplotype_matches_individuals
Rcpp::IntegerVector meiosis_dist_haplotype_matches_individuals(const Rcpp::XPtr<Individual> suspect, const Rcpp::List individuals);
RcppExport SEXP _malan_meiosis_dist_haplotype_matches_individuals(SEXP suspectSEXP, SEXP individualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr<Individual> >::type suspect(suspectSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type individuals(individualsSEXP);
    rcpp_result_gen = Rcpp::wrap(meiosis_dist_haplotype_matches_individuals(suspect, individuals));
    return rcpp_result_gen;
END_RCPP
}
// pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists
Rcpp::IntegerMatrix pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(const Rcpp::XPtr<Individual> suspect, int generation_upper_bound_in_result);
RcppExport SEXP _malan_pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(SEXP suspectSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr<Individual> >::type suspect(suspectSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(suspect, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// meiotic_dist
int meiotic_dist(Rcpp::XPtr<Individual> ind1, Rcpp::XPtr<Individual> ind2);
RcppExport SEXP _malan_meiotic_dist(SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(meiotic_dist(ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// count_haplotype_occurrences_pedigree
int count_haplotype_occurrences_pedigree(Rcpp::XPtr<Pedigree> pedigree, const Rcpp::IntegerVector haplotype, int generation_upper_bound_in_result);
RcppExport SEXP _malan_count_haplotype_occurrences_pedigree(SEXP pedigreeSEXP, SEXP haplotypeSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type haplotype(haplotypeSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(count_haplotype_occurrences_pedigree(pedigree, haplotype, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// get_individual
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid);
RcppExport SEXP _malan_get_individual(SEXP populationSEXP, SEXP pidSEXP) {
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
RcppExport SEXP _malan_get_pid(SEXP individualSEXP) {
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
RcppExport SEXP _malan_print_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    print_individual(individual);
    return R_NilValue;
END_RCPP
}
// get_generation
int get_generation(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_get_generation(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_generation(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_from_individual
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_get_pedigree_from_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_from_individual(individual));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_id_from_pid
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids);
RcppExport SEXP _malan_get_pedigree_id_from_pid(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pids(pidsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_id_from_pid(population, pids));
    return rcpp_result_gen;
END_RCPP
}
// count_brothers
int count_brothers(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_count_brothers(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(count_brothers(individual));
    return rcpp_result_gen;
END_RCPP
}
// brothers_matching
int brothers_matching(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_brothers_matching(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(brothers_matching(individual));
    return rcpp_result_gen;
END_RCPP
}
// father_matches
bool father_matches(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_father_matches(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(father_matches(individual));
    return rcpp_result_gen;
END_RCPP
}
// grandfather_matches
bool grandfather_matches(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_grandfather_matches(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(grandfather_matches(individual));
    return rcpp_result_gen;
END_RCPP
}
// count_uncles
int count_uncles(Rcpp::XPtr<Individual> individual);
RcppExport SEXP _malan_count_uncles(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    rcpp_result_gen = Rcpp::wrap(count_uncles(individual));
    return rcpp_result_gen;
END_RCPP
}
// malan_test
void malan_test();
RcppExport SEXP _malan_malan_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    malan_test();
    return R_NilValue;
END_RCPP
}
// pop_size
int pop_size(Rcpp::XPtr<Population> population);
RcppExport SEXP _malan_pop_size(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    rcpp_result_gen = Rcpp::wrap(pop_size(population));
    return rcpp_result_gen;
END_RCPP
}
// meioses_generation_distribution
Rcpp::IntegerMatrix meioses_generation_distribution(Rcpp::XPtr<Individual> individual, int generation_upper_bound_in_result);
RcppExport SEXP _malan_meioses_generation_distribution(SEXP individualSEXP, SEXP generation_upper_bound_in_resultSEXP) {
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
RcppExport SEXP _malan_population_size_generation(SEXP populationSEXP, SEXP generation_upper_bound_in_resultSEXP) {
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
RcppExport SEXP _malan_pedigree_size_generation(SEXP pedigreeSEXP, SEXP generation_upper_bound_in_resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< int >::type generation_upper_bound_in_result(generation_upper_bound_in_resultSEXP);
    rcpp_result_gen = Rcpp::wrap(pedigree_size_generation(pedigree, generation_upper_bound_in_result));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_id
int get_pedigree_id(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _malan_get_pedigree_id(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_id(ped));
    return rcpp_result_gen;
END_RCPP
}
// pedigrees_count
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP _malan_pedigrees_count(SEXP pedigreesSEXP) {
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
RcppExport SEXP _malan_pedigree_size(SEXP pedSEXP) {
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
RcppExport SEXP _malan_pedigrees_table(SEXP pedigreesSEXP) {
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
RcppExport SEXP _malan_get_pedigree(SEXP pedigreesSEXP, SEXP indexSEXP) {
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
RcppExport SEXP _malan_print_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    print_pedigree(ped);
    return R_NilValue;
END_RCPP
}
// get_pids_in_pedigree
Rcpp::IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _malan_get_pids_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pids_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotypes_in_pedigree
Rcpp::List get_haplotypes_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _malan_get_haplotypes_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotypes_in_pedigree(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_edgelist
Rcpp::CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _malan_get_pedigree_edgelist(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_edgelist(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_as_graph
Rcpp::List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP _malan_get_pedigree_as_graph(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_as_graph(ped));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_malan_wipe_pedigrees", (DL_FUNC) &_malan_wipe_pedigrees, 1},
    {"_malan_build_pedigrees", (DL_FUNC) &_malan_build_pedigrees, 2},
    {"_malan_sample_geneology", (DL_FUNC) &_malan_sample_geneology, 9},
    {"_malan_sample_geneology_varying_size", (DL_FUNC) &_malan_sample_geneology_varying_size, 7},
    {"_malan_indices_in_mixture", (DL_FUNC) &_malan_indices_in_mixture, 3},
    {"_malan_pedigree_get_haplotypes_pids", (DL_FUNC) &_malan_pedigree_get_haplotypes_pids, 2},
    {"_malan_individuals_get_haplotypes", (DL_FUNC) &_malan_individuals_get_haplotypes, 1},
    {"_malan_pedigree_populate_haplotypes", (DL_FUNC) &_malan_pedigree_populate_haplotypes, 3},
    {"_malan_pedigrees_all_populate_haplotypes", (DL_FUNC) &_malan_pedigrees_all_populate_haplotypes, 4},
    {"_malan_pedigrees_all_populate_haplotypes_ladder_bounded", (DL_FUNC) &_malan_pedigrees_all_populate_haplotypes_ladder_bounded, 5},
    {"_malan_get_haplotype", (DL_FUNC) &_malan_get_haplotype, 1},
    {"_malan_count_haplotype_occurrences_individuals", (DL_FUNC) &_malan_count_haplotype_occurrences_individuals, 2},
    {"_malan_meiosis_dist_haplotype_matches_individuals", (DL_FUNC) &_malan_meiosis_dist_haplotype_matches_individuals, 2},
    {"_malan_pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists", (DL_FUNC) &_malan_pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists, 2},
    {"_malan_meiotic_dist", (DL_FUNC) &_malan_meiotic_dist, 2},
    {"_malan_count_haplotype_occurrences_pedigree", (DL_FUNC) &_malan_count_haplotype_occurrences_pedigree, 3},
    {"_malan_get_individual", (DL_FUNC) &_malan_get_individual, 2},
    {"_malan_get_pid", (DL_FUNC) &_malan_get_pid, 1},
    {"_malan_print_individual", (DL_FUNC) &_malan_print_individual, 1},
    {"_malan_get_generation", (DL_FUNC) &_malan_get_generation, 1},
    {"_malan_get_pedigree_from_individual", (DL_FUNC) &_malan_get_pedigree_from_individual, 1},
    {"_malan_get_pedigree_id_from_pid", (DL_FUNC) &_malan_get_pedigree_id_from_pid, 2},
    {"_malan_count_brothers", (DL_FUNC) &_malan_count_brothers, 1},
    {"_malan_brothers_matching", (DL_FUNC) &_malan_brothers_matching, 1},
    {"_malan_father_matches", (DL_FUNC) &_malan_father_matches, 1},
    {"_malan_grandfather_matches", (DL_FUNC) &_malan_grandfather_matches, 1},
    {"_malan_count_uncles", (DL_FUNC) &_malan_count_uncles, 1},
    {"_malan_malan_test", (DL_FUNC) &_malan_malan_test, 0},
    {"_malan_pop_size", (DL_FUNC) &_malan_pop_size, 1},
    {"_malan_meioses_generation_distribution", (DL_FUNC) &_malan_meioses_generation_distribution, 2},
    {"_malan_population_size_generation", (DL_FUNC) &_malan_population_size_generation, 2},
    {"_malan_pedigree_size_generation", (DL_FUNC) &_malan_pedigree_size_generation, 2},
    {"_malan_get_pedigree_id", (DL_FUNC) &_malan_get_pedigree_id, 1},
    {"_malan_pedigrees_count", (DL_FUNC) &_malan_pedigrees_count, 1},
    {"_malan_pedigree_size", (DL_FUNC) &_malan_pedigree_size, 1},
    {"_malan_pedigrees_table", (DL_FUNC) &_malan_pedigrees_table, 1},
    {"_malan_get_pedigree", (DL_FUNC) &_malan_get_pedigree, 2},
    {"_malan_print_pedigree", (DL_FUNC) &_malan_print_pedigree, 1},
    {"_malan_get_pids_in_pedigree", (DL_FUNC) &_malan_get_pids_in_pedigree, 1},
    {"_malan_get_haplotypes_in_pedigree", (DL_FUNC) &_malan_get_haplotypes_in_pedigree, 1},
    {"_malan_get_pedigree_edgelist", (DL_FUNC) &_malan_get_pedigree_edgelist, 1},
    {"_malan_get_pedigree_as_graph", (DL_FUNC) &_malan_get_pedigree_as_graph, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_malan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
