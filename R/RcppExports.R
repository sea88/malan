# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Build pedigrees from (individuals in) a population.
#' 
#' In a newly simulated population, each individual only knows its father and children. 
#' Using this information, this function builds pedigrees.
#' This makes it easier to e.g. population haplotypes, find path between two individuals 
#' (if they are not in the same pedigree, they are not connected).
#' 
#' @param population Population generated by [sample_geneology()] or [sample_geneology_varying_size()].
#' @param progress Show progress.
#' 
#' @return An object with class `malan_pedigreelist` (an internal list of external pointers to pedigrees).
#' 
#' @seealso [sample_geneology()] and [sample_geneology_varying_size()] for simulating populations.
#'
#' @export
build_pedigrees <- function(population, progress = TRUE) {
    .Call('_malan_build_pedigrees', PACKAGE = 'malan', population, progress)
}

#' Simulate a geneology with constant population size.
#' 
#' This function simulates a geneology where the last generation has `population_size` individuals. 
#' 
#' By the backwards simulating process of the Wright-Fisher model, 
#' individuals with no descendants in the end population are not simulated. 
#' If for some reason additional full generations should be simulated, 
#' the number can be specified via the `extra_generations_full` parameter.
#' This can for example be useful if one wants to simulate the 
#' final 3 generations although some of these may not get (male) children.
#' 
#' Let \eqn{\alpha} be the parameter of a symmetric Dirichlet distribution 
#' specifying each man's probability to be the father of an arbitrary 
#' male in the next generation. When \eqn{\alpha = 5}, a man's relative probability 
#' to be the father has 95\% probability to lie between 0.32 and 2.05, compared with a 
#' constant 1 under the standard Wright-Fisher model and the standard deviation in 
#' the number of male offspring per man is 1.10 (standard Wright-Fisher = 1).
#' 
#' This symmetric Dirichlet distribution is implemented by drawing 
#' father (unscaled) probabilities from a Gamma distribution with 
#' parameters `gamma_parameter_shape` and `gamma_parameter_scale` 
#' that are then normalised to sum to 1. 
#' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
#' the following must be used:
#' \eqn{`gamma_parameter_shape` = \alpha}
#' and 
#' \eqn{`gamma_parameter_scale` = 1/\alpha}.
#' 
#' @param population_size The size of the population.
#' @param generations The number of generations to simulate: 
#'        \itemize{
#'           \item -1 for simulate to 1 founder
#'           \item else simulate this number of generations.
#'        }
#' @param extra_generations_full Additional full generations to be simulated.
#' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
#' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
#' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
#' @param progress Show progress.
#' @param individuals_generations_return How many generations back to return (pointers to) individuals for.
#' @param verbose_result Verbose result.
#' 
#' @return A list with the following entries:
#' \itemize{
#'   \item `population`. An external pointer to the population.
#'   \item `generations`. Generations actually simulated, mostly useful when parameter `generations = -1`.
#'   \item `founders`. Number of founders after the simulated `generations`.
#'   \item `growth_type`. Growth type model.
#'   \item `sdo_type`. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on `enable_gamma_variance_extension`.
#'   \item `end_generation_individuals`. Pointers to individuals in end generation.
#'   \item `individuals_generations`. Pointers to individuals in end generation in addition to the previous `individuals_generations_return`.
#' }
#' If `verbose_result` is true, then these additional components are also returned:
#' \itemize{
#'   \item `individual_pids`. A matrix with pid (person id) for each individual.
#'   \item `father_pids`. A matrix with pid (person id) for each individual's father.
#'   \item `father_indices`. A matrix with indices for fathers.
#' }
#' 
#' @seealso [sample_geneology_varying_size()].
#' 
#' @import Rcpp
#' @import RcppProgress
#' @import RcppArmadillo
#' @export
sample_geneology <- function(population_size, generations, extra_generations_full = 0L, gamma_parameter_shape = 5.0, gamma_parameter_scale = 1.0/5.0, enable_gamma_variance_extension = FALSE, progress = TRUE, individuals_generations_return = 2L, verbose_result = FALSE) {
    .Call('_malan_sample_geneology', PACKAGE = 'malan', population_size, generations, extra_generations_full, gamma_parameter_shape, gamma_parameter_scale, enable_gamma_variance_extension, progress, individuals_generations_return, verbose_result)
}

#' Simulate a geneology with varying population size.
#' 
#' This function simulates a geneology with varying population size specified
#' by a vector of population sizes, one for each generation. 
#' 
#' By the backwards simulating process of the Wright-Fisher model, 
#' individuals with no descendants in the end population are not simulated 
#' If for some reason additional full generations should be simulated, 
#' the number can be specified via the \code{extra_generations_full} parameter.
#' This can for example be useful if one wants to simulate the 
#' final 3 generations although some of these may not get (male) children.
#' 
#' Let \eqn{\alpha} be the parameter of a symmetric Dirichlet distribution 
#' specifying each man's probability to be the father of an arbitrary 
#' male in the next generation. When \eqn{\alpha = 5}, a man's relative probability 
#' to be the father has 95\% probability to lie between 0.32 and 2.05, compared with a 
#' constant 1 under the standard Wright-Fisher model and the standard deviation in 
#' the number of male offspring per man is 1.10 (standard Wright-Fisher = 1).
#' 
#' This symmetric Dirichlet distribution is implemented by drawing 
#' father (unscaled) probabilities from a Gamma distribution with 
#' parameters `gamma_parameter_shape` and `gamma_parameter_scale` 
#' that are then normalised to sum to 1. 
#' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
#' the following must be used:
#' \eqn{`gamma_parameter_shape` = \alpha}
#' and 
#' \eqn{`gamma_parameter_scale` = 1/\alpha}.
#' 
#' @param population_sizes The size of the population at each generation, `g`. 
#'        `population_sizes[g]` is the population size at generation `g`.
#'        The length of population_sizes is the number of generations being simulated.
#' @param extra_generations_full Additional full generations to be simulated.
#' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
#' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
#' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
#' @param progress Show progress.
#' @param individuals_generations_return How many generations back to return (pointers to) individuals for.
#' 
#' @return A malan_simulation / list with the following entries:
#' \itemize{
#'   \item \code{population}. An external pointer to the population.
#'   \item \code{generations}. Generations actually simulated, mostly useful when parameter \code{generations = -1}.
#'   \item \code{founders}. Number of founders after the simulated \code{generations}.
#'   \item \code{growth_type}. Growth type model.
#'   \item \code{sdo_type}. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on \code{enable_gamma_variance_extension}.
#'   \item \code{end_generation_individuals}. Pointers to individuals in end generation.
#'   \item \code{individuals_generations}. Pointers to individuals in end generation in addition to the previous \code{individuals_generations_return}.
#' }
#'
#' @seealso [sample_geneology()].
#' 
#' @import Rcpp
#' @import RcppProgress
#' @import RcppArmadillo
#' @export
sample_geneology_varying_size <- function(population_sizes, extra_generations_full = 0L, gamma_parameter_shape = 5.0, gamma_parameter_scale = 1.0/5.0, enable_gamma_variance_extension = FALSE, progress = TRUE, individuals_generations_return = 2L) {
    .Call('_malan_sample_geneology_varying_size', PACKAGE = 'malan', population_sizes, extra_generations_full, gamma_parameter_shape, gamma_parameter_scale, enable_gamma_variance_extension, progress, individuals_generations_return)
}

#' Calculate genotype probabilities with theta
#' 
#' @param allele_dist Allele distribution (probabilities) -- gets normalised
#' @param theta Theta correction between 0 and 1 (both included)
#'
#' @export
calc_autosomal_genotype_probs <- function(allele_dist, theta) {
    .Call('_malan_calc_autosomal_genotype_probs', PACKAGE = 'malan', allele_dist, theta)
}

#' Calculate conditional genotype cumulative probabilities with theta
#' 
#' @param allele_dist Allele distribution (probabilities) -- gets normalised
#' @param theta Theta correction between 0 and 1 (both included)
#' 
#' @return Matrix: row i: conditional cumulative distribution of alleles given allele i
#'
#' @export
calc_autosomal_genotype_conditional_cumdist <- function(allele_dist, theta) {
    .Call('_malan_calc_autosomal_genotype_conditional_cumdist', PACKAGE = 'malan', allele_dist, theta)
}

#' Sample genotype with theta
#' 
#' @param allele_dist Allele distribution (probabilities) -- gets normalised
#' @param theta Theta correction between 0 and 1 (both included)
#'
#' @export
sample_autosomal_genotype <- function(allele_dist, theta) {
    .Call('_malan_sample_autosomal_genotype', PACKAGE = 'malan', allele_dist, theta)
}

#' Populate 1-locus autosomal DNA profile in pedigrees.
#' 
#' Populate 1-locus autosomal DNA profile from founder and down in all pedigrees.
#' Note, that only alleles from ladder is assigned and 
#' that all founders draw type randomly.
#' 
#' Note, that pedigrees must first have been inferred by [build_pedigrees()].
#' 
#' @param pedigrees Pedigree list in which to populate haplotypes
#' @param allele_dist Allele distribution (probabilities) -- gets normalised
#' @param theta Theta correction between 0 and 1 (both included)
#' @param mutation_rate Mutation rate between 0 and 1 (both included)
#' @param progress Show progress
#'
#' @seealso [pedigrees_all_populate_haplotypes_custom_founders()] and 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @export
pedigrees_all_populate_autosomal <- function(pedigrees, allele_dist, theta, mutation_rate, progress = TRUE) {
    invisible(.Call('_malan_pedigrees_all_populate_autosomal', PACKAGE = 'malan', pedigrees, allele_dist, theta, mutation_rate, progress))
}

hash_colisions <- function(p) {
    .Call('_malan_hash_colisions', PACKAGE = 'malan', p)
}

#' Estimate theta from genetypes
#' 
#' Estimate theta for one subpopulation given a sample of genotypes.
#' 
#' @param x Matrix of genotypes: two columns (allele1 and allele2) and a row per individual
#' 
#' @return List:
#' * estimate: Vector of length 1 containing estimate of theta or NA if it could not be estimated
#' * error: true if an error happened, false otherwise
#' * details: contains description if an error happened
#' 
#' @export
estimate_theta_1subpop_sample <- function(x) {
    .Call('_malan_estimate_theta_1subpop_sample', PACKAGE = 'malan', x)
}

#' Estimate theta from individuals
#' 
#' Estimate theta for one subpopulation given a sample of genotypes.
#' 
#' @param individuals Individuals to get haplotypes for.
#' 
#' @return List:
#' * estimate: Vector of length 1 containing estimate of theta or NA if it could not be estimated
#' * error: true if an error happened, false otherwise
#' * details: contains description if an error happened
#' 
#' @export
estimate_theta_1subpop_individuals <- function(individuals) {
    .Call('_malan_estimate_theta_1subpop_individuals', PACKAGE = 'malan', individuals)
}

#' Populate haplotypes in pedigrees (0-founder/unbounded).
#' 
#' Populate haplotypes from founder and down in all pedigrees.
#' Note, that haplotypes are unbounded and 
#' that all founders get haplotype `rep(0L, loci)`.
#' 
#' Note, that pedigrees must first have been inferred by [build_pedigrees()].
#' 
#' @param pedigrees Pedigree list in which to populate haplotypes
#' @param loci Number of loci
#' @param mutation_rates Vector with mutation rates, length `loci`
#' @param progress Show progress
#'
#' @seealso [pedigrees_all_populate_haplotypes_custom_founders()] and 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @export
pedigrees_all_populate_haplotypes <- function(pedigrees, loci, mutation_rates, progress = TRUE) {
    invisible(.Call('_malan_pedigrees_all_populate_haplotypes', PACKAGE = 'malan', pedigrees, loci, mutation_rates, progress))
}

#' Populate haplotypes in pedigrees (custom founder/unbounded).
#' 
#' Populate haplotypes from founder and down in all pedigrees.
#' Note, that haplotypes are unbounded.
#' All founders get a haplotype from calling the user 
#' provided function `get_founder_haplotype()`.
#' 
#' Note, that pedigrees must first have been inferred by [build_pedigrees()].
#' 
#' @param pedigrees Pedigree list in which to populate haplotypes
#' @param mutation_rates Vector with mutation rates
#' @param get_founder_haplotype Function taking no arguments returning a haplotype of `length(mutation_rates)`
#' @param progress Show progress
#'
#' @seealso [pedigrees_all_populate_haplotypes()] and 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @export
pedigrees_all_populate_haplotypes_custom_founders <- function(pedigrees, mutation_rates, get_founder_haplotype = NULL, progress = TRUE) {
    invisible(.Call('_malan_pedigrees_all_populate_haplotypes_custom_founders', PACKAGE = 'malan', pedigrees, mutation_rates, get_founder_haplotype, progress))
}

#' Populate haplotypes in pedigrees (custom founder/bounded).
#' 
#' Populate haplotypes from founder and down in all pedigrees.
#' Note, that haplotypes are bounded by `ladder_min` and `ladder_max`.
#' All founders get a haplotype from calling the user 
#' provided function `get_founder_haplotype()`.
#' 
#' Note, that pedigrees must first have been inferred by [build_pedigrees()].
#' 
#' @param pedigrees Pedigree list in which to populate haplotypes
#' @param mutation_rates Vector with mutation rates
#' @param ladder_min Lower bounds for haplotypes, same length as `mutation_rates`
#' @param ladder_max Upper bounds for haplotypes, same length as `mutation_rates`; all entries must be strictly greater than `ladder_min`
#' @param get_founder_haplotype Function taking no arguments returning a haplotype of `length(mutation_rates)`
#' @param progress Show progress
#'
#' @seealso [pedigrees_all_populate_haplotypes()] and 
#' [pedigrees_all_populate_haplotypes_custom_founders()].
#' 
#' @export
pedigrees_all_populate_haplotypes_ladder_bounded <- function(pedigrees, mutation_rates, ladder_min, ladder_max, get_founder_haplotype = NULL, progress = TRUE) {
    invisible(.Call('_malan_pedigrees_all_populate_haplotypes_ladder_bounded', PACKAGE = 'malan', pedigrees, mutation_rates, ladder_min, ladder_max, get_founder_haplotype, progress))
}

#' Get haplotype from an individual
#' 
#' Requires that haplotypes are first populated, e.g. 
#' with [pedigrees_all_populate_haplotypes()], 
#' [pedigrees_all_populate_haplotypes_custom_founders()], or 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @param individual Individual to get haplotypes for.
#' @return Haplotype for `individual`.
#' 
#' @seealso [get_haplotypes_individuals()] and [get_haplotypes_pids()].
#' 
#' @export
get_haplotype <- function(individual) {
    .Call('_malan_get_haplotype', PACKAGE = 'malan', individual)
}

#' Get haplotype matrix from list of individuals
#' 
#' Requires that haplotypes are first populated, e.g. 
#' with [pedigrees_all_populate_haplotypes()], 
#' [pedigrees_all_populate_haplotypes_custom_founders()], or 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @param individuals Individuals to get haplotypes for.
#' @return Matrix of haplotypes where row `i` is the haplotype of `individuals[[i]]`.
#' 
#' @seealso [get_haplotypes_pids()].
#' 
#' @export
get_haplotypes_individuals <- function(individuals) {
    .Call('_malan_get_haplotypes_individuals', PACKAGE = 'malan', individuals)
}

#' Get haplotypes from a vector of pids.
#' 
#' Requires that haplotypes are first populated, e.g. 
#' with [pedigrees_all_populate_haplotypes()], 
#' [pedigrees_all_populate_haplotypes_custom_founders()], or 
#' [pedigrees_all_populate_haplotypes_ladder_bounded()].
#' 
#' @param population Population
#' @param pids Vector of pids to get haplotypes for.
#' 
#' @return Matrix of haplotypes where row `i` is the haplotype of `individuals[[i]]`.
#' 
#' @seealso [get_haplotypes_individuals()].
#' 
#' @export
get_haplotypes_pids <- function(population, pids) {
    .Call('_malan_get_haplotypes_pids', PACKAGE = 'malan', population, pids)
}

#' Count haplotypes occurrences in list of individuals
#' 
#' Counts the number of types `haplotype` appears in `individuals`.
#' 
#' @param individuals List of individuals to count occurrences in.
#' @param haplotype Haplotype to count occurrences of.
#' 
#' @return Number of times that `haplotype` occurred amongst `individuals`.
#' 
#' @seealso [pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists()].
#' 
#' @export
count_haplotype_occurrences_individuals <- function(individuals, haplotype) {
    .Call('_malan_count_haplotype_occurrences_individuals', PACKAGE = 'malan', individuals, haplotype)
}

#' Get individuals matching from list of individuals
#' 
#' Get the indvididuals that matches `haplotype` in `individuals`.
#' 
#' @param individuals List of individuals to count occurrences in.
#' @param haplotype Haplotype to count occurrences of.
#' 
#' @return List of individuals that matches `haplotype` amongst `individuals`.
#' 
#' @seealso [pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists()].
#' 
#' @export
haplotype_matches_individuals <- function(individuals, haplotype) {
    .Call('_malan_haplotype_matches_individuals', PACKAGE = 'malan', individuals, haplotype)
}

#' Count haplotypes occurrences in pedigree
#' 
#' Counts the number of types `haplotype` appears in `pedigree`.
#' 
#' @param pedigree Pedigree to count occurrences in.
#' @param haplotype Haplotype to count occurrences of.
#' @param generation_upper_bound_in_result Only consider matches in 
#' generation 0, 1, ... generation_upper_bound_in_result.
#' -1 means disabled, consider all generations.
#' End generation is generation 0.
#' Second last generation is 1. 
#' And so on.
#' 
#' @return Number of times that `haplotype` occurred in `pedigree`.
#' 
#' @seealso [pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists()].
#' 
#' @export
count_haplotype_occurrences_pedigree <- function(pedigree, haplotype, generation_upper_bound_in_result = -1L) {
    .Call('_malan_count_haplotype_occurrences_pedigree', PACKAGE = 'malan', pedigree, haplotype, generation_upper_bound_in_result)
}

#' Information about matching individuals
#' 
#' Gives information about all individuals in pedigree that matches an individual.
#' Just as [count_haplotype_occurrences_individuals()] counts the number of 
#' occurrences amongst a list of individuals, 
#' this gives detailed information about matching individuals in the pedigree, 
#' e.g. meiotic distances and maximum L1 distance on the path as some of these 
#' matches may have (back)mutations between in between them (but often this will be 0).
#' 
#' @param suspect Individual that others must match the profile of.
#' @param generation_upper_bound_in_result Only consider matches in 
#' generation 0, 1, ... generation_upper_bound_in_result.
#' -1 means disabled, consider all generations.
#' End generation is generation 0.
#' Second last generation is 1. 
#' And so on.
#' 
#' @return Matrix with information about matching individuals. 
#' Columns in order: meioses (meiotic distance to `suspect`), 
#' max_L1 (on the path between the matching individual and `suspect`, 
#' what is the maximum L1 distance between the `suspect`'s profile and the 
#' profiles of the individuals on the path), 
#' pid (pid of matching individual)
#' 
#' @seealso [count_haplotype_occurrences_individuals()].
#'
#' @export
pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists <- function(suspect, generation_upper_bound_in_result = -1L) {
    .Call('_malan_pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists', PACKAGE = 'malan', suspect, generation_upper_bound_in_result)
}

#' Meiotic distance between two individuals
#' 
#' Get the number of meioses between two individuals.
#' Note, that pedigrees must first have been inferred by [build_pedigrees()].
#' 
#' @param ind1 Individual 1
#' @param ind2 Individual 2
#' 
#' @return Number of meioses between `ind1` and `ind2` if they are in the same pedigree, else -1.
#' 
#' @export
meiotic_dist <- function(ind1, ind2) {
    .Call('_malan_meiotic_dist', PACKAGE = 'malan', ind1, ind2)
}

#' Get individual by pid
#' 
#' @param population Population
#' @param pid pid
#' 
#' @return Individual
#' 
#' @export
get_individual <- function(population, pid) {
    .Call('_malan_get_individual', PACKAGE = 'malan', population, pid)
}

#' Get pid from individual
#' 
#' @param individual Individual to get pid of
#' 
#' @return pid
#' 
#' @export
get_pid <- function(individual) {
    .Call('_malan_get_pid', PACKAGE = 'malan', individual)
}

#' Print individual
#' 
#' @param individual Individual
#' 
#' @export
print_individual <- function(individual) {
    invisible(.Call('_malan_print_individual', PACKAGE = 'malan', individual))
}

#' Get individual's generation number
#' 
#' Note that generation 0 is final, end generation. 
#' 1 is second last generation etc.
#' 
#' @param individual Individual
#' 
#' @return generation
#' 
#' @export
get_generation <- function(individual) {
    .Call('_malan_get_generation', PACKAGE = 'malan', individual)
}

#' Get pedigree from individual
#' 
#' @param individual Individual
#' 
#' @return pedigree
#' 
#' @export
get_pedigree_from_individual <- function(individual) {
    .Call('_malan_get_pedigree_from_individual', PACKAGE = 'malan', individual)
}

#' Get pedigree ids from pids
#'
#' @param population Population
#' @param pids Pids
#' 
#' @return Vector with pedigree ids
#' 
#' @export
get_pedigree_id_from_pid <- function(population, pids) {
    .Call('_malan_get_pedigree_id_from_pid', PACKAGE = 'malan', population, pids)
}

#' Get individual's family information
#'
#' @param individual individual
#' 
#' @return List with family information
#' 
#' @export
get_family_info <- function(individual) {
    .Call('_malan_get_family_info', PACKAGE = 'malan', individual)
}

#' Number of brothes
#' 
#' Get individual's number of brothes
#'
#' @param individual individual
#' 
#' @return Number of brothers
#' 
#' @export
count_brothers <- function(individual) {
    .Call('_malan_count_brothers', PACKAGE = 'malan', individual)
}

#' Number of brothes with matching haplotype
#' 
#' Get individual's number of brothes that matches `individual`'s haplotype
#'
#' @param individual individual
#' 
#' @return Number of brothers that matches `individual`'s haplotype
#' 
#' @export
brothers_matching <- function(individual) {
    .Call('_malan_brothers_matching', PACKAGE = 'malan', individual)
}

#' Father matches
#' 
#' Does the father have the same profile as `individual`?
#'
#' @param individual individual
#' 
#' @return Whether father has the same profile as `individual` or not
#' 
#' @export
father_matches <- function(individual) {
    .Call('_malan_father_matches', PACKAGE = 'malan', individual)
}

#' Grandfather matches
#' 
#' Does the frandfather have the same profile as `individual`?
#'
#' @param individual individual
#' 
#' @return Whether grandfather has the same profile as `individual` or not
#' 
#' @export
grandfather_matches <- function(individual) {
    .Call('_malan_grandfather_matches', PACKAGE = 'malan', individual)
}

#' Number of uncles
#' 
#' Get individual's number of uncles
#'
#' @param individual individual
#' 
#' @return Number of uncles
#' 
#' @export
count_uncles <- function(individual) {
    .Call('_malan_count_uncles', PACKAGE = 'malan', individual)
}

pop_size <- function(population) {
    .Call('_malan_pop_size', PACKAGE = 'malan', population)
}

#' Get all individuals in population
#' 
#' @param population Population
#'
#' @export
get_individuals <- function(population) {
    .Call('_malan_get_individuals', PACKAGE = 'malan', population)
}

#' Meiotic distribution
#' 
#' Get the distribution of number of meioses from `individual` 
#' to all individuals in `individual`'s pedigree.
#' Note the `generation_upper_bound_in_result` parameter.
#' 
#' @param individual Individual to calculate all meiotic distances from
#' @param generation_upper_bound_in_result Limit on distribution; -1 means no limit. 
#' 0 is the final generation. 1 second last generation etc.
#' 
#' @export
meioses_generation_distribution <- function(individual, generation_upper_bound_in_result = -1L) {
    .Call('_malan_meioses_generation_distribution', PACKAGE = 'malan', individual, generation_upper_bound_in_result)
}

#' Size of population
#' 
#' Get the size of the population.
#' Note the `generation_upper_bound_in_result` parameter.
#' 
#' @param population Population to get size of
#' @param generation_upper_bound_in_result Limit on generation to include in count; -1 means no limit. 
#' 0 only include the final generation. 1 only second last generation etc.
#' 
#' @export
population_size_generation <- function(population, generation_upper_bound_in_result = -1L) {
    .Call('_malan_population_size_generation', PACKAGE = 'malan', population, generation_upper_bound_in_result)
}

#' Size of pedigree
#' 
#' Get the size of the pedigree.
#' Note the `generation_upper_bound_in_result` parameter.
#' 
#' @param pedigree Pedigree to get size of
#' @param generation_upper_bound_in_result Limit on generation to include in count; -1 means no limit. 
#' 0 only include the final generation. 1 only second last generation etc.
#' 
#' @export
pedigree_size_generation <- function(pedigree, generation_upper_bound_in_result = -1L) {
    .Call('_malan_pedigree_size_generation', PACKAGE = 'malan', pedigree, generation_upper_bound_in_result)
}

#' Mixture information about 2 persons' mixture of donor1 and donor2.
#' 
#' @param individuals Individuals to consider as possible contributors and thereby get information from.
#' @param donor1 Contributor1/donor 1
#' @param donor2 Contributor2/donor 2
#' @return A list with mixture information about the mixture \code{donor1}+\code{donor2}+\code{donor3} from \code{individuals}
#' 
#' @seealso \code{\link{mixture_info_by_individuals_3pers}}
#' 
#' @export
mixture_info_by_individuals <- function(individuals, donor1, donor2) {
    .Call('_malan_mixture_info_by_individuals', PACKAGE = 'malan', individuals, donor1, donor2)
}

#' Mixture information about 3 persons' mixture of donor1, donor2 and donor3.
#' 
#' @inherit mixture_info_by_individuals
#' @param donor3 Contributor2/donor 3
#' 
#' @seealso \code{\link{mixture_info_by_individuals}}
#' 
#' @export
mixture_info_by_individuals_3pers <- function(individuals, donor1, donor2, donor3) {
    .Call('_malan_mixture_info_by_individuals_3pers', PACKAGE = 'malan', individuals, donor1, donor2, donor3)
}

#' Get pedigree id
#' 
#' @param ped Pedigree
#' 
#' @export
get_pedigree_id <- function(ped) {
    .Call('_malan_get_pedigree_id', PACKAGE = 'malan', ped)
}

#' Get number of pedigrees
#' 
#' @param pedigrees Pedigrees
#' 
#' @export
pedigrees_count <- function(pedigrees) {
    .Call('_malan_pedigrees_count', PACKAGE = 'malan', pedigrees)
}

#' Get pedigree size
#' 
#' @param ped Pedigree
#' 
#' @export
pedigree_size <- function(ped) {
    .Call('_malan_pedigree_size', PACKAGE = 'malan', ped)
}

#' Get distribution of pedigree sizes
#' 
#' @param pedigrees Pedigrees
#' 
#' @export
pedigrees_table <- function(pedigrees) {
    .Call('_malan_pedigrees_table', PACKAGE = 'malan', pedigrees)
}

get_pedigree <- function(pedigrees, index) {
    .Call('_malan_get_pedigree', PACKAGE = 'malan', pedigrees, index)
}

print_pedigree <- function(ped) {
    invisible(.Call('_malan_print_pedigree', PACKAGE = 'malan', ped))
}

#' Get pids in pedigree
#' 
#' @param ped Pedigree
#' 
#' @export
get_pids_in_pedigree <- function(ped) {
    .Call('_malan_get_pids_in_pedigree', PACKAGE = 'malan', ped)
}

#' Get haplotypes in pedigree
#' 
#' @param ped Pedigree
#' 
#' @return List with haplotypes
#' 
#' @export
get_haplotypes_in_pedigree <- function(ped) {
    .Call('_malan_get_haplotypes_in_pedigree', PACKAGE = 'malan', ped)
}

get_pedigree_edgelist <- function(ped) {
    .Call('_malan_get_pedigree_edgelist', PACKAGE = 'malan', ped)
}

#' Get pedigree information as graph (mainly intended for plotting)
#' 
#' @param ped Pedigree
#' 
#' @export
get_pedigree_as_graph <- function(ped) {
    .Call('_malan_get_pedigree_as_graph', PACKAGE = 'malan', ped)
}

#' Get pedigrees information in tidy format
#' 
#' @param pedigrees Pedigrees
#' 
get_pedigrees_tidy <- function(pedigrees) {
    .Call('_malan_get_pedigrees_tidy', PACKAGE = 'malan', pedigrees)
}

#' Generate test population
#' 
#' @return An external pointer to the population.
test_create_population <- function() {
    .Call('_malan_test_create_population', PACKAGE = 'malan')
}

