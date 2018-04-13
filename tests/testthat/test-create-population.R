context("Pedigrees and haplotypes")

test_pop <- test_create_population()

test_that("test_create_population works", {
  expect_failure(expect_null(test_pop))
  expect_output(print(test_pop), regexp = "^Population with 11 individuals$")
  expect_equal(pop_size(test_pop), 11)
})

indvs <- get_individuals(test_pop)
test_that("get_individuals works", {
  expect_failure(expect_null(indvs))
  expect_equal(length(indvs), 11L)
})

peds <- build_pedigrees(test_pop, progress = FALSE)
test_that("build_pedigrees works", {
  expect_output(print(peds), regexp = "^List of 1 pedigrees \\(of size 11\\)$")
  expect_equal(pedigrees_count(peds), 1L)
})
ped <- peds[[1L]]
pids <- sort(get_pids_in_pedigree(ped))

test_that("pedigree pids works", {
  expect_equal(length(pids), 11L)
  expect_true(all(pids == 1L:11L))
})

LOCI <- 5L

pedigrees_all_populate_haplotypes(peds, loci = LOCI, mutation_rates = rep(0, LOCI), progress = FALSE)
test_that("pedigrees_all_populate_haplotypes works", {
  expect_output(print(peds), regexp = "^List of 1 pedigrees \\(of size 11\\)$")
})

haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes haplotypes works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 11L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_ped) == 0L))
  
  #haps_from_pids
  expect_true(is.matrix(haps_from_pids))
  expect_equal(nrow(haps_from_pids), 11L)
  expect_equal(ncol(haps_from_pids), LOCI)
  
  #haps_from_indvs
  expect_equal(haps_from_indvs, haps_from_pids)
  
  #hap_from_indv
  expect_equal(haps_from_ped, hap_from_indv)
})




test_that("pedigrees_all_populate_haplotypes_ladder_bounded error works", {
  #haps_from_ped
  expect_error(
    pedigrees_all_populate_haplotypes_ladder_bounded(peds, 
                                                     mutation_rates = rep(1, LOCI), 
                                                     ladder_min = rep(10L, LOCI), 
                                                     ladder_max = rep(10L, LOCI), 
                                                     get_founder_haplotype = function() rep(10L, LOCI),
                                                     progress = FALSE)
  )
})


pedigrees_all_populate_haplotypes_ladder_bounded(peds, 
                                                 mutation_rates = rep(0, LOCI), 
                                                 ladder_min = rep(10L, LOCI), 
                                                 ladder_max = rep(11L, LOCI), 
                                                 get_founder_haplotype = function() rep(10L, LOCI),
                                                 progress = FALSE)

haps_from_ped <- get_haplotypes_in_pedigree(ped)
haps_from_pids <- get_haplotypes_pids(test_pop, pids)
haps_from_indvs <- get_haplotypes_individuals(indvs)
hap_from_indv <- lapply(pids, function(pid) get_haplotype(get_individual(test_pop, pid)))

test_that("pedigrees_all_populate_haplotypes_ladder_bounded haplotypes works", {
  #haps_from_ped
  expect_true(is.list(haps_from_ped))
  expect_equal(length(haps_from_ped), 11L)
  expect_equal(length(haps_from_ped[[1L]]), LOCI)
  expect_true(all(unlist(haps_from_ped) == 10L))
  
  #haps_from_pids
  expect_true(is.matrix(haps_from_pids))
  expect_equal(nrow(haps_from_pids), 11L)
  expect_equal(ncol(haps_from_pids), LOCI)
  
  #haps_from_indvs
  expect_equal(haps_from_indvs, haps_from_pids)
  
  #hap_from_indv
  expect_equal(haps_from_ped, hap_from_indv)
})


