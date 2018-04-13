context("Graph")

test_pop <- test_create_population()
indvs <- get_individuals(test_pop)
peds <- build_pedigrees(test_pop, progress = FALSE)
ped <- peds[[1L]]

if (FALSE) {
  plot(ped)
}

pid1_meiotic_dist_calc <- meioses_generation_distribution(get_individual(test_pop, pid = 1L))
lines <- "generation meioses count
          0          0       1
          0          4       2
          0          6       2
          1          1       1
          1          3       1
          1          5       1
          2          2       1
          2          4       1
          3          3       1"
lines <- gsub("[ ]+", " ", lines)
lines <- gsub("\n ", "\n", lines)
con <- textConnection(lines)
pid1_meiotic_dist_known <- as.matrix(read.csv(con, sep = " "))

test_that("meioses_generation_distribution works", {
  expect_equal(pid1_meiotic_dist_calc, pid1_meiotic_dist_known)
})
