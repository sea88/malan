context("Autosomal")

ESTIMATION_TOL <- 0.01

set.seed(1)
sim_res_fixed <- sample_geneology(population_size = 1e3, 
                                  generations = 100, 
                                  extra_generations_full = 2,
                                  individuals_generations_return = 2, # default value
                                  progress = FALSE)

peds <- build_pedigrees(sim_res_fixed$population, progress = FALSE)

# alleles         A    B    C
allele_prob <- c(0.7, 0.2, 0.1)
theta <- 0.1
L <- length(allele_prob)

x <- replicate(10000, sample_autosomal_genotype(allele_prob, theta))
x_p_tab <- prop.table(table(c(x)))
x_p <- as.numeric(x_p_tab)

test_that("sample_autosomal_genotype works", {
  expect_equal(x_p, allele_prob, tol = ESTIMATION_TOL)
})


geno_probs_R_mat <- matrix(0, nrow = L, ncol = L)  
for (i in 1L:L) {
  for (j in i:L) {
    #P_AA = F*p_A + (1 - F) * p_A^2
    #X != A: P_AX = (1 - F) * 2 * p_A*p_X
    prob <- if (i == j) # homzyg:
      theta*allele_prob[i] + (1-theta)*allele_prob[i]^2
    else
      (1-theta)*2*allele_prob[i]*allele_prob[j]
    
    # only assign (i, j);
    # if (j, i) should also be assigned, the factor 2 above in
    # the else case must be removed
    geno_probs_R_mat[i, j] <- prob
  }
}
geno_probs_R <- geno_probs_R_mat[upper.tri(geno_probs_R_mat, diag = TRUE)]
geno_probs_rcpp <- calc_autosomal_genotype_probs(allele_prob, theta)
test_that("calc_autosomal_genotype_probs works", {
  expect_equal(geno_probs_R, geno_probs_rcpp)
})

pedigrees_all_populate_autosomal(peds, allele_prob, theta, mutation_rate = 0, FALSE)

geno_mat <- do.call(rbind, lapply(seq_len(pedigrees_count(peds)), function(i) {
  hs <- do.call(rbind, get_haplotypes_in_pedigree(peds[[i]]))
  hs
}))
geno_vec <- c(geno_mat)
geno_vec_p <- as.numeric(prop.table(table(geno_vec)))

test_that("pedigrees_all_populate_autosomal works", {
  expect_equal(geno_vec_p, allele_prob, tol = ESTIMATION_TOL)
})



y <- t(x)
x_unique_geno <- y[!duplicated(y), ]
ols_quantities <- apply(x_unique_geno, 1,
                        function(as) {
                          geno_prob <- mean(apply(y, 1, function(as2) all(as == as2)))
                          
                          if (length(unique(as)) == 1L) {
                            p <- x_p_tab[as.character(as[1])]
                            
                            x <- p - p^2
                            y <- geno_prob - p^2
                            return(c(x, y))
                          } else {
                            p <- x_p_tab[as.character(as[1])]
                            q <- x_p_tab[as.character(as[2])]
                            
                            x <- -2*p*q
                            y <- geno_prob - 2*p*q
                            return(c(x, y))
                          }
                        })
ols_x <- as.matrix(ols_quantities[1L, ])
ols_y <- ols_quantities[2L, ]

test_that("estimate_theta_1subpop works", {
  expect_equal(qr.solve(ols_x, ols_y), estimate_theta_1subpop(y))
})


# bootstrap:
theta_boot <- replicate(100, {
  yboot <- y[sample(nrow(y), replace = TRUE), ]
  estimate_theta_1subpop(yboot)
})
# Expanding a bit...
theta_boot_rng <- c(0.9, 1.1) * range(theta_boot)
theta_boot_rng

test_that("estimate_theta_1subpop boot contains true", {
  expect_true(all(theta_boot >= 0 & theta_boot <= 1))
  expect_true(theta >= theta_boot_rng[1L] & theta <= theta_boot_rng[2L])
})

