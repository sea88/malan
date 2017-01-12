#' @export
print.malan_population_abort <-
  function(x, ...) {
    if (!is(x, "malan_population_abort")) stop("x must be a malan_population_abort object")
    cat("Operation was cancelled, hence the assembly was not finished.\n")    
    return(invisible(NULL))
  }

#' @export
print.malan_population <-
  function(x, ...) {
    if (!is(x, "malan_population")) stop("x must be a malan_population object")
    
    cat("Population with ", formatC(pop_size(x), big.mark = ","), " individuals\n", sep = "")
    
    return(invisible(NULL))
  }

#' @export
print.malan_pedigreelist <-
  function(x, ...) {
    if (!is(x, "malan_pedigreelist")) stop("x must be a malan_pedigreelist object")
    
    sizes <- unlist(lapply(1L:pedigrees_count(x), function(i) pedigree_size(x[[i]])))
    sizes_str <- ""
    
    max_print <- 6L
    
    if (length(sizes) > 0L) {
      if (length(sizes) <= max_print) {
        sizes_str <- paste0(" (of size ", paste0(sizes, collapse = ", "), ")")
      } else {
        sizes_str <- paste0(" (of size ", paste0(sizes[1L:max_print], collapse = ", "), ", ...)")
      }
    }
    
    cat("List of ", formatC(pedigrees_count(x), big.mark = ","), " pedigrees", sizes_str, "\n", sep = "")
    
    return(invisible(NULL))
  }

  
#' @export
`[[.malan_pedigreelist` <- function(x, ...) {
  i <- ..1
  if (length(i) != 1L || !is.integer(i) || i[1L] <= 0L || i > pedigrees_count(x)) {
    stop("Wrong pedigree selected or invalid element selection criteria")
  }
  
  p <- get_pedigree(x, i - 1L) # -1 to go to 0-based indexing
  return(p)
}

#' @export
`[[.malan_population` <- function(x, ...) {
  pid <- ..1
  if (length(pid) != 1L || !is.numeric(pid)) {
    stop("Wrong individual selected or invalid element selection criteria")
  }
  
  if (!is.integer(pid)) {
    pid <- as.integer(pid)
    warning("Converting to integer explicitely (remember L postfix)")
  }
  
  p <- get_individual(x, pid)
  return(p)
}

#' @export
print.malan_pedigree <-
  function(x, ...) {
    if (!is(x, "malan_pedigree")) stop("x must be a malan_pedigree object")
    
    print_pedigree(x)
    
    return(invisible(NULL))
  }

#' @export
pedigree_as_igraph <-
  function(x, ...) {
    if (!is(x, "malan_pedigree")) stop("x must be a malan_pedigree object")
    
    ginfo <- get_pedigree_as_graph(x)
    g <- igraph::graph_from_data_frame(ginfo$edgelist, directed = TRUE, vertices = ginfo$nodes)
    
    #co <- igraph::layout_nicely(g, dim = 2)
    co <- igraph::layout_as_tree(g, mode = "out")
    attr(g, "layout") <- co

    return(g)
  }

#' @export  
plot.malan_pedigree <-
  function(x, ...) {
    if (!is(x, "malan_pedigree")) stop("x must be a malan_pedigree object")
    
    g <- pedigree_as_igraph(x)    
    igraph::plot.igraph(g)
    
    return(invisible(NULL))
    #eturn(g)
  }
  
#' @export  
plot_with_haplotypes <-
  function(x, population, ...) {
    if (!is(x, "malan_pedigree")) stop("x must be a malan_pedigree object")
    
    x_pids <- get_pids_in_pedigree(x)
    haps <- pedigree_get_father_haplotypes_pids(population, x_pids)
    
    #haps_str <- unlist(lapply(haps, paste0, collapse = ","))
    #haps_str <- unlist(lapply(haps, function(h) {
    #  paste0(h, collapse = ",")
    #}))
    haps_str <- unlist(lapply(haps, function(h) {
      paste0(strwrap(paste0(h, collapse = " "), width = 15), collapse = "\n")
      #paste0(h, collapse = ",")
    }))
        
    g <- pedigree_as_igraph(x)
    #igraph::plot.igraph(g, vertex.label = haps_str)
    igraph::plot.igraph(g, vertex.label = haps_str, vertex.label.cex = 0.75, layout = igraph::layout_as_tree(graph = g))
    
    return(invisible(NULL))
    #eturn(g)
  }
  
#' @export  
plot_with_haplotypes_mark_pids <-
  function(x, population, mark_pids = NULL, ...) {
    if (!is(x, "malan_pedigree")) stop("x must be a malan_pedigree object")
    
    x_pids <- get_pids_in_pedigree(x)
    haps <- pedigree_get_father_haplotypes_pids(population, x_pids)
    
    haps_str <- unlist(lapply(seq_along(haps), function(h_i) {
      h <- haps[[h_i]]
      paste0(strwrap(paste0(x_pids[h_i], ":", paste0(h, collapse = " ")), width = 15), collapse = "\n")
    }))
    
    vertex_colors <- rep("orange", length(haps_str))
    if (!is.null(mark_pids)) {
      vertex_colors[x_pids %in% mark_pids] <- "red"
    }
    
    g <- pedigree_as_igraph(x)
    igraph::V(g)$color <- vertex_colors
    igraph::plot.igraph(g, vertex.label = haps_str, vertex.label.cex = 0.75, layout = igraph::layout_as_tree(graph = g))
    
    return(invisible(NULL))
    #eturn(g)
  }
  


