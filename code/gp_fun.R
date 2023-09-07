###################################
## function to generate proposal ##
## n: number of times to sample iteratively
generate_proposal <- function(x) {
  rand.data <- data.frame(celltype = sample(x = celltype, size = length(loci), replace = TRUE),
                          locus = sample(x = loci, size = length(loci), replace = FALSE))
  ## Randomly select gene name by cytoband
  gene <- loci2gene %>% group_by(cytoband) %>% slice_sample(n = 1) %>% dplyr::rename("gene"="gene_name")
  ## Append gene name to rand.data
  rand.data <- rand.data %>% left_join(x = ., y = gene, by = c("locus"="cytoband"))
}


######################################
## function to progressive sampling ##
progressive_quantile <- function(group.prob, n_proposal, proposal) {
  ## Probabilistic
  ## Sample row number based on probs of ranked df
  rand.proposal <- sample(x = 1:nrow(proposal), ## Proposal should be ranked high to low
                          size = 100, ## Select 100 random proposals based on the prob below
                          prob = rep(rev(group.prob), each = nrow(proposal)/5))
  proposal <- proposal[rand.proposal,]
  }

mating.pair <- function(sel.proposal, n.pairs = 100, n_proposal = 1000, df, khan.method, fitness.df, top.df, generation, n.loci) {
  if(khan.method == TRUE) {
    ## Get khan info
    khan <- fitness.df[1,] ## Which gen and proposal is khan?
    khan.df <- top.df %>% dplyr::filter(n == khan$n, ## Get khan proposals
                                        gen == khan$gen) %>%
      dplyr::select(-gen) %>%
      mutate(n = 1000 + generation)
    ############################
    ## Preparing proposal set ##
    ## Get row.num of random group 5 to kick out
    rand.kick <- sel.proposal %>% rownames_to_column(var = "row") %>% dplyr::filter(quantile == 5) %>% slice_sample(n = 1) %>% .$row
    proposal.set <- sel.proposal %>% dplyr::filter(rownames(.) != rand.kick)
    ## Add khan into proposal set
    khan <- khan %>% dplyr::rename("quantile"="gen") %>% mutate(n = 1000 + generation)
    proposal.set <- rbind(khan, proposal.set)
    ######################
    ## Preparing df set ##
    df.set <- rbind(df, khan.df)
  } else {
    proposal.set <- sel.proposal
    df.set <- df
  }
  ## Get set of mating pairs
  mates <- replicate(n = n.pairs, sample(x = proposal.set$n, size = 2, replace = FALSE)) %>% as.data.frame() %>% as.list()

  mating <- mclapply(X = mates, FUN = function(x) {
    child.list <- list()
    for(n in 1:(n_proposal/n.pairs)) {
      ## Parent1: Randomly select 50% proposal elements
      p1 <- df.set %>% dplyr::filter(n == x[1]) %>% .[sample(x = nrow(.), size = ceiling(n.loci*0.5)),] %>%
        dplyr::select(-n)
      ## Parent2: Remove locus selected in P1
      p2 <- df.set %>% dplyr::filter(n == x[2],
                                 !locus %in% p1$locus) %>%
        dplyr::select(-n)
      ## Combine P1 and P2 to create child
      child <- rbind(p1,p2) %>% remove_rownames()
      child.list[[n]] <- child
    }
    child.list <- as.data.frame(do.call(rbind, child.list))
  },
  mc.cores = 10
  )
}

add.mutations <- function(new.proposal, mut.rate = 0.05) {
  ## Cell type mutation
  ct.mutation <- sample(x = c(0,1), size = nrow(new.proposal), prob = c(1-mut.rate, mut.rate), replace = TRUE) ## 0 = no mutation 1 = mutation
  new.proposal[which(ct.mutation ==1), "celltype"] <- sample(x = celltype, size = sum(ct.mutation), replace = TRUE) ## Add mutation
  ## Gene mutation
  gene.mutation <- sample(x = c(0,1), size = nrow(new.proposal), prob = c(1-mut.rate, mut.rate), replace = TRUE)
  
  mutated.loci <- new.proposal[which(gene.mutation == 1), "locus"]
  loci.gene.sample <- mclapply(X = as.list(mutated.loci), FUN = function(x) {
    new.gene <- loci2gene %>% dplyr::filter(cytoband == x) %>% slice_sample(n = 1) %>% dplyr::rename("gene"="gene_name")
  },
  mc.cores = 10)
  loci.gene.sample <- as.data.frame(do.call(rbind, loci.gene.sample))
  gene <- loci2gene %>% dplyr::filter(cytoband %in% mutated.loci)
  new.proposal[which(gene.mutation ==1), "gene"] <- loci.gene.sample$gene
  return(new.proposal)
}


## Jaccard similarity
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
