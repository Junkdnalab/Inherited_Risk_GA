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
                          size = n_proposal, ## Select 100 random proposals based on the prob below
                          prob = rep(rev(group.prob), each = nrow(proposal)/5))
  proposal <- proposal[rand.proposal,]
  }

mating.pair <- function(sel.proposal, n.pairs, n_proposal, df, khan.method, fitness.df, top.df, generation, n.loci, n.cores) {
  if(khan.method == TRUE) {
    ## Get khan info
    khan <- fitness.df[1,] ## Which gen and proposal is khan?
    khan.df <- top.df %>% dplyr::filter(n == khan$n, ## Get khan proposals
                                        gen == khan$gen) %>%
      dplyr::select(-gen) %>%
      mutate(n = n_proposal + generation)
    ############################
    ## Preparing proposal set ##
    ## Get row.num of random group 5 to kick out
    rand.kick <- sel.proposal %>% rownames_to_column(var = "row") %>% dplyr::filter(quantile == 5) %>% slice_sample(n = 1) %>% .$row
    proposal.set <- sel.proposal %>% dplyr::filter(rownames(.) != rand.kick)
    ## Add khan into proposal set
    khan <- khan %>% dplyr::rename("quantile"="gen") %>% mutate(n = n_proposal + generation)
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
  mc.cores = n.cores
  )
}

add.mutations <- function(new.proposal, mut.rate = 0.05, n.cores) {
  ## Cell type mutation
  ct.mutation <- sample(x = c(0,1), size = nrow(new.proposal), prob = c(1-mut.rate, mut.rate), replace = TRUE) ## 0 = no mutation 1 = mutation
  new.proposal[which(ct.mutation ==1), "celltype"] <- sample(x = celltype, size = sum(ct.mutation), replace = TRUE) ## Add mutation
  ## Gene mutation
  gene.mutation <- sample(x = c(0,1), size = nrow(new.proposal), prob = c(1-mut.rate, mut.rate), replace = TRUE)
  
  mutated.loci <- new.proposal[which(gene.mutation == 1), "locus"]
  loci.gene.sample <- mclapply(X = as.list(mutated.loci), FUN = function(x) {
    new.gene <- loci2gene %>% dplyr::filter(cytoband == x) %>% slice_sample(n = 1) %>% dplyr::rename("gene"="gene_name")
  },
  mc.cores = n.cores)
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

#########################
#########################
## genetic programming ##
#########################
#########################

gwas_gp <- function(init.proposal, generations, obj.fun, khan.method, selection.method, proposal.size, parent.size, n.cores) {
  ## Store result for fitness function
  fitness.res.list <- list() 
  ## Store initial proposal as df
  df <- init.proposal
  n.loci <- length(unique(df$locus))
  ## Set initial proposal as top proposal. This df will be updated to only include the top 1k proposals by avg of avg using current and new proposal
  top.df <- init.proposal
  #################################################################
  ## Loop through proposal for each generation using the obj fun ##
  for(generation in 1:generations) {
    message(paste0("Working on F", generation-1))
    
    ## When we're on the last generation, save the result of the proposal so we can output it. 
    if(generation == max(generations)) {
      last.gen <- df ## Record the last generation
      last.gen$gen <- generation -1
    }
    
    ## Storing results from obj fun for fitness function
    obj.fun.res <- list()
    
    #############################
    ## CEVM Objective Function ##
    if("cevm" %in% obj.fun) { ## If cevm then do the following
      cevm <- cevm.genes %>% mutate(cevm = 1)
  
      cevm.obj <- df %>% left_join(x = ., y = cevm, by = c("gene", "celltype")) %>%
        mutate(cevm = replace_na(cevm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cevm = sum(cevm) / length(cevm))
      obj.fun.res[["cevm"]] <- as.data.frame(cevm.obj)
    }
    
    ####################################
    ## FindMarkers Objective Function ##
    if("findmarkers" %in% obj.fun) {
      fm <- marker.hits %>% mutate(fm = 1)
        
      fm.obj <- df %>% left_join(x = ., y = fm, by = c("gene", "celltype")) %>%
        mutate(fm = replace_na(fm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(fm = sum(fm) / length(fm))
      obj.fun.res[["findmarkers"]] <- as.data.frame(fm.obj)
    }
    
    ########################################
    ## Combined gene marker - CEVM and FM ##
    if("marker.gene" %in% obj.fun) {
      ## marker gene data
      marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
      
      ## Data process
      marker.gene.obj <- df %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>% ## join marker.gene to proposal
        mutate(marker.gene = replace_na(marker.gene, replace = 0)) %>%  ## replace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(marker.gene = sum(marker.gene) / length(marker.gene)) ## get % pos hit
      obj.fun.res[["marker.gene"]] <- as.data.frame(marker.gene.obj) ## store in list
    }
    ##############################
    ## Magma Objective Function ##
    if("magma" %in% obj.fun) {
      magma.obj <- df %>% left_join(x = ., y = magma, by = "gene") %>% ## join magma to proposal
        mutate(magma = replace_na(magma, replace = 0)) %>%  ## replace NA with 0 
        dplyr::group_by(n) %>% dplyr::summarise(magma = sum(magma) / length(magma)) ## get % pos hit
      obj.fun.res[["magma"]] <- as.data.frame(magma.obj) ## store in list
    }
    
    ####################################
    ## Cancer Gene Objective Function ## 
    if("cancer.gene" %in% obj.fun) {
      hi.mut.obj <- df %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>% ## join cancer.gene to proposal
        mutate(cancer.gene = replace_na(cancer.gene, replace = 0)) %>%  ## rpelace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(cancer.gene = sum(cancer.gene) / length(cancer.gene)) ## get % pos hit
      obj.fun.res[["cancer.gene"]] <- as.data.frame(hi.mut.obj) ## store in list
    }
    
    ############################
    ## PPI Objective Function ##
    if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
      ppi.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- ppi %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                       relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x) 
        
      }, mc.cores = n.cores)
      ppi.obj <- as.data.frame(do.call(rbind, ppi.obj))
      
      #######################
      ## Intracellular PPI ##
      if("intra.ppi" %in% obj.fun) {
        intra.df <-  ppi.obj %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(intra.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        intra.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
          intra.ppi.genes <- intra.df %>% dplyr::filter(n == x)
          if(nrow(intra.ppi.genes) > 0) {
            intra.ppi.genes <- data.frame(gene = unique(c(intra.ppi.genes$gene1, intra.ppi.genes$gene2)),
                                          intra.ppi = 1)
            ## Get specific proposal set 
            intra.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = intra.ppi.genes, by = "gene") %>%
              mutate(intra.ppi = replace_na(intra.ppi, replace = 0)) %>%
              summarise(intra.ppi = sum(intra.ppi) /length(intra.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "intra.ppi")]
          } else {
            intra.res <- data.frame(n = x,
                                    intra.ppi = 0)
          }
        }, mc.cores = n.cores)
        intra.obj <- as.data.frame(do.call(rbind, intra.obj))
        obj.fun.res[["intra.ppi"]] <- intra.obj
      } ## End of intra ppi
      
      #############
      ## All PPI ##
      if("ppi" %in% obj.fun) {
        all.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          ppi.df <- ppi.obj %>% dplyr::filter(n == x)
          ppi.genes <- data.frame(gene = unique(c(ppi.df$gene1, ppi.df$gene2)),
                                  ppi = 1)
          ppi.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = ppi.genes, by = "gene") %>%
            mutate(ppi = replace_na(ppi, replace = 0)) %>%
            dplyr::summarise(ppi = sum(ppi) / length(ppi)) %>%
            mutate(n = x) %>% .[, c("n", "ppi")]
        }, mc.cores = n.cores)
        all.ppi.obj <- as.data.frame(do.call(rbind, all.ppi.obj))
        ppi.res <- data.frame(n = unique(df$n), 
                              ppi = 0)
        for(r in unique(all.ppi.obj$n)) {
          ppi.res[which(ppi.res$n == r), "ppi"] <- all.ppi.obj[which(all.ppi.obj$n == r), "ppi"]
        }
        obj.fun.res[["ppi"]] <- ppi.res
      } ## End of ppi
      
      ################
      ## Marker PPI ##
      if("marker.ppi" %in% obj.fun) {
        all.markers <- combined.gene.marker %>% mutate(marker.gene = 1)
        marker.obj <- ppi.obj %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() #%>% ## Keep positive hits.
        
        if(nrow(marker.obj) == 0) {
          marker.ppi.obj <- data.frame(n = unique(df$n),
                                       marker.ppi = 0)
        } else {
          marker.ppi.obj <- mclapply(X = as.list(unique(marker.obj$n)), FUN = function(x) {
            marker.df <- marker.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
            
            marker.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
              mutate(marker.ppi = replace_na(marker.ppi, replace = 0)) %>%
              dplyr::summarise(marker.ppi = sum(marker.ppi) / length(marker.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "marker.ppi")]
          }, mc.cores = n.cores)
          marker.ppi.obj <- as.data.frame(do.call(rbind, marker.ppi.obj))
          
          marker.ppi.obj <- data.frame(n = unique(df$n)) %>%
            left_join(x = ., y = marker.ppi.obj, by = "n") %>%
            mutate(marker.ppi = replace_na(marker.ppi, replace = 0))
        }
        
        obj.fun.res[["marker.ppi"]] <- marker.ppi.obj
      } ## End of marker.ppi
      
      #######################
      ## Intercellular PPI ##
      if("inter.ppi" %in% obj.fun) {
        marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        
        marker.obj <- ppi.obj %>% left_join(x = ., y = marker.gene, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = marker.gene, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() %>% ## Keep positive hits.
          dplyr::select(-marker.gene.x, -marker.gene.y) ## Get rid of columns no longer needed
        
        if(nrow(marker.obj) == 0) {
          inter.ppi.obj <- data.frame(n = unique(df$n),
                                      inter.ppi = 0)
        } else {
          ## First identity all source genes in celltalkdb. We will check if the source gene is found in 
          ## gene1 or gene2. If no source gene in either then remove. If source gene is found in both 
          ## Then remove. We will be filtering for rowsums == 1 on source_genesymbol1 and 2.
          inter.obj <- marker.obj %>% left_join(x = ., 
                                                y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                                by = c("gene1"="source_genesymbol"),
                                                relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("source_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                      by = c("gene2"="source_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>% 
            dplyr::rename("source_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
            
            ## Now trying to find targets the same way we did for source genes
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene1"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene2"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
            ## Removing target and source identity column
            dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
          
          inter.ppi.obj <- mclapply(X = as.list(unique(inter.obj$n)), FUN = function(x) {
            marker.df <- inter.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
            
            inter.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0)) %>%
              dplyr::summarise(inter.ppi = sum(inter.ppi) / length(inter.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "inter.ppi")]
          }, mc.cores = n.cores)
          inter.ppi.obj <- as.data.frame(do.call(rbind, inter.ppi.obj)) 
          
          if(nrow(inter.ppi.obj) == 0) {
            inter.ppi.obj <- data.frame(n = unique(df$n),
                                        inter.ppi = 0)
          } else {
            inter.ppi.obj <- data.frame(n = unique(df$n)) %>%
              left_join(x = ., y = inter.ppi.obj, by = "n") %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0))
          }
        }
        
        obj.fun.res[["inter.ppi"]] <- inter.ppi.obj
      } ## End of inter.ppi
      
      ###############
      ## Multi PPI ##
      if("multi.ppi" %in% obj.fun) {
        multi.obj <- ppi.obj %>% dplyr::select(-celltype1,-celltype2) %>% distinct()
        multi.res <- mclapply(X = as.list(unique(multi.obj$n)), FUN = function(x) {
          multi.df <- multi.obj %>% dplyr::filter(n %in% x)
          df1 <- multi.df %>% dplyr::select(gene1, locus1) %>% dplyr::rename("gene"="gene1", "locus"="locus1")
          df2 <- multi.df %>% dplyr::select(gene2, locus2) %>% dplyr::rename("gene"="gene2", "locus"="locus2")
          gene.locus.set <- rbind(df1, df2) %>% mutate(set = paste(gene, locus, sep = "_")) %>%  ## Combine df1 and df2 and create set name
            group_by(set) %>% ## Group by set
            tally() %>% ## ## Count up how many times the set appears
            dplyr::filter(n > 1) %>% nrow() ## Select gene locus combination with more than one interaction
          res <- data.frame(n = x, multi.ppi = gene.locus.set / n.loci)
        }, mc.cores = n.cores)
        multi.res <- as.data.frame(do.call(rbind, multi.res))
        obj.fun.res[["multi.ppi"]] <- multi.res
      } ## End of multi.ppi
      
      if("go.cc.ppi" %in% obj.fun) {
        go.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          go.ppi.df <- ppi.obj %>% dplyr::filter(n == x) %>% dplyr::select(-celltype1, -celltype2) %>%
            left_join(x = ., ## Add cellular component info for all genes in gene1 (NA means no info on gene)
                      y = map.go.gene[,c("external_gene_name", "name_1006")], 
                      by = c("gene1"="external_gene_name"), 
                      relationship = "many-to-many") %>%
            dplyr::rename("cc1"="name_1006") %>% ## Rename cc column as cc1 to link to gene1
            left_join(x = ., ## Add cc info for gene2
                      y = map.go.gene[,c("external_gene_name", "name_1006")],
                      by = c("gene2"="external_gene_name"),
                      relationship = "many-to-many") %>%
            dplyr::rename("cc2"="name_1006") %>% ## Rename cc column as cc2 to link to gene2
            .[which(.$cc1 == .$cc2),] ## Filter for cases where gene1 and gene2 are in same cellular compartment
          ## Get all genes from go.ppi.df
          go.genes <- data.frame(gene = unique(c(go.ppi.df$gene1, go.ppi.df$gene2)),
                                 go.cc = 1)
          go.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = go.genes, by = "gene") %>%
            mutate(go.cc = replace_na(go.cc, replace = 0)) %>%
            dplyr::summarise(go.cc = sum(go.cc) / length(go.cc)) %>%
            mutate(n = x) %>% .[, c("n", "go.cc")]
        }, mc.cores = n.cores)
        go.ppi.obj <- as.data.frame(do.call(rbind, go.ppi.obj))
        
        obj.fun.res[["go.cc.ppi"]] <- go.ppi.obj
      }
      
    } ## End of broad ppi umbrella
    
    #################################
    ## Promoter Objective Function ##
    if("promoter" %in% obj.fun) {
      promoter.obj <- df %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
        mutate(promoter = replace_na(promoter, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(promoter = sum(promoter)/length(promoter))
      
      obj.fun.res[["promoter"]] <- promoter.obj
    }
    
    #############################
    ## ATAC Objective Function ##
    ## Marker atac peak 
    if("marker.atac" %in% obj.fun) {
      m.atac.obj <- df %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
        mutate(marker.atac = replace_na(marker.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(marker.atac = sum(marker.atac)/length(marker.atac))
      obj.fun.res[["marker.atac"]] <- m.atac.obj
    }
    ## Common atac peak
    if("common.atac" %in% obj.fun) {
      c.atac.obj <- df %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
        mutate(common.atac = replace_na(common.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(common.atac = sum(common.atac)/length(common.atac))
      obj.fun.res[["common.atac"]] <- c.atac.obj
    }
    ########################
    ## lncrna interaction ##
    if("lncrna" %in% obj.fun) {
      lncrna.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- rna.protein.map %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                                   relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x)
      }, mc.cores = n.cores)
      lncrna.obj <- as.data.frame(do.call(rbind, lncrna.obj))
      
      all.lncrna.obj <- mclapply(X = as.list(unique(lncrna.obj$n)), FUN = function(x) {
        lncrna.df <- lncrna.obj %>% dplyr::filter(n == x)
        lncrna.genes <- data.frame(gene = unique(c(lncrna.df$gene1, lncrna.df$gene2)),
                                   lncrna = 1)
        lncrna.res <- df %>% dplyr::filter(n == x) %>%
          left_join(x = ., y = lncrna.genes, by = "gene") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0)) %>%
          dplyr::summarise(lncrna = sum(lncrna) / length(lncrna)) %>%
          mutate(n = x) %>% .[, c("n", "lncrna")]
      }, mc.cores = n.cores)
      all.lncrna.obj <- as.data.frame(do.call(rbind, all.lncrna.obj))
      
      if(nrow(all.lncrna.obj) == 0) { ## If no lncrna interaction then do this
        lncrna.res <- data.frame(n = unique(df$n),
                                 lncrna = 0)
      } else { ## If lncrna interaction present
        lncrna.res <- data.frame(n = unique(df$n)) %>%
          left_join(x = ., y = all.lncrna.obj, by = "n") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0))
      }
      obj.fun.res[["lncrna"]] <- lncrna.res
    } ## end of lncrna
    
    ######################
    ## tfbs interaction ##
    if("tfbs" %in% obj.fun) {

      tfbs.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n) %>%
          left_join(x = ., y = ct_exp, by = c("celltype"="type"), relationship = "many-to-many") %>% ## Add celltype and tf gene info
          left_join(x = ., y = mb.tfbs.disrupt, by = c("gene_name"="geneSymbol", "locus"="LEAD_SNP")) %>%
          dplyr::select(-gene_name) %>%
          dplyr::group_by(locus, celltype, gene) %>%
          mutate(tfbs = replace_na(tfbs, replace = 0)) %>% ## Add 0 where NA is
          dplyr::summarise(tfbs = sum(tfbs)) %>% ## sum all tfbs disrupt hits
          mutate(tfbs = case_when(tfbs > 0 ~ 1, ## boolean value conversion
                                  tfbs == 0 ~ 0),
                 n = x) %>%
          dplyr::select(celltype, locus, gene, n, tfbs) %>%
          dplyr::ungroup() 
        
      }, mc.cores = n.cores)
      tfbs.obj <- as.data.frame(do.call(rbind, tfbs.obj))
      tfbs.obj <- tfbs.obj %>% dplyr::group_by(n) %>% dplyr::summarise(tfbs = sum(tfbs)/length(tfbs))
      obj.fun.res[["tfbs"]] <- tfbs.obj
    } ## end of tfbs
    
    if("tfbs.marker" %in% obj.fun) {
      tfbs.marker.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with TF expressed data
        proposal.set <- ct_exp %>% mutate(tf.hit = 1) %>% 
          left_join(x = proposal.set, y = ., by = c("gene"="gene_name", "celltype"="type")) %>%
          mutate(tf.hit = replace_na(tf.hit, replace = 0))
        
        ## If there are any TF expressed in correct cell type then add in motifbreakR info
        if(any(proposal.set$tf.hit == 1)) {
          ## What tf.gene are a hit
          tf.genes <- proposal.set %>% dplyr::filter(tf.hit == 1) %>% pull(gene)
          ## What SNPs have motif broken?
          for(tf.gene in tf.genes) {
            tfbs.disrupt <- mb.tfbs.disrupt %>% dplyr::filter(geneSymbol %in% tf.gene)
            ## Now add a point for all lead SNPs with a MB disruption for that tf.gene
            for(tfbs in unique(tfbs.disrupt$LEAD_SNP)) {
              proposal.set[which(proposal.set$locus == tfbs), "tf.hit"] <- 1
            }
          }
        } ## end of IF condition 
        proposal.set <- proposal.set %>% mutate(n = x) %>% ## add back generation info
          dplyr::select(celltype, locus, gene, n, tf.hit) %>%
          dplyr::group_by(n) %>%
          dplyr::summarise(tfbs = sum(tf.hit)/length(tf.hit))
          
      }, mc.cores = n.cores)
      tfbs.marker.obj <- as.data.frame(do.call(rbind, tfbs.marker.obj))
    }
    
    #####################################
    ## Fitness mean objective function ##
    
    for(n in 1:length(obj.fun)) {
      if(n == 1) {
        fitness.obj <- obj.fun.res[[obj.fun[n]]]
      } else {
        add_test <- obj.fun.res[[obj.fun[n]]]
        fitness.obj <- left_join(x = fitness.obj, y = add_test, by = "n")
      }
    }
    fitness.obj[is.na(fitness.obj)] <- 0
    fitness.obj$fitness <- apply(X = fitness.obj %>% dplyr::select(-n), MARGIN = 1, FUN = function(x) {
      mean(x)
    })
    fitness.obj$gen <- generation - 1 
    
    ######################################
    ## Rank the proposals into 5 groups ##
    
    ranked.proposals <- fitness.obj %>% dplyr::select(n, fitness) %>% arrange(., desc(fitness)) %>%
      mutate(quantile = rep(x = 5:1, each = nrow(.)/5))
    ## Prob of quantiles
    group.prob <- c(0.02,0.03,0.05,0.1,0.8)
    names(group.prob) <- 1:max(ranked.proposals$quantile)
    
    #####################
    ## Updating top.df ##
    
    if(generation-1 == 0) { 
      top.df$gen <- generation - 1
      fitness.df <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
    } else {
      new.fitness <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
      ## Update top rank proposal df with new generation data. Basically, replace proposals that are weaker than the new proposals
      fitness.df <- rbind(fitness.df, new.fitness) %>% arrange(desc(fitness)) %>% .[1:1000,]
      ## Update top df proposal elements based on fitness.df
      keep.list <- fitness.df %>% dplyr::select(n, gen) %>% mutate(keep = 1)
      current.proposal <- df %>% mutate(gen = generation - 1)
      top.df <- rbind(top.df, current.proposal) %>% left_join(x = ., y = keep.list, by = c("n", "gen")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep)
    }
    
    ###################################
    ## Randomly Select 100 Proposals ##
    
    if("progressive" %in% selection.method) {
      
      ## Select 100 parents based on progressive scale
      sel.proposal <- progressive_quantile(group.prob = group.prob, n_proposal = parent.size, proposal = ranked.proposals)
    }
    
    if("lexicase" %in% selection.method) {
      sel.proposal <- list()
      
      for(lexicase.run in seq(parent.size)) { ## Find parents based on parent.size
        shuff.of <- sample(x = obj.fun, size = length(obj.fun), replace = FALSE) ## Shuffle obj func
        pool <- fitness.obj ## Store fitness.obj as pool 
        shuff.counter <- 1 ## Counter used for moving through the shuff.of
        while(nrow(pool) > 0 & shuff.counter <= length(shuff.of)) { ## While pool is greater than 0 and shuff.counter less than shuff.of
          message(paste0("Working on parent ", lexicase.run, ": Shuffle OF ", shuff.counter, " is ", shuff.of[[shuff.counter]]), ". Pool is ", nrow(pool))
          col.index <- which(colnames(pool) == shuff.of[[shuff.counter]])
          max.score <- pool %>% pull(shuff.of[[shuff.counter]]) %>% max() 
          pool <- pool %>% filter_at(col.index, all_vars(. == max.score)) 
          
          ## If only one parent then stop selection
          if(nrow(pool) == 1) { 
            sel.proposal[[lexicase.run]] <- pool %>% dplyr::select(n, gen)
            message(paste0("Parent ", lexicase.run, " found. Only one left. Stopped at shuff.counter ", shuff.counter))
            break
          } else if(nrow(pool) > 1) { ## If pool is still greater than 1
            if(shuff.counter == length(shuff.of)) {
              sel.proposal[[lexicase.run]] <- pool %>% .[sample(x = seq(nrow(pool)), size = 1),] %>% dplyr::select(n, gen)
              message(paste0("Parent ", lexicase.run, " found. More than one. Random select Stopped at shuff.counter ", shuff.counter))
              break
            } else {
              shuff.counter <- shuff.counter + 1
            }
          } 
        } ## end of while loop
      } ## end of lexicase selection
      sel.proposal <- as.data.frame(do.call(rbind, sel.proposal)) ## Process parent selection result
      sel.proposal <- sel.proposal %>% left_join(x = ., y = ranked.proposals, by = "n") %>% dplyr::select(-gen)
    } ## end of lexicase function
    
    
    
    ############
    ## Mating ##
    new.proposal <- mating.pair(sel.proposal = sel.proposal, n_proposal = proposal.size, n.pairs = parent.size, df = df, khan.method = khan.method, fitness.df = fitness.df, top.df = top.df,
                                generation = generation, n.loci = n.loci, n.cores = n.cores)
    new.proposal <- as.data.frame(do.call(rbind,new.proposal)) %>% remove_rownames() %>%
      mutate(n = rep(1:proposal.size, each = n.loci))
    ##############
    ## Mutation ##
    new.proposal <- add.mutations(new.proposal = new.proposal, n.cores = n.cores)
    
    ## Insert khan to the proposal
    rand.kick <- new.proposal %>% rownames_to_column(var = "row") %>% slice_sample(n = 1) %>% .$n
    new.proposal <- new.proposal %>% dplyr::filter(n != rand.kick)
    
    khan <- fitness.df[1,] ## Which gen and proposal is khan?
    khan.df <- top.df %>% dplyr::filter(n == khan$n, ## Get khan proposals
                                        gen == khan$gen) %>%
      dplyr::select(-gen) %>%
      mutate(n = proposal.size + generation)
    
    new.proposal <- rbind(khan.df, new.proposal)
    
    ###########################################
    ## Remove one random and add khan inside ##
    
    df <- new.proposal ## Overwrite multi_proposal with new.proposal
    
    ###############
    ## Store res ##
    fitness.res.list[[generation]] <- fitness.obj
    
    
  } ## End loop for generations
  fitness.res.list[["proposal"]] <- last.gen
  fitness.res.list[["top"]] <- top.df
  fitness.res.list[["fitness"]] <-fitness.df
  fitness.res.list[["selection"]] <- selection.method
  return(fitness.res.list)
} ## End of wrapper

#######################
## NO OMICS GA MODEL ##
#######################

gwas_gp_no_omics <- function(init.proposal, generations, obj.fun, exclude.celltype = c("ERGpos_Tumor", "ERGneg_Tumor"), khan.method) {
  ## Store result for fitness function
  fitness.res.list <- list() 
  ## Store initial proposal as df
  df <- init.proposal
  n.loci <- length(unique(df$locus))
  ## Set initial proposal as top proposal. This df will be updated to only include the top 1k proposals by avg of avg using current and new proposal
  top.df <- init.proposal
  #################################################################
  ## Loop through proposal for each generation using the obj fun ##
  for(generation in 1:generations) {
    message(paste0("Working on F", generation-1))
    
    ## When we're on the last generation, save the result of the proposal so we can output it. 
    if(generation == max(generations)) {
      last.gen <- df ## Record the last generation
      last.gen$gen <- generation -1
    }
    
    ## Storing results from obj fun for fitness function
    obj.fun.res <- list()
    
    #############################
    ## CEVM Objective Function ##
    if("cevm" %in% obj.fun) { ## If cevm then do the following
      if(is.null(exclude.celltype)) { 
        cevm <- cevm.genes %>% dplyr::filter(!celltype %in% exclude.celltype) %>% mutate(cevm = 1)
      } else {
        cevm <- cevm.genes %>% mutate(cevm = 1)
      }
      cevm.obj <- df %>% left_join(x = ., y = cevm, by = c("gene", "celltype")) %>%
        mutate(cevm = replace_na(cevm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cevm = sum(cevm) / length(cevm))
      obj.fun.res[["cevm"]] <- as.data.frame(cevm.obj)
    }
    
    ####################################
    ## FindMarkers Objective Function ##
    if("findmarkers" %in% obj.fun) {
      if(is.null(exclude.celltype)) {
        fm <- marker.hits %>% dplyr::filter(!celltype %in% exclude.celltype) %>% mutate(fm = 1)
      } else {
        fm <- marker.hits %>% mutate(fm = 1)
      }
      fm.obj <- df %>% left_join(x = ., y = fm, by = c("gene", "celltype")) %>%
        mutate(fm = replace_na(fm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(fm = sum(fm) / length(fm))
      obj.fun.res[["findmarkers"]] <- as.data.frame(fm.obj)
    }
    
    ########################################
    ## Combined gene marker - CEVM and FM ##
    if("marker.gene" %in% obj.fun) {
      if(is.null(exclude.celltype)) {
        marker.gene <- combined.gene.marker %>% dplyr::filter(!celltype %in% exclude.celltype) %>% mutate(marker.gene = 1)
      } else {
        marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
      }
      marker.gene.obj <- df %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>%
        mutate(marker.gene = replace_na(marker.gene, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(marker.gene = sum(marker.gene) / length(marker.gene))
      obj.fun.res[["marker.gene"]] <- as.data.frame(marker.gene.obj)
    }
    ##############################
    ## Magma Objective Function ##
    if("magma" %in% obj.fun) {
      magma.obj <- df %>% left_join(x = ., y = magma, by = "gene") %>%
        mutate(magma = replace_na(magma, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(magma = sum(magma) / length(magma))
      obj.fun.res[["magma"]] <- as.data.frame(magma.obj)
    }
    
    if("cancer.gene" %in% obj.fun) {
      hi.mut.obj <- df %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>%
        mutate(cancer.gene = replace_na(cancer.gene, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cancer.gene = sum(cancer.gene) / length(cancer.gene))
      obj.fun.res[["cancer.gene"]] <- as.data.frame(hi.mut.obj)
    }
    
    ############################
    ## PPI Objective Function ##
    if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
      ppi.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- ppi %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                       relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x) 
        
      }, mc.cores = 10)
      ppi.obj <- as.data.frame(do.call(rbind, ppi.obj))
      
      if("intra.ppi" %in% obj.fun) {
        intra.df <-  ppi.obj %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(intra.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        intra.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
          intra.ppi.genes <- intra.df %>% dplyr::filter(n == x)
          if(nrow(intra.ppi.genes) > 0) {
            intra.ppi.genes <- data.frame(gene = unique(c(intra.ppi.genes$gene1, intra.ppi.genes$gene2)),
                                          intra.ppi = 1)
            ## Get specific proposal set 
            intra.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = intra.ppi.genes, by = "gene") %>%
              mutate(intra.ppi = replace_na(intra.ppi, replace = 0)) %>%
              summarise(intra.ppi = sum(intra.ppi) /length(intra.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "intra.ppi")]
          } else {
            intra.res <- data.frame(n = x,
                                    intra.ppi = 0)
          }
        }, mc.cores = 10)
        intra.obj <- as.data.frame(do.call(rbind, intra.obj))
        obj.fun.res[["intra.ppi"]] <- intra.obj
      } ## End of intra ppi
      
      if("ppi" %in% obj.fun) {
        all.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          ppi.df <- ppi.obj %>% dplyr::filter(n == x)
          ppi.genes <- data.frame(gene = unique(c(ppi.df$gene1, ppi.df$gene2)),
                                  ppi = 1)
          ppi.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = ppi.genes, by = "gene") %>%
            mutate(ppi = replace_na(ppi, replace = 0)) %>%
            dplyr::summarise(ppi = sum(ppi) / length(ppi)) %>%
            mutate(n = x) %>% .[, c("n", "ppi")]
        }, mc.cores = 10)
        all.ppi.obj <- as.data.frame(do.call(rbind, all.ppi.obj))
        ppi.res <- data.frame(n = unique(df$n), 
                              ppi = 0)
        for(r in unique(all.ppi.obj$n)) {
          ppi.res[which(ppi.res$n == r), "ppi"] <- all.ppi.obj[which(all.ppi.obj$n == r), "ppi"]
        }
        obj.fun.res[["ppi"]] <- ppi.res
      } ## End of ppi
      
      if("marker.ppi" %in% obj.fun) {
        if(is.null(exclude.celltype)) {
          # fm <- marker.hits %>% dplyr::filter(!celltype %in% exclude.celltype) %>% mutate(marker = 1)
          # cevm <- cevm.genes %>% dplyr::filter(!celltype %in% exclude.celltype) %>% mutate(marker = 1)
          marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        } else {
          marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        }
        all.markers <- marker.gene %>% distinct() #rbind(fm, cevm) %>% distinct()
        marker.obj <- ppi.obj %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() #%>% ## Keep positive hits.
        
        if(nrow(marker.obj) == 0) {
          marker.ppi.obj <- data.frame(n = unique(df$n),
                                       marker.ppi = 0)
        } else {
          marker.ppi.obj <- mclapply(X = as.list(unique(marker.obj$n)), FUN = function(x) {
            marker.df <- marker.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
            
            marker.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
              mutate(marker.ppi = replace_na(marker.ppi, replace = 0)) %>%
              dplyr::summarise(marker.ppi = sum(marker.ppi) / length(marker.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "marker.ppi")]
          }, mc.cores = 10)
          marker.ppi.obj <- as.data.frame(do.call(rbind, marker.ppi.obj))
        }
        
        obj.fun.res[["marker.ppi"]] <- marker.ppi.obj
      } ## End of marker.ppi
      
      if("inter.ppi" %in% obj.fun) {
        inter.df <-  ppi.obj %>% mutate(inter.commun = case_when(celltype1 != celltype2 ~ 1,
                                                                 celltype1 == celltype2 ~ 0)) %>% 
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(inter.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        
        if(nrow(inter.df) == 0) {
          inter.ppi.obj <- data.frame(n = unique(df$n),
                                      inter.ppi = 0)
        } else {
          ## First identity all source genes in celltalkdb. We will check if the source gene is found in 
          ## gene1 or gene2. If no source gene in either then remove. If source gene is found in both 
          ## Then remove. We will be filtering for rowsums == 1 on source_genesymbol1 and 2.
          inter.obj <- inter.df %>% left_join(x = ., 
                                                y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                                by = c("gene1"="source_genesymbol"),
                                                relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("source_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                      by = c("gene2"="source_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>% 
            dplyr::rename("source_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
            
            ## Now trying to find targets the same way we did for source genes
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene1"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene2"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
            ## Removing target and source identity column
            dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
          
          inter.ppi.obj <- mclapply(X = as.list(unique(inter.obj$n)), FUN = function(x) {
            marker.df <- inter.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
            
            inter.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0)) %>%
              dplyr::summarise(inter.ppi = sum(inter.ppi) / length(inter.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "inter.ppi")]
          }, mc.cores = 10)
          inter.ppi.obj <- as.data.frame(do.call(rbind, inter.ppi.obj)) 
          
          if(nrow(inter.ppi.obj) == 0) {
            inter.ppi.obj <- data.frame(n = unique(df$n),
                                        inter.ppi = 0)
          }
        }
        
        obj.fun.res[["inter.ppi"]] <- inter.ppi.obj
      } ## End of inter.ppi
      
      if("multi.ppi" %in% obj.fun) {
        multi.obj <- ppi.obj %>% dplyr::select(-celltype1,-celltype2) %>% distinct()
        multi.res <- mclapply(X = as.list(unique(multi.obj$n)), FUN = function(x) {
          multi.df <- multi.obj %>% dplyr::filter(n %in% x)
          df1 <- multi.df %>% dplyr::select(gene1, locus1) %>% dplyr::rename("gene"="gene1", "locus"="locus1")
          df2 <- multi.df %>% dplyr::select(gene2, locus2) %>% dplyr::rename("gene"="gene2", "locus"="locus2")
          gene.locus.set <- rbind(df1, df2) %>% mutate(set = paste(gene, locus, sep = "_")) %>%  ## Combine df1 and df2 and create set name
            group_by(set) %>% ## Group by set
            tally() %>% ## ## Count up how many times the set appears
            dplyr::filter(n > 1) %>% nrow() ## Select gene locus combination with more than one interaction
          res <- data.frame(n = x, multi.ppi = gene.locus.set / n.loci)
        }, mc.cores = 10)
        multi.res <- as.data.frame(do.call(rbind, multi.res))
        obj.fun.res[["multi.ppi"]] <- multi.res
      } ## End of multi.ppi
      
      if("go.cc.ppi" %in% obj.fun) {
        go.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          go.ppi.df <- ppi.obj %>% dplyr::filter(n == x) %>% dplyr::select(-celltype1, -celltype2) %>%
            left_join(x = ., ## Add cellular component info for all genes in gene1 (NA means no info on gene)
                      y = map.go.gene[,c("external_gene_name", "name_1006")], 
                      by = c("gene1"="external_gene_name"), 
                      relationship = "many-to-many") %>%
            dplyr::rename("cc1"="name_1006") %>% ## Rename cc column as cc1 to link to gene1
            left_join(x = ., ## Add cc info for gene2
                      y = map.go.gene[,c("external_gene_name", "name_1006")],
                      by = c("gene2"="external_gene_name"),
                      relationship = "many-to-many") %>%
            dplyr::rename("cc2"="name_1006") %>% ## Rename cc column as cc2 to link to gene2
            .[which(.$cc1 == .$cc2),] ## Filter for cases where gene1 and gene2 are in same cellular compartment
          ## Get all genes from go.ppi.df
          go.genes <- data.frame(gene = unique(c(go.ppi.df$gene1, go.ppi.df$gene2)),
                                 go.cc = 1)
          go.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = go.genes, by = "gene") %>%
            mutate(go.cc = replace_na(go.cc, replace = 0)) %>%
            dplyr::summarise(go.cc = sum(go.cc) / length(go.cc)) %>%
            mutate(n = x) %>% .[, c("n", "go.cc")]
        }, mc.cores = 10)
        go.ppi.obj <- as.data.frame(do.call(rbind, go.ppi.obj))
        
        obj.fun.res[["go.cc.ppi"]] <- go.ppi.obj
      }
      
    } ## End of broad ppi umbrella
    
    #################################
    ## Promoter Objective Function ##
    if("promoter" %in% obj.fun) {
      promoter.obj <- df %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
        mutate(promoter = replace_na(promoter, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(promoter = sum(promoter)/length(promoter))
      
      obj.fun.res[["promoter"]] <- promoter.obj
    }
    
    #############################
    ## ATAC Objective Function ##
    ## Marker atac peak 
    if("marker.atac" %in% obj.fun) {
      m.atac.obj <- df %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
        mutate(marker.atac = replace_na(marker.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(marker.atac = sum(marker.atac)/length(marker.atac))
      obj.fun.res[["marker.atac"]] <- m.atac.obj
    }
    ## Common atac peak
    if("common.atac" %in% obj.fun) {
      c.atac.obj <- df %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
        mutate(common.atac = replace_na(common.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(common.atac = sum(common.atac)/length(common.atac))
      obj.fun.res[["common.atac"]] <- c.atac.obj
    }
    ########################
    ## lncrna interaction ##
    if("lncrna.ppi" %in% obj.fun) {
      lncrna.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- rna.protein.map %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                                   relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x)
      }, mc.cores = 10)
      lncrna.obj <- as.data.frame(do.call(rbind, lncrna.obj))
      
      all.lncrna.obj <- mclapply(X = as.list(unique(lncrna.obj$n)), FUN = function(x) {
        lncrna.df <- lncrna.obj %>% dplyr::filter(n == x)
        lncrna.genes <- data.frame(gene = unique(c(lncrna.df$gene1, lncrna.df$gene2)),
                                   lncrna = 1)
        lncrna.res <- df %>% dplyr::filter(n == x) %>%
          left_join(x = ., y = lncrna.genes, by = "gene") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0)) %>%
          dplyr::summarise(lncrna = sum(lncrna) / length(lncrna)) %>%
          mutate(n = x) %>% .[, c("n", "lncrna")]
      }, mc.cores = 10)
      all.lncrna.obj <- as.data.frame(do.call(rbind, all.lncrna.obj))
      
      lncrna.res <- data.frame(n = unique(df$n)) %>%
        left_join(x = ., y = all.lncrna.obj, by = "n") %>%
        mutate(lncrna = replace_na(lncrna, replace = 0))
      obj.fun.res[["lncrna.ppi"]] <- lncrna.res
    } ## end of lncrna
    
    ########################
    ## lncrna interaction ##
    if("tfbs" %in% obj.fun) {
      
      tfbs.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n) %>%
          left_join(x = ., y = ct_exp, by = c("celltype"="type"), relationship = "many-to-many") %>% ## Add celltype and tf gene info
          left_join(x = ., y = mb.tfbs.disrupt, by = c("gene_name"="geneSymbol", "locus"="LEAD_SNP")) %>%
          dplyr::select(-gene_name) %>%
          dplyr::group_by(locus, celltype, gene) %>%
          mutate(tfbs = replace_na(tfbs, replace = 0)) %>% ## Add 0 where NA is
          dplyr::summarise(tfbs = sum(tfbs)) %>% ## sum all tfbs disrupt hits
          mutate(tfbs = case_when(tfbs > 0 ~ 1, ## boolean value conversion
                                  tfbs == 0 ~ 0),
                 n = x) %>%
          dplyr::select(celltype, locus, gene, n, tfbs) %>%
          dplyr::ungroup() 
        
      }, mc.cores = 10)
      tfbs.obj <- as.data.frame(do.call(rbind, tfbs.obj))
      tfbs.obj <- tfbs.obj %>% dplyr::group_by(n) %>% dplyr::summarise(tfbs = sum(tfbs)/length(tfbs))
      obj.fun.res[["tfbs"]] <- tfbs.obj
    } ## end of tfbs
    
    
    #####################################
    ## Fitness mean objective function ##
    
    for(n in 1:length(obj.fun)) {
      if(n == 1) {
        fitness.obj <- obj.fun.res[[obj.fun[n]]]
      } else {
        add_test <- obj.fun.res[[obj.fun[n]]]
        fitness.obj <- left_join(x = fitness.obj, y = add_test, by = "n")
      }
    }
    fitness.obj[is.na(fitness.obj)] <- 0
    fitness.obj$fitness <- apply(X = fitness.obj %>% dplyr::select(-n), MARGIN = 1, FUN = function(x) {
      mean(x)
    })
    fitness.obj$gen <- generation - 1 
    
    
    ##########
    ## Rank ##
    ## Rank the proposals into 5 groups
    ranked.proposals <- fitness.obj %>% dplyr::select(n, fitness) %>% arrange(., desc(fitness)) %>%
      mutate(quantile = rep(x = 5:1, each = nrow(.)/5))
    ## Prob of quantiles
    group.prob <- c(0.02,0.03,0.05,0.1,0.8)
    names(group.prob) <- 1:max(ranked.proposals$quantile)
    
    #####################
    ## Updating top.df ##
    
    if(generation-1 == 0) { 
      top.df$gen <- generation - 1
      fitness.df <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
    } else {
      new.fitness <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
      ## Update top rank proposal df with new generation data. Basically, replace proposals that are weaker than the new proposals
      fitness.df <- rbind(fitness.df, new.fitness) %>% arrange(desc(fitness)) %>% .[1:1000,]
      ## Update top df proposal elements based on fitness.df
      keep.list <- fitness.df %>% dplyr::select(n, gen) %>% mutate(keep = 1)
      current.proposal <- df %>% mutate(gen = generation - 1)
      top.df <- rbind(top.df, current.proposal) %>% left_join(x = ., y = keep.list, by = c("n", "gen")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep)
    }
    
    ###################################
    ## Randomly Select 100 Proposals ##
    sel.proposal <- progressive_quantile(group.prob = group.prob, n_proposal = 1e2, proposal = ranked.proposals)
    
    ############
    ## Mating ##
    new.proposal <- mating.pair(sel.proposal = sel.proposal, df = df, khan.method = khan.method, fitness.df = fitness.df, top.df = top.df,
                                generation = generation, n.loci = n.loci)
    new.proposal <- as.data.frame(do.call(rbind,new.proposal)) %>% remove_rownames() %>%
      mutate(n = rep(1:1000, each = n.loci))
    ##############
    ## Mutation ##
    new.proposal <- add.mutations(new.proposal = new.proposal)
    
    ## Insert khan to the proposal
    rand.kick <- new.proposal %>% rownames_to_column(var = "row") %>% slice_sample(n = 1) %>% .$n
    new.proposal <- new.proposal %>% dplyr::filter(n != rand.kick)
    
    khan <- fitness.df[1,] ## Which gen and proposal is khan?
    khan.df <- top.df %>% dplyr::filter(n == khan$n, ## Get khan proposals
                                        gen == khan$gen) %>%
      dplyr::select(-gen) %>%
      mutate(n = 1000 + generation)
    
    new.proposal <- rbind(khan.df, new.proposal)
    
    ###########################################
    ## Remove one random and add khan inside ##
    
    df <- new.proposal ## Overwrite multi_proposal with new.proposal
    
    ###############
    ## Store res ##
    fitness.res.list[[generation]] <- fitness.obj
    
    
  } ## End loop for generations
  fitness.res.list[["proposal"]] <- last.gen
  fitness.res.list[["top"]] <- top.df
  fitness.res.list[["fitness"]] <-fitness.df
  return(fitness.res.list)
} ## End of wrapper

##############################
## GA model with OF scaling ##
##############################

gwas_gp_scaled <- function(init.proposal, generations, obj.fun, khan.method, selection.method, proposal.size, parent.size, n.cores) {
  ## Store result for fitness function
  fitness.res.list <- list() 
  ## Store initial proposal as df
  df <- init.proposal
  n.loci <- length(unique(df$locus))
  ## Set initial proposal as top proposal. This df will be updated to only include the top 1k proposals by avg of avg using current and new proposal
  top.df <- init.proposal
  #################################################################
  ## Loop through proposal for each generation using the obj fun ##
  for(generation in 1:generations) {
    message(paste0("Working on F", generation-1))
    
    ## When we're on the last generation, save the result of the proposal so we can output it. 
    if(generation == max(generations)) {
      last.gen <- df ## Record the last generation
      last.gen$gen <- generation -1
    }
    
    ## Storing results from obj fun for fitness function
    obj.fun.res <- list()
    
    #############################
    ## CEVM Objective Function ##
    if("cevm" %in% obj.fun) { ## If cevm then do the following
      cevm <- cevm.genes %>% mutate(cevm = 1)
      
      cevm.obj <- df %>% left_join(x = ., y = cevm, by = c("gene", "celltype")) %>%
        mutate(cevm = replace_na(cevm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cevm = sum(cevm) / as.numeric(max.hit[which(max.hit$of == "cevm"), "max_loci"]))
      obj.fun.res[["cevm"]] <- as.data.frame(cevm.obj)
    }
    
    ####################################
    ## FindMarkers Objective Function ##
    if("findmarkers" %in% obj.fun) {
      fm <- marker.hits %>% mutate(fm = 1)
      
      fm.obj <- df %>% left_join(x = ., y = fm, by = c("gene", "celltype")) %>%
        mutate(fm = replace_na(fm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(fm = sum(fm) / as.numeric(max.hit[which(max.hit$of == "fm"), "max_loci"]))
      obj.fun.res[["findmarkers"]] <- as.data.frame(fm.obj)
    }
    
    ########################################
    ## Combined gene marker - CEVM and FM ##
    if("marker.gene" %in% obj.fun) {
      ## marker gene data
      marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
      
      ## Data process
      marker.gene.obj <- df %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>% ## join marker.gene to proposal
        mutate(marker.gene = replace_na(marker.gene, replace = 0)) %>%  ## replace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(marker.gene = sum(marker.gene) / as.numeric(max.hit[which(max.hit$of == "marker.gene"), "max_loci"])) ## get % pos hit
      obj.fun.res[["marker.gene"]] <- as.data.frame(marker.gene.obj) ## store in list
    }
    ##############################
    ## Magma Objective Function ##
    if("magma" %in% obj.fun) {
      magma.obj <- df %>% left_join(x = ., y = magma, by = "gene") %>% ## join magma to proposal
        mutate(magma = replace_na(magma, replace = 0)) %>%  ## replace NA with 0 
        dplyr::group_by(n) %>% dplyr::summarise(magma = sum(magma) / as.numeric(max.hit[which(max.hit$of == "magma"), "max_loci"])) ## get % pos hit
      obj.fun.res[["magma"]] <- as.data.frame(magma.obj) ## store in list
    }
    
    ####################################
    ## Cancer Gene Objective Function ## 
    if("cancer.gene" %in% obj.fun) {
      hi.mut.obj <- df %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>% ## join cancer.gene to proposal
        mutate(cancer.gene = replace_na(cancer.gene, replace = 0)) %>%  ## rpelace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(cancer.gene = sum(cancer.gene) / as.numeric(max.hit[which(max.hit$of == "cancer.gene"), "max_loci"])) ## get % pos hit
      obj.fun.res[["cancer.gene"]] <- as.data.frame(hi.mut.obj) ## store in list
    }
    
    ############################
    ## PPI Objective Function ##
    if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
      ppi.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- ppi %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                       relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x) 
        
      }, mc.cores = n.cores)
      ppi.obj <- as.data.frame(do.call(rbind, ppi.obj))
      
      #######################
      ## Intracellular PPI ##
      if("intra.ppi" %in% obj.fun) {
        intra.df <-  ppi.obj %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(intra.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        intra.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
          intra.ppi.genes <- intra.df %>% dplyr::filter(n == x)
          if(nrow(intra.ppi.genes) > 0) {
            intra.ppi.genes <- data.frame(gene = unique(c(intra.ppi.genes$gene1, intra.ppi.genes$gene2)),
                                          intra.ppi = 1)
            ## Get specific proposal set 
            intra.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = intra.ppi.genes, by = "gene") %>%
              mutate(intra.ppi = replace_na(intra.ppi, replace = 0)) %>%
              summarise(intra.ppi = sum(intra.ppi) /as.numeric(max.hit[which(max.hit$of == "intra.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "intra.ppi")]
          } else {
            intra.res <- data.frame(n = x,
                                    intra.ppi = 0)
          }
        }, mc.cores = n.cores)
        intra.obj <- as.data.frame(do.call(rbind, intra.obj))
        obj.fun.res[["intra.ppi"]] <- intra.obj
      } ## End of intra ppi
      
      #############
      ## All PPI ##
      if("ppi" %in% obj.fun) {
        all.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          ppi.df <- ppi.obj %>% dplyr::filter(n == x)
          ppi.genes <- data.frame(gene = unique(c(ppi.df$gene1, ppi.df$gene2)),
                                  ppi = 1)
          ppi.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = ppi.genes, by = "gene") %>%
            mutate(ppi = replace_na(ppi, replace = 0)) %>%
            dplyr::summarise(ppi = sum(ppi) / as.numeric(max.hit[which(max.hit$of == "ppi"), "max_loci"])) %>%
            mutate(n = x) %>% .[, c("n", "ppi")]
        }, mc.cores = n.cores)
        all.ppi.obj <- as.data.frame(do.call(rbind, all.ppi.obj))
        ppi.res <- data.frame(n = unique(df$n), 
                              ppi = 0)
        for(r in unique(all.ppi.obj$n)) {
          ppi.res[which(ppi.res$n == r), "ppi"] <- all.ppi.obj[which(all.ppi.obj$n == r), "ppi"]
        }
        obj.fun.res[["ppi"]] <- ppi.res
      } ## End of ppi
      
      ################
      ## Marker PPI ##
      if("marker.ppi" %in% obj.fun) {
        all.markers <- combined.gene.marker %>% mutate(marker.gene = 1)
        marker.obj <- ppi.obj %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() #%>% ## Keep positive hits.
        
        if(nrow(marker.obj) == 0) {
          marker.ppi.obj <- data.frame(n = unique(df$n),
                                       marker.ppi = 0)
        } else {
          marker.ppi.obj <- mclapply(X = as.list(unique(marker.obj$n)), FUN = function(x) {
            marker.df <- marker.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
            
            marker.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
              mutate(marker.ppi = replace_na(marker.ppi, replace = 0)) %>%
              dplyr::summarise(marker.ppi = sum(marker.ppi) / as.numeric(max.hit[which(max.hit$of == "marker.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "marker.ppi")]
          }, mc.cores = n.cores)
          marker.ppi.obj <- as.data.frame(do.call(rbind, marker.ppi.obj))
          
          marker.ppi.obj <- data.frame(n = unique(df$n)) %>%
            left_join(x = ., y = marker.ppi.obj, by = "n") %>%
            mutate(marker.ppi = replace_na(marker.ppi, replace = 0))
        }
        
        obj.fun.res[["marker.ppi"]] <- marker.ppi.obj
      } ## End of marker.ppi
      
      #######################
      ## Intercellular PPI ##
      if("inter.ppi" %in% obj.fun) {
        marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        
        marker.obj <- ppi.obj %>% left_join(x = ., y = marker.gene, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = marker.gene, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() %>% ## Keep positive hits.
          dplyr::select(-marker.gene.x, -marker.gene.y) ## Get rid of columns no longer needed
        
        if(nrow(marker.obj) == 0) {
          inter.ppi.obj <- data.frame(n = unique(df$n),
                                      inter.ppi = 0)
        } else {
          ## First identity all source genes in celltalkdb. We will check if the source gene is found in 
          ## gene1 or gene2. If no source gene in either then remove. If source gene is found in both 
          ## Then remove. We will be filtering for rowsums == 1 on source_genesymbol1 and 2.
          inter.obj <- marker.obj %>% left_join(x = ., 
                                                y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                                by = c("gene1"="source_genesymbol"),
                                                relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("source_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                      by = c("gene2"="source_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>% 
            dplyr::rename("source_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
            
            ## Now trying to find targets the same way we did for source genes
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene1"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene2"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
            ## Removing target and source identity column
            dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
          
          inter.ppi.obj <- mclapply(X = as.list(unique(inter.obj$n)), FUN = function(x) {
            marker.df <- inter.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
            
            inter.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0)) %>%
              dplyr::summarise(inter.ppi = sum(inter.ppi) / as.numeric(max.hit[which(max.hit$of == "inter.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "inter.ppi")]
          }, mc.cores = n.cores)
          inter.ppi.obj <- as.data.frame(do.call(rbind, inter.ppi.obj)) 
          
          if(nrow(inter.ppi.obj) == 0) {
            inter.ppi.obj <- data.frame(n = unique(df$n),
                                        inter.ppi = 0)
          } else {
            inter.ppi.obj <- data.frame(n = unique(df$n)) %>%
              left_join(x = ., y = inter.ppi.obj, by = "n") %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0))
          }
        }
        
        obj.fun.res[["inter.ppi"]] <- inter.ppi.obj
      } ## End of inter.ppi
      
      ###############
      ## Multi PPI ##
      if("multi.ppi" %in% obj.fun) {
        multi.obj <- ppi.obj %>% dplyr::select(-celltype1,-celltype2) %>% distinct()
        multi.res <- mclapply(X = as.list(unique(multi.obj$n)), FUN = function(x) {
          multi.df <- multi.obj %>% dplyr::filter(n %in% x)
          df1 <- multi.df %>% dplyr::select(gene1, locus1) %>% dplyr::rename("gene"="gene1", "locus"="locus1")
          df2 <- multi.df %>% dplyr::select(gene2, locus2) %>% dplyr::rename("gene"="gene2", "locus"="locus2")
          gene.locus.set <- rbind(df1, df2) %>% mutate(set = paste(gene, locus, sep = "_")) %>%  ## Combine df1 and df2 and create set name
            group_by(set) %>% ## Group by set
            tally() %>% ## ## Count up how many times the set appears
            dplyr::filter(n > 1) %>% nrow() ## Select gene locus combination with more than one interaction
          res <- data.frame(n = x, multi.ppi = gene.locus.set / n.loci)
        }, mc.cores = 60)
        multi.res <- as.data.frame(do.call(rbind, multi.res))
        obj.fun.res[["multi.ppi"]] <- multi.res
      } ## End of multi.ppi
      
      if("go.cc.ppi" %in% obj.fun) {
        go.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          go.ppi.df <- ppi.obj %>% dplyr::filter(n == x) %>% dplyr::select(-celltype1, -celltype2) %>%
            left_join(x = ., ## Add cellular component info for all genes in gene1 (NA means no info on gene)
                      y = map.go.gene[,c("external_gene_name", "name_1006")], 
                      by = c("gene1"="external_gene_name"), 
                      relationship = "many-to-many") %>%
            dplyr::rename("cc1"="name_1006") %>% ## Rename cc column as cc1 to link to gene1
            left_join(x = ., ## Add cc info for gene2
                      y = map.go.gene[,c("external_gene_name", "name_1006")],
                      by = c("gene2"="external_gene_name"),
                      relationship = "many-to-many") %>%
            dplyr::rename("cc2"="name_1006") %>% ## Rename cc column as cc2 to link to gene2
            .[which(.$cc1 == .$cc2),] ## Filter for cases where gene1 and gene2 are in same cellular compartment
          ## Get all genes from go.ppi.df
          go.genes <- data.frame(gene = unique(c(go.ppi.df$gene1, go.ppi.df$gene2)),
                                 go.cc = 1)
          go.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = go.genes, by = "gene") %>%
            mutate(go.cc = replace_na(go.cc, replace = 0)) %>%
            dplyr::summarise(go.cc = sum(go.cc) / as.numeric(max.hit[which(max.hit$of == "go.cc"), "max_loci"])) %>%
            mutate(n = x) %>% .[, c("n", "go.cc")]
        }, mc.cores = n.cores)
        go.ppi.obj <- as.data.frame(do.call(rbind, go.ppi.obj))
        
        obj.fun.res[["go.cc.ppi"]] <- go.ppi.obj
      }
      
    } ## End of broad ppi umbrella
    
    #################################
    ## Promoter Objective Function ##
    if("promoter" %in% obj.fun) {
      promoter.obj <- df %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
        mutate(promoter = replace_na(promoter, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(promoter = sum(promoter) / as.numeric(max.hit[which(max.hit$of == "promoter"), "max_loci"]))
      
      obj.fun.res[["promoter"]] <- promoter.obj
    }
    
    #############################
    ## ATAC Objective Function ##
    ## Marker atac peak 
    if("marker.atac" %in% obj.fun) {
      m.atac.obj <- df %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
        mutate(marker.atac = replace_na(marker.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(marker.atac = sum(marker.atac) / as.numeric(max.hit[which(max.hit$of == "marker.atac"), "max_loci"]))
      obj.fun.res[["marker.atac"]] <- m.atac.obj
    }
    ## Common atac peak
    if("common.atac" %in% obj.fun) {
      c.atac.obj <- df %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
        mutate(common.atac = replace_na(common.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(common.atac = sum(common.atac) / as.numeric(max.hit[which(max.hit$of == "common.atac"), "max_loci"]))
      obj.fun.res[["common.atac"]] <- c.atac.obj
    }
    ########################
    ## lncrna interaction ##
    if("lncrna" %in% obj.fun) {
      lncrna.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- rna.protein.map %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                                   relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x)
      }, mc.cores = n.cores)
      lncrna.obj <- as.data.frame(do.call(rbind, lncrna.obj))
      
      all.lncrna.obj <- mclapply(X = as.list(unique(lncrna.obj$n)), FUN = function(x) {
        lncrna.df <- lncrna.obj %>% dplyr::filter(n == x)
        lncrna.genes <- data.frame(gene = unique(c(lncrna.df$gene1, lncrna.df$gene2)),
                                   lncrna = 1)
        lncrna.res <- df %>% dplyr::filter(n == x) %>%
          left_join(x = ., y = lncrna.genes, by = "gene") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0)) %>%
          dplyr::summarise(lncrna = sum(lncrna) / as.numeric(max.hit[which(max.hit$of == "lncrna"), "max_loci"])) %>%
          mutate(n = x) %>% .[, c("n", "lncrna")]
      }, mc.cores = n.cores)
      all.lncrna.obj <- as.data.frame(do.call(rbind, all.lncrna.obj))
      
      if(nrow(all.lncrna.obj) == 0) { ## If no lncrna interaction then do this
        lncrna.res <- data.frame(n = unique(df$n),
                                 lncrna = 0)
      } else { ## If lncrna interaction present
        lncrna.res <- data.frame(n = unique(df$n)) %>%
          left_join(x = ., y = all.lncrna.obj, by = "n") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0))
      }
      obj.fun.res[["lncrna"]] <- lncrna.res
    } ## end of lncrna
    
    ######################
    ## tfbs interaction ##
    if("tfbs" %in% obj.fun) {
      
      tfbs.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n) %>%
          left_join(x = ., y = ct_exp, by = c("celltype"="type"), relationship = "many-to-many") %>% ## Add celltype and tf gene info
          left_join(x = ., y = mb.tfbs.disrupt, by = c("gene_name"="geneSymbol", "locus"="LEAD_SNP")) %>%
          dplyr::select(-gene_name) %>%
          dplyr::group_by(locus, celltype, gene) %>%
          mutate(tfbs = replace_na(tfbs, replace = 0)) %>% ## Add 0 where NA is
          dplyr::summarise(tfbs = sum(tfbs)) %>% ## sum all tfbs disrupt hits
          mutate(tfbs = case_when(tfbs > 0 ~ 1, ## boolean value conversion
                                  tfbs == 0 ~ 0),
                 n = x) %>%
          dplyr::select(celltype, locus, gene, n, tfbs) %>%
          dplyr::ungroup() 
        
      }, mc.cores = n.cores)
      tfbs.obj <- as.data.frame(do.call(rbind, tfbs.obj))
      tfbs.obj <- tfbs.obj %>% dplyr::group_by(n) %>% dplyr::summarise(tfbs = sum(tfbs) / as.numeric(max.hit[which(max.hit$of == "tfbs"), "max_loci"]))
      obj.fun.res[["tfbs"]] <- tfbs.obj
    } ## end of tfbs
    
    if("tfbs.marker" %in% obj.fun) {
      tfbs.marker.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with TF expressed data
        proposal.set <- ct_exp %>% mutate(tf.hit = 1) %>% 
          left_join(x = proposal.set, y = ., by = c("gene"="gene_name", "celltype"="type")) %>%
          mutate(tf.hit = replace_na(tf.hit, replace = 0))
        
        ## If there are any TF expressed in correct cell type then add in motifbreakR info
        if(any(proposal.set$tf.hit == 1)) {
          ## What tf.gene are a hit
          tf.genes <- proposal.set %>% dplyr::filter(tf.hit == 1) %>% pull(gene)
          ## What SNPs have motif broken?
          for(tf.gene in tf.genes) {
            tfbs.disrupt <- mb.tfbs.disrupt %>% dplyr::filter(geneSymbol %in% tf.gene)
            ## Now add a point for all lead SNPs with a MB disruption for that tf.gene
            for(tfbs in unique(tfbs.disrupt$LEAD_SNP)) {
              proposal.set[which(proposal.set$locus == tfbs), "tf.hit"] <- 1
            }
          }
        } ## end of IF condition 
        proposal.set <- proposal.set %>% mutate(n = x) %>% ## add back generation info
          dplyr::select(celltype, locus, gene, n, tf.hit) %>%
          dplyr::group_by(n) %>%
          dplyr::summarise(tfbs = sum(tf.hit) / as.numeric(max.hit[which(max.hit$of == "tf.hit"), "max_loci"]))
        
      }, mc.cores = n.cores)
      tfbs.marker.obj <- as.data.frame(do.call(rbind, tfbs.marker.obj))
    }
    
    #####################################
    ## Fitness mean objective function ##
    
    for(n in 1:length(obj.fun)) {
      if(n == 1) {
        fitness.obj <- obj.fun.res[[obj.fun[n]]]
      } else {
        add_test <- obj.fun.res[[obj.fun[n]]]
        fitness.obj <- left_join(x = fitness.obj, y = add_test, by = "n")
      }
    }
    fitness.obj[is.na(fitness.obj)] <- 0
    fitness.obj$fitness <- apply(X = fitness.obj %>% dplyr::select(-n), MARGIN = 1, FUN = function(x) {
      mean(x)
    })
    fitness.obj$gen <- generation - 1 
    
    ######################################
    ## Rank the proposals into 5 groups ##
    
    ranked.proposals <- fitness.obj %>% dplyr::select(n, fitness) %>% arrange(., desc(fitness)) %>%
      mutate(quantile = rep(x = 5:1, each = nrow(.)/5))
    ## Prob of quantiles
    group.prob <- c(0.02,0.03,0.05,0.1,0.8)
    names(group.prob) <- 1:max(ranked.proposals$quantile)
    
    #####################
    ## Updating top.df ##
    
    if(generation-1 == 0) { 
      top.df$gen <- generation - 1
      fitness.df <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
    } else {
      new.fitness <- ranked.proposals %>% dplyr::select(-quantile) %>% mutate(gen = generation - 1)
      ## Update top rank proposal df with new generation data. Basically, replace proposals that are weaker than the new proposals
      fitness.df <- rbind(fitness.df, new.fitness) %>% arrange(desc(fitness)) %>% .[1:1000,]
      ## Update top df proposal elements based on fitness.df
      keep.list <- fitness.df %>% dplyr::select(n, gen) %>% mutate(keep = 1)
      current.proposal <- df %>% mutate(gen = generation - 1)
      top.df <- rbind(top.df, current.proposal) %>% left_join(x = ., y = keep.list, by = c("n", "gen")) %>%
        dplyr::filter(keep == 1) %>%
        dplyr::select(-keep)
    }
    
    ###################################
    ## Randomly Select 100 Proposals ##
    
    if("progressive" %in% selection.method) {
      
      ## Select 100 parents based on progressive scale
      sel.proposal <- progressive_quantile(group.prob = group.prob, n_proposal = parent.size, proposal = ranked.proposals)
    }
    
    if("lexicase" %in% selection.method) {
      sel.proposal <- list()
      
      for(lexicase.run in seq(parent.size)) { ## Find parents based on parent.size
        shuff.of <- sample(x = obj.fun, size = length(obj.fun), replace = FALSE) ## Shuffle obj func
        pool <- fitness.obj ## Store fitness.obj as pool 
        shuff.counter <- 1 ## Counter used for moving through the shuff.of
        while(nrow(pool) > 0 & shuff.counter <= length(shuff.of)) { ## While pool is greater than 0 and shuff.counter less than shuff.of
          message(paste0("Working on parent ", lexicase.run, ": Shuffle OF ", shuff.counter, " is ", shuff.of[[shuff.counter]]), ". Pool is ", nrow(pool))
          col.index <- which(colnames(pool) == shuff.of[[shuff.counter]])
          max.score <- pool %>% pull(shuff.of[[shuff.counter]]) %>% max() 
          pool <- pool %>% filter_at(col.index, all_vars(. == max.score)) 
          
          ## If only one parent then stop selection
          if(nrow(pool) == 1) { 
            sel.proposal[[lexicase.run]] <- pool %>% dplyr::select(n, gen)
            message(paste0("Parent ", lexicase.run, " found. Only one left. Stopped at shuff.counter ", shuff.counter))
            break
          } else if(nrow(pool) > 1) { ## If pool is still greater than 1
            if(shuff.counter == length(shuff.of)) {
              sel.proposal[[lexicase.run]] <- pool %>% .[sample(x = seq(nrow(pool)), size = 1),] %>% dplyr::select(n, gen)
              message(paste0("Parent ", lexicase.run, " found. More than one. Random select Stopped at shuff.counter ", shuff.counter))
              break
            } else {
              shuff.counter <- shuff.counter + 1
            }
          } 
        } ## end of while loop
      } ## end of lexicase selection
      sel.proposal <- as.data.frame(do.call(rbind, sel.proposal)) ## Process parent selection result
      sel.proposal <- sel.proposal %>% left_join(x = ., y = ranked.proposals, by = "n") %>% dplyr::select(-gen)
    } ## end of lexicase function
    
    
    
    ############
    ## Mating ##
    new.proposal <- mating.pair(sel.proposal = sel.proposal, n_proposal = proposal.size, n.pairs = parent.size, df = df, khan.method = khan.method, fitness.df = fitness.df, top.df = top.df,
                                generation = generation, n.loci = n.loci, n.cores = n.cores)
    new.proposal <- as.data.frame(do.call(rbind,new.proposal)) %>% remove_rownames() %>%
      mutate(n = rep(1:proposal.size, each = n.loci))
    ##############
    ## Mutation ##
    new.proposal <- add.mutations(new.proposal = new.proposal, n.cores = n.cores)
    
    ## Insert khan to the proposal
    rand.kick <- new.proposal %>% rownames_to_column(var = "row") %>% slice_sample(n = 1) %>% .$n
    new.proposal <- new.proposal %>% dplyr::filter(n != rand.kick)
    
    khan <- fitness.df[1,] ## Which gen and proposal is khan?
    khan.df <- top.df %>% dplyr::filter(n == khan$n, ## Get khan proposals
                                        gen == khan$gen) %>%
      dplyr::select(-gen) %>%
      mutate(n = proposal.size + generation)
    
    new.proposal <- rbind(khan.df, new.proposal)
    
    ###########################################
    ## Remove one random and add khan inside ##
    
    df <- new.proposal ## Overwrite multi_proposal with new.proposal
    
    ###############
    ## Store res ##
    fitness.res.list[[generation]] <- fitness.obj
    
    
  } ## End loop for generations
  fitness.res.list[["proposal"]] <- last.gen
  fitness.res.list[["top"]] <- top.df
  fitness.res.list[["fitness"]] <-fitness.df
  return(fitness.res.list)
} ## End of wrapper

#######################
## Calculate fitness ##

calc_fitness <- function(df, obj.fun, n.cores, method) {
  ## Storing results from obj fun for fitness function
  obj.fun.res <- list()
  
  if(method == "scaled") {
    #############################
    ## CEVM Objective Function ##
    if("cevm" %in% obj.fun) { ## If cevm then do the following
      cevm <- cevm.genes %>% mutate(cevm = 1)
      
      cevm.obj <- df %>% left_join(x = ., y = cevm, by = c("gene", "celltype")) %>%
        mutate(cevm = replace_na(cevm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cevm = sum(cevm) / as.numeric(max.hit[which(max.hit$of == "cevm"), "max_loci"]))
      obj.fun.res[["cevm"]] <- as.data.frame(cevm.obj)
    }
    
    ####################################
    ## FindMarkers Objective Function ##
    if("findmarkers" %in% obj.fun) {
      fm <- marker.hits %>% mutate(fm = 1)
      
      fm.obj <- df %>% left_join(x = ., y = fm, by = c("gene", "celltype")) %>%
        mutate(fm = replace_na(fm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(fm = sum(fm) / as.numeric(max.hit[which(max.hit$of == "fm"), "max_loci"]))
      obj.fun.res[["findmarkers"]] <- as.data.frame(fm.obj)
    }
    
    ########################################
    ## Combined gene marker - CEVM and FM ##
    if("marker.gene" %in% obj.fun) {
      ## marker gene data
      marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
      
      ## Data process
      marker.gene.obj <- df %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>% ## join marker.gene to proposal
        mutate(marker.gene = replace_na(marker.gene, replace = 0)) %>%  ## replace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(marker.gene = sum(marker.gene) / as.numeric(max.hit[which(max.hit$of == "marker.gene"), "max_loci"])) ## get % pos hit
      obj.fun.res[["marker.gene"]] <- as.data.frame(marker.gene.obj) ## store in list
    }
    ##############################
    ## Magma Objective Function ##
    if("magma" %in% obj.fun) {
      magma.obj <- df %>% left_join(x = ., y = magma, by = "gene") %>% ## join magma to proposal
        mutate(magma = replace_na(magma, replace = 0)) %>%  ## replace NA with 0 
        dplyr::group_by(n) %>% dplyr::summarise(magma = sum(magma) / as.numeric(max.hit[which(max.hit$of == "magma"), "max_loci"])) ## get % pos hit
      obj.fun.res[["magma"]] <- as.data.frame(magma.obj) ## store in list
    }
    
    ####################################
    ## Cancer Gene Objective Function ## 
    if("cancer.gene" %in% obj.fun) {
      hi.mut.obj <- df %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>% ## join cancer.gene to proposal
        mutate(cancer.gene = replace_na(cancer.gene, replace = 0)) %>%  ## rpelace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(cancer.gene = sum(cancer.gene) / as.numeric(max.hit[which(max.hit$of == "cancer.gene"), "max_loci"])) ## get % pos hit
      obj.fun.res[["cancer.gene"]] <- as.data.frame(hi.mut.obj) ## store in list
    }
    
    ############################
    ## PPI Objective Function ##
    if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
      ppi.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- ppi %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                       relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x) 
        
      }, mc.cores = n.cores)
      ppi.obj <- as.data.frame(do.call(rbind, ppi.obj))
      
      #######################
      ## Intracellular PPI ##
      if("intra.ppi" %in% obj.fun) {
        intra.df <-  ppi.obj %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(intra.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        intra.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
          intra.ppi.genes <- intra.df %>% dplyr::filter(n == x)
          if(nrow(intra.ppi.genes) > 0) {
            intra.ppi.genes <- data.frame(gene = unique(c(intra.ppi.genes$gene1, intra.ppi.genes$gene2)),
                                          intra.ppi = 1)
            ## Get specific proposal set 
            intra.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = intra.ppi.genes, by = "gene") %>%
              mutate(intra.ppi = replace_na(intra.ppi, replace = 0)) %>%
              summarise(intra.ppi = sum(intra.ppi) /as.numeric(max.hit[which(max.hit$of == "intra.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "intra.ppi")]
          } else {
            intra.res <- data.frame(n = x,
                                    intra.ppi = 0)
          }
        }, mc.cores = n.cores)
        intra.obj <- as.data.frame(do.call(rbind, intra.obj))
        obj.fun.res[["intra.ppi"]] <- intra.obj
      } ## End of intra ppi
      
      #############
      ## All PPI ##
      if("ppi" %in% obj.fun) {
        all.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          ppi.df <- ppi.obj %>% dplyr::filter(n == x)
          ppi.genes <- data.frame(gene = unique(c(ppi.df$gene1, ppi.df$gene2)),
                                  ppi = 1)
          ppi.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = ppi.genes, by = "gene") %>%
            mutate(ppi = replace_na(ppi, replace = 0)) %>%
            dplyr::summarise(ppi = sum(ppi) / as.numeric(max.hit[which(max.hit$of == "ppi"), "max_loci"])) %>%
            mutate(n = x) %>% .[, c("n", "ppi")]
        }, mc.cores = n.cores)
        all.ppi.obj <- as.data.frame(do.call(rbind, all.ppi.obj))
        ppi.res <- data.frame(n = unique(df$n), 
                              ppi = 0)
        for(r in unique(all.ppi.obj$n)) {
          ppi.res[which(ppi.res$n == r), "ppi"] <- all.ppi.obj[which(all.ppi.obj$n == r), "ppi"]
        }
        obj.fun.res[["ppi"]] <- ppi.res
      } ## End of ppi
      
      ################
      ## Marker PPI ##
      if("marker.ppi" %in% obj.fun) {
        all.markers <- combined.gene.marker %>% mutate(marker.gene = 1)
        marker.obj <- ppi.obj %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() #%>% ## Keep positive hits.
        
        if(nrow(marker.obj) == 0) {
          marker.ppi.obj <- data.frame(n = unique(df$n),
                                       marker.ppi = 0)
        } else {
          marker.ppi.obj <- mclapply(X = as.list(unique(marker.obj$n)), FUN = function(x) {
            marker.df <- marker.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
            
            marker.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
              mutate(marker.ppi = replace_na(marker.ppi, replace = 0)) %>%
              dplyr::summarise(marker.ppi = sum(marker.ppi) / as.numeric(max.hit[which(max.hit$of == "marker.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "marker.ppi")]
          }, mc.cores = n.cores)
          marker.ppi.obj <- as.data.frame(do.call(rbind, marker.ppi.obj))
          
          marker.ppi.obj <- data.frame(n = unique(df$n)) %>%
            left_join(x = ., y = marker.ppi.obj, by = "n") %>%
            mutate(marker.ppi = replace_na(marker.ppi, replace = 0))
        }
        
        obj.fun.res[["marker.ppi"]] <- marker.ppi.obj
      } ## End of marker.ppi
      
      #######################
      ## Intercellular PPI ##
      if("inter.ppi" %in% obj.fun) {
        marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        
        marker.obj <- ppi.obj %>% left_join(x = ., y = marker.gene, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = marker.gene, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() %>% ## Keep positive hits.
          dplyr::select(-marker.gene.x, -marker.gene.y) ## Get rid of columns no longer needed
        
        if(nrow(marker.obj) == 0) {
          inter.ppi.obj <- data.frame(n = unique(df$n),
                                      inter.ppi = 0)
        } else {
          ## First identity all source genes in celltalkdb. We will check if the source gene is found in 
          ## gene1 or gene2. If no source gene in either then remove. If source gene is found in both 
          ## Then remove. We will be filtering for rowsums == 1 on source_genesymbol1 and 2.
          inter.obj <- marker.obj %>% left_join(x = ., 
                                                y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                                by = c("gene1"="source_genesymbol"),
                                                relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("source_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                      by = c("gene2"="source_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>% 
            dplyr::rename("source_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
            
            ## Now trying to find targets the same way we did for source genes
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene1"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene2"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
            ## Removing target and source identity column
            dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
          
          inter.ppi.obj <- mclapply(X = as.list(unique(inter.obj$n)), FUN = function(x) {
            marker.df <- inter.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
            
            inter.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0)) %>%
              dplyr::summarise(inter.ppi = sum(inter.ppi) / as.numeric(max.hit[which(max.hit$of == "inter.ppi"), "max_loci"])) %>%
              mutate(n = x) %>% .[, c("n", "inter.ppi")]
          }, mc.cores = n.cores)
          inter.ppi.obj <- as.data.frame(do.call(rbind, inter.ppi.obj)) 
          
          if(nrow(inter.ppi.obj) == 0) {
            inter.ppi.obj <- data.frame(n = unique(df$n),
                                        inter.ppi = 0)
          } else {
            inter.ppi.obj <- data.frame(n = unique(df$n)) %>%
              left_join(x = ., y = inter.ppi.obj, by = "n") %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0))
          }
        }
        
        obj.fun.res[["inter.ppi"]] <- inter.ppi.obj
      } ## End of inter.ppi
      
      ###############
      ## Multi PPI ##
      if("multi.ppi" %in% obj.fun) {
        multi.obj <- ppi.obj %>% dplyr::select(-celltype1,-celltype2) %>% distinct()
        multi.res <- mclapply(X = as.list(unique(multi.obj$n)), FUN = function(x) {
          multi.df <- multi.obj %>% dplyr::filter(n %in% x)
          df1 <- multi.df %>% dplyr::select(gene1, locus1) %>% dplyr::rename("gene"="gene1", "locus"="locus1")
          df2 <- multi.df %>% dplyr::select(gene2, locus2) %>% dplyr::rename("gene"="gene2", "locus"="locus2")
          gene.locus.set <- rbind(df1, df2) %>% mutate(set = paste(gene, locus, sep = "_")) %>%  ## Combine df1 and df2 and create set name
            group_by(set) %>% ## Group by set
            tally() %>% ## ## Count up how many times the set appears
            dplyr::filter(n > 1) %>% nrow() ## Select gene locus combination with more than one interaction
          res <- data.frame(n = x, multi.ppi = gene.locus.set / n.loci)
        }, mc.cores = n.cores)
        multi.res <- as.data.frame(do.call(rbind, multi.res))
        obj.fun.res[["multi.ppi"]] <- multi.res
      } ## End of multi.ppi
      
      if("go.cc.ppi" %in% obj.fun) {
        go.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          go.ppi.df <- ppi.obj %>% dplyr::filter(n == x) %>% dplyr::select(-celltype1, -celltype2) %>%
            left_join(x = ., ## Add cellular component info for all genes in gene1 (NA means no info on gene)
                      y = map.go.gene[,c("external_gene_name", "name_1006")], 
                      by = c("gene1"="external_gene_name"), 
                      relationship = "many-to-many") %>%
            dplyr::rename("cc1"="name_1006") %>% ## Rename cc column as cc1 to link to gene1
            left_join(x = ., ## Add cc info for gene2
                      y = map.go.gene[,c("external_gene_name", "name_1006")],
                      by = c("gene2"="external_gene_name"),
                      relationship = "many-to-many") %>%
            dplyr::rename("cc2"="name_1006") %>% ## Rename cc column as cc2 to link to gene2
            .[which(.$cc1 == .$cc2),] ## Filter for cases where gene1 and gene2 are in same cellular compartment
          ## Get all genes from go.ppi.df
          go.genes <- data.frame(gene = unique(c(go.ppi.df$gene1, go.ppi.df$gene2)),
                                 go.cc = 1)
          go.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = go.genes, by = "gene") %>%
            mutate(go.cc = replace_na(go.cc, replace = 0)) %>%
            dplyr::summarise(go.cc = sum(go.cc) / as.numeric(max.hit[which(max.hit$of == "go.cc"), "max_loci"])) %>%
            mutate(n = x) %>% .[, c("n", "go.cc")]
        }, mc.cores = n.cores)
        go.ppi.obj <- as.data.frame(do.call(rbind, go.ppi.obj))
        
        obj.fun.res[["go.cc.ppi"]] <- go.ppi.obj
      }
      
    } ## End of broad ppi umbrella
    
    #################################
    ## Promoter Objective Function ##
    if("promoter" %in% obj.fun) {
      promoter.obj <- df %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
        mutate(promoter = replace_na(promoter, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(promoter = sum(promoter)/as.numeric(max.hit[which(max.hit$of == "promoter"), "max_loci"]))
      
      obj.fun.res[["promoter"]] <- promoter.obj
    }
    
    #############################
    ## ATAC Objective Function ##
    ## Marker atac peak 
    if("marker.atac" %in% obj.fun) {
      m.atac.obj <- df %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
        mutate(marker.atac = replace_na(marker.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(marker.atac = sum(marker.atac)/as.numeric(max.hit[which(max.hit$of == "marker.atac"), "max_loci"]))
      obj.fun.res[["marker.atac"]] <- m.atac.obj
    }
    ## Common atac peak
    if("common.atac" %in% obj.fun) {
      c.atac.obj <- df %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
        mutate(common.atac = replace_na(common.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(common.atac = sum(common.atac)/as.numeric(max.hit[which(max.hit$of == "common.atac"), "max_loci"]))
      obj.fun.res[["common.atac"]] <- c.atac.obj
    }
    ########################
    ## lncrna interaction ##
    if("lncrna" %in% obj.fun) {
      lncrna.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- rna.protein.map %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                                   relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x)
      }, mc.cores = n.cores)
      lncrna.obj <- as.data.frame(do.call(rbind, lncrna.obj))
      
      all.lncrna.obj <- mclapply(X = as.list(unique(lncrna.obj$n)), FUN = function(x) {
        lncrna.df <- lncrna.obj %>% dplyr::filter(n == x)
        lncrna.genes <- data.frame(gene = unique(c(lncrna.df$gene1, lncrna.df$gene2)),
                                   lncrna = 1)
        lncrna.res <- df %>% dplyr::filter(n == x) %>%
          left_join(x = ., y = lncrna.genes, by = "gene") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0)) %>%
          dplyr::summarise(lncrna = sum(lncrna) / as.numeric(max.hit[which(max.hit$of == "lncrna"), "max_loci"])) %>%
          mutate(n = x) %>% .[, c("n", "lncrna")]
      }, mc.cores = n.cores)
      all.lncrna.obj <- as.data.frame(do.call(rbind, all.lncrna.obj))
      
      if(nrow(all.lncrna.obj) == 0) { ## If no lncrna interaction then do this
        lncrna.res <- data.frame(n = unique(df$n),
                                 lncrna = 0)
      } else { ## If lncrna interaction present
        lncrna.res <- data.frame(n = unique(df$n)) %>%
          left_join(x = ., y = all.lncrna.obj, by = "n") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0))
      }
      
      obj.fun.res[["lncrna"]] <- lncrna.res
    } ## end of lncrna
    
    ######################
    ## tfbs interaction ##
    if("tfbs" %in% obj.fun) {
      
      tfbs.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n) %>%
          left_join(x = ., y = ct_exp, by = c("celltype"="type"), relationship = "many-to-many") %>% ## Add celltype and tf gene info
          left_join(x = ., y = mb.tfbs.disrupt, by = c("gene_name"="geneSymbol", "locus"="LEAD_SNP")) %>%
          dplyr::select(-gene_name) %>%
          dplyr::group_by(locus, celltype, gene) %>%
          mutate(tfbs = replace_na(tfbs, replace = 0)) %>% ## Add 0 where NA is
          dplyr::summarise(tfbs = sum(tfbs)) %>% ## sum all tfbs disrupt hits
          mutate(tfbs = case_when(tfbs > 0 ~ 1, ## boolean value conversion
                                  tfbs == 0 ~ 0),
                 n = x) %>%
          dplyr::select(celltype, locus, gene, n, tfbs) %>%
          dplyr::ungroup() 
        
      }, mc.cores = n.cores)
      tfbs.obj <- as.data.frame(do.call(rbind, tfbs.obj))
      tfbs.obj <- tfbs.obj %>% dplyr::group_by(n) %>% dplyr::summarise(tfbs = sum(tfbs)/as.numeric(max.hit[which(max.hit$of == "tfbs"), "max_loci"]))
      obj.fun.res[["tfbs"]] <- tfbs.obj
    } ## end of tfbs
    
    if("tfbs.marker" %in% obj.fun) {
      tfbs.marker.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with TF expressed data
        proposal.set <- ct_exp %>% mutate(tf.hit = 1) %>% 
          left_join(x = proposal.set, y = ., by = c("gene"="gene_name", "celltype"="type")) %>%
          mutate(tf.hit = replace_na(tf.hit, replace = 0))
        
        ## If there are any TF expressed in correct cell type then add in motifbreakR info
        if(any(proposal.set$tf.hit == 1)) {
          ## What tf.gene are a hit
          tf.genes <- proposal.set %>% dplyr::filter(tf.hit == 1) %>% pull(gene)
          ## What SNPs have motif broken?
          for(tf.gene in tf.genes) {
            tfbs.disrupt <- mb.tfbs.disrupt %>% dplyr::filter(geneSymbol %in% tf.gene)
            ## Now add a point for all lead SNPs with a MB disruption for that tf.gene
            for(tfbs in unique(tfbs.disrupt$LEAD_SNP)) {
              proposal.set[which(proposal.set$locus == tfbs), "tf.hit"] <- 1
            }
          }
        } ## end of IF condition 
        proposal.set <- proposal.set %>% mutate(n = x) %>% ## add back generation info
          dplyr::select(celltype, locus, gene, n, tf.hit) %>%
          dplyr::group_by(n) %>%
          dplyr::summarise(tfbs = sum(tf.hit)/as.numeric(max.hit[which(max.hit$of == "tf.hit"), "max_loci"]))
        
      }, mc.cores = n.cores)
      tfbs.marker.obj <- as.data.frame(do.call(rbind, tfbs.marker.obj))
    }
  }
  
  if(method == "unscaled") {
    #############################
    ## CEVM Objective Function ##
    if("cevm" %in% obj.fun) { ## If cevm then do the following
      cevm <- cevm.genes %>% mutate(cevm = 1)
      
      cevm.obj <- df %>% left_join(x = ., y = cevm, by = c("gene", "celltype")) %>%
        mutate(cevm = replace_na(cevm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(cevm = sum(cevm) / length(cevm))
      obj.fun.res[["cevm"]] <- as.data.frame(cevm.obj)
    }
    
    ####################################
    ## FindMarkers Objective Function ##
    if("findmarkers" %in% obj.fun) {
      fm <- marker.hits %>% mutate(fm = 1)
      
      fm.obj <- df %>% left_join(x = ., y = fm, by = c("gene", "celltype")) %>%
        mutate(fm = replace_na(fm, replace = 0)) %>% 
        dplyr::group_by(n) %>% dplyr::summarise(fm = sum(fm) / length(fm))
      obj.fun.res[["findmarkers"]] <- as.data.frame(fm.obj)
    }
    
    ########################################
    ## Combined gene marker - CEVM and FM ##
    if("marker.gene" %in% obj.fun) {
      ## marker gene data
      marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
      
      ## Data process
      marker.gene.obj <- df %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>% ## join marker.gene to proposal
        mutate(marker.gene = replace_na(marker.gene, replace = 0)) %>%  ## replace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(marker.gene = sum(marker.gene) / length(marker.gene)) ## get % pos hit
      obj.fun.res[["marker.gene"]] <- as.data.frame(marker.gene.obj) ## store in list
    }
    ##############################
    ## Magma Objective Function ##
    if("magma" %in% obj.fun) {
      magma.obj <- df %>% left_join(x = ., y = magma, by = "gene") %>% ## join magma to proposal
        mutate(magma = replace_na(magma, replace = 0)) %>%  ## replace NA with 0 
        dplyr::group_by(n) %>% dplyr::summarise(magma = sum(magma) / length(magma)) ## get % pos hit
      obj.fun.res[["magma"]] <- as.data.frame(magma.obj) ## store in list
    }
    
    ####################################
    ## Cancer Gene Objective Function ## 
    if("cancer.gene" %in% obj.fun) {
      hi.mut.obj <- df %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>% ## join cancer.gene to proposal
        mutate(cancer.gene = replace_na(cancer.gene, replace = 0)) %>%  ## rpelace NA with 0
        dplyr::group_by(n) %>% dplyr::summarise(cancer.gene = sum(cancer.gene) / length(cancer.gene)) ## get % pos hit
      obj.fun.res[["cancer.gene"]] <- as.data.frame(hi.mut.obj) ## store in list
    }
    
    ############################
    ## PPI Objective Function ##
    if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
      ppi.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- ppi %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                       relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x) 
        
      }, mc.cores = n.cores)
      ppi.obj <- as.data.frame(do.call(rbind, ppi.obj))
      
      #######################
      ## Intracellular PPI ##
      if("intra.ppi" %in% obj.fun) {
        intra.df <-  ppi.obj %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
          ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
          dplyr:::filter(intra.commun == 1) %>%
          group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                               pmin(gene1, gene2), sep = "_")) %>%
          dplyr::filter(!grp %in% celltalk_gene$grp) %>%
          ungroup() %>%
          dplyr::select(-grp) %>%
          distinct()
        intra.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
          intra.ppi.genes <- intra.df %>% dplyr::filter(n == x)
          if(nrow(intra.ppi.genes) > 0) {
            intra.ppi.genes <- data.frame(gene = unique(c(intra.ppi.genes$gene1, intra.ppi.genes$gene2)),
                                          intra.ppi = 1)
            ## Get specific proposal set 
            intra.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = intra.ppi.genes, by = "gene") %>%
              mutate(intra.ppi = replace_na(intra.ppi, replace = 0)) %>%
              summarise(intra.ppi = sum(intra.ppi) /length(intra.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "intra.ppi")]
          } else {
            intra.res <- data.frame(n = x,
                                    intra.ppi = 0)
          }
        }, mc.cores = n.cores)
        intra.obj <- as.data.frame(do.call(rbind, intra.obj))
        obj.fun.res[["intra.ppi"]] <- intra.obj
      } ## End of intra ppi
      
      #############
      ## All PPI ##
      if("ppi" %in% obj.fun) {
        all.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          ppi.df <- ppi.obj %>% dplyr::filter(n == x)
          ppi.genes <- data.frame(gene = unique(c(ppi.df$gene1, ppi.df$gene2)),
                                  ppi = 1)
          ppi.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = ppi.genes, by = "gene") %>%
            mutate(ppi = replace_na(ppi, replace = 0)) %>%
            dplyr::summarise(ppi = sum(ppi) / length(ppi)) %>%
            mutate(n = x) %>% .[, c("n", "ppi")]
        }, mc.cores = n.cores)
        all.ppi.obj <- as.data.frame(do.call(rbind, all.ppi.obj))
        ppi.res <- data.frame(n = unique(df$n), 
                              ppi = 0)
        for(r in unique(all.ppi.obj$n)) {
          ppi.res[which(ppi.res$n == r), "ppi"] <- all.ppi.obj[which(all.ppi.obj$n == r), "ppi"]
        }
        obj.fun.res[["ppi"]] <- ppi.res
      } ## End of ppi
      
      ################
      ## Marker PPI ##
      if("marker.ppi" %in% obj.fun) {
        all.markers <- combined.gene.marker %>% mutate(marker.gene = 1)
        marker.obj <- ppi.obj %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() #%>% ## Keep positive hits.
        
        if(nrow(marker.obj) == 0) {
          marker.ppi.obj <- data.frame(n = unique(df$n),
                                       marker.ppi = 0)
        } else {
          marker.ppi.obj <- mclapply(X = as.list(unique(marker.obj$n)), FUN = function(x) {
            marker.df <- marker.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
            
            marker.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
              mutate(marker.ppi = replace_na(marker.ppi, replace = 0)) %>%
              dplyr::summarise(marker.ppi = sum(marker.ppi) / length(marker.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "marker.ppi")]
          }, mc.cores = n.cores)
          marker.ppi.obj <- as.data.frame(do.call(rbind, marker.ppi.obj))
          
          marker.ppi.obj <- data.frame(n = unique(df$n)) %>%
            left_join(x = ., y = marker.ppi.obj, by = "n") %>%
            mutate(marker.ppi = replace_na(marker.ppi, replace = 0))
        }
        
        obj.fun.res[["marker.ppi"]] <- marker.ppi.obj
      } ## End of marker.ppi
      
      #######################
      ## Intercellular PPI ##
      if("inter.ppi" %in% obj.fun) {
        marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
        
        marker.obj <- ppi.obj %>% left_join(x = ., y = marker.gene, by = c("gene1"="gene", "celltype1"="celltype")) %>%
          left_join(x = ., y = marker.gene, by = c("gene2"="gene", "celltype2"="celltype")) %>%
          na.omit() %>% ## Keep positive hits.
          dplyr::select(-marker.gene.x, -marker.gene.y) ## Get rid of columns no longer needed
        
        if(nrow(marker.obj) == 0) {
          inter.ppi.obj <- data.frame(n = unique(df$n),
                                      inter.ppi = 0)
        } else {
          ## First identity all source genes in celltalkdb. We will check if the source gene is found in 
          ## gene1 or gene2. If no source gene in either then remove. If source gene is found in both 
          ## Then remove. We will be filtering for rowsums == 1 on source_genesymbol1 and 2.
          inter.obj <- marker.obj %>% left_join(x = ., 
                                                y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                                by = c("gene1"="source_genesymbol"),
                                                relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("source_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                      by = c("gene2"="source_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>% 
            dplyr::rename("source_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
            
            ## Now trying to find targets the same way we did for source genes
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene1"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene1" = "inter.ppi") %>%
            left_join(x = ., 
                      y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                      by = c("gene2"="target_genesymbol"),
                      relationship = "many-to-many") %>%
            distinct() %>%
            dplyr::rename("target_gene2" = "inter.ppi") %>%
            .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
            ## Removing target and source identity column
            dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
          
          inter.ppi.obj <- mclapply(X = as.list(unique(inter.obj$n)), FUN = function(x) {
            marker.df <- inter.obj %>% dplyr::filter(n == x)
            df1 <- marker.df %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            df2 <- marker.df %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
            inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
            
            inter.res <- df %>% dplyr::filter(n == x) %>%
              left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0)) %>%
              dplyr::summarise(inter.ppi = sum(inter.ppi) / length(inter.ppi)) %>%
              mutate(n = x) %>% .[, c("n", "inter.ppi")]
          }, mc.cores = n.cores)
          inter.ppi.obj <- as.data.frame(do.call(rbind, inter.ppi.obj)) 
          
          if(nrow(inter.ppi.obj) == 0) {
            inter.ppi.obj <- data.frame(n = unique(df$n),
                                        inter.ppi = 0)
          } else {
            inter.ppi.obj <- data.frame(n = unique(df$n)) %>%
              left_join(x = ., y = inter.ppi.obj, by = "n") %>%
              mutate(inter.ppi = replace_na(inter.ppi, replace = 0))
          }
        }
        
        obj.fun.res[["inter.ppi"]] <- inter.ppi.obj
      } ## End of inter.ppi
      
      ###############
      ## Multi PPI ##
      if("multi.ppi" %in% obj.fun) {
        multi.obj <- ppi.obj %>% dplyr::select(-celltype1,-celltype2) %>% distinct()
        multi.res <- mclapply(X = as.list(unique(multi.obj$n)), FUN = function(x) {
          multi.df <- multi.obj %>% dplyr::filter(n %in% x)
          df1 <- multi.df %>% dplyr::select(gene1, locus1) %>% dplyr::rename("gene"="gene1", "locus"="locus1")
          df2 <- multi.df %>% dplyr::select(gene2, locus2) %>% dplyr::rename("gene"="gene2", "locus"="locus2")
          gene.locus.set <- rbind(df1, df2) %>% mutate(set = paste(gene, locus, sep = "_")) %>%  ## Combine df1 and df2 and create set name
            group_by(set) %>% ## Group by set
            tally() %>% ## ## Count up how many times the set appears
            dplyr::filter(n > 1) %>% nrow() ## Select gene locus combination with more than one interaction
          res <- data.frame(n = x, multi.ppi = gene.locus.set / n.loci)
        }, mc.cores = n.cores)
        multi.res <- as.data.frame(do.call(rbind, multi.res))
        obj.fun.res[["multi.ppi"]] <- multi.res
      } ## End of multi.ppi
      
      if("go.cc.ppi" %in% obj.fun) {
        go.ppi.obj <- mclapply(X = as.list(unique(ppi.obj$n)), FUN = function(x) {
          go.ppi.df <- ppi.obj %>% dplyr::filter(n == x) %>% dplyr::select(-celltype1, -celltype2) %>%
            left_join(x = ., ## Add cellular component info for all genes in gene1 (NA means no info on gene)
                      y = map.go.gene[,c("external_gene_name", "name_1006")], 
                      by = c("gene1"="external_gene_name"), 
                      relationship = "many-to-many") %>%
            dplyr::rename("cc1"="name_1006") %>% ## Rename cc column as cc1 to link to gene1
            left_join(x = ., ## Add cc info for gene2
                      y = map.go.gene[,c("external_gene_name", "name_1006")],
                      by = c("gene2"="external_gene_name"),
                      relationship = "many-to-many") %>%
            dplyr::rename("cc2"="name_1006") %>% ## Rename cc column as cc2 to link to gene2
            .[which(.$cc1 == .$cc2),] ## Filter for cases where gene1 and gene2 are in same cellular compartment
          ## Get all genes from go.ppi.df
          go.genes <- data.frame(gene = unique(c(go.ppi.df$gene1, go.ppi.df$gene2)),
                                 go.cc = 1)
          go.res <- df %>% dplyr::filter(n == x) %>%
            left_join(x = ., y = go.genes, by = "gene") %>%
            mutate(go.cc = replace_na(go.cc, replace = 0)) %>%
            dplyr::summarise(go.cc = sum(go.cc) / length(go.cc)) %>%
            mutate(n = x) %>% .[, c("n", "go.cc")]
        }, mc.cores = n.cores)
        go.ppi.obj <- as.data.frame(do.call(rbind, go.ppi.obj))
        
        obj.fun.res[["go.cc.ppi"]] <- go.ppi.obj
      }
      
    } ## End of broad ppi umbrella
    
    #################################
    ## Promoter Objective Function ##
    if("promoter" %in% obj.fun) {
      promoter.obj <- df %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
        mutate(promoter = replace_na(promoter, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(promoter = sum(promoter)/length(promoter))
      
      obj.fun.res[["promoter"]] <- promoter.obj
    }
    
    #############################
    ## ATAC Objective Function ##
    ## Marker atac peak 
    if("marker.atac" %in% obj.fun) {
      m.atac.obj <- df %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
        mutate(marker.atac = replace_na(marker.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(marker.atac = sum(marker.atac)/length(marker.atac))
      obj.fun.res[["marker.atac"]] <- m.atac.obj
    }
    ## Common atac peak
    if("common.atac" %in% obj.fun) {
      c.atac.obj <- df %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
        mutate(common.atac = replace_na(common.atac, replace = 0)) %>%
        dplyr::group_by(n) %>% dplyr::summarise(common.atac = sum(common.atac)/length(common.atac))
      obj.fun.res[["common.atac"]] <- c.atac.obj
    }
    ########################
    ## lncrna interaction ##
    if("lncrna" %in% obj.fun) {
      lncrna.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with ppi
        match.set <- rna.protein.map %>% left_join(x = ., y = proposal.set, by = c("gene1"="gene"),
                                                   relationship = "many-to-many") %>%
          dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
          left_join(x = ., y = proposal.set, by = c("gene2"="gene"),
                    relationship = "many-to-many") %>%
          dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
          na.omit() %>%
          mutate(n = x)
      }, mc.cores = n.cores)
      lncrna.obj <- as.data.frame(do.call(rbind, lncrna.obj))
      
      all.lncrna.obj <- mclapply(X = as.list(unique(lncrna.obj$n)), FUN = function(x) {
        lncrna.df <- lncrna.obj %>% dplyr::filter(n == x)
        lncrna.genes <- data.frame(gene = unique(c(lncrna.df$gene1, lncrna.df$gene2)),
                                   lncrna = 1)
        lncrna.res <- df %>% dplyr::filter(n == x) %>%
          left_join(x = ., y = lncrna.genes, by = "gene") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0)) %>%
          dplyr::summarise(lncrna = sum(lncrna) / length(lncrna)) %>%
          mutate(n = x) %>% .[, c("n", "lncrna")]
      }, mc.cores = n.cores)
      all.lncrna.obj <- as.data.frame(do.call(rbind, all.lncrna.obj))
      
      if(nrow(all.lncrna.obj) == 0) { ## If no lncrna interaction then do this
        lncrna.res <- data.frame(n = unique(df$n),
                                 lncrna = 0)
      } else { ## If lncrna interaction present
        lncrna.res <- data.frame(n = unique(df$n)) %>%
          left_join(x = ., y = all.lncrna.obj, by = "n") %>%
          mutate(lncrna = replace_na(lncrna, replace = 0))
      }
      
      obj.fun.res[["lncrna"]] <- lncrna.res
    } ## end of lncrna
    
    ######################
    ## tfbs interaction ##
    if("tfbs" %in% obj.fun) {
      
      tfbs.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n) %>%
          left_join(x = ., y = ct_exp, by = c("celltype"="type"), relationship = "many-to-many") %>% ## Add celltype and tf gene info
          left_join(x = ., y = mb.tfbs.disrupt, by = c("gene_name"="geneSymbol", "locus"="LEAD_SNP")) %>%
          dplyr::select(-gene_name) %>%
          dplyr::group_by(locus, celltype, gene) %>%
          mutate(tfbs = replace_na(tfbs, replace = 0)) %>% ## Add 0 where NA is
          dplyr::summarise(tfbs = sum(tfbs)) %>% ## sum all tfbs disrupt hits
          mutate(tfbs = case_when(tfbs > 0 ~ 1, ## boolean value conversion
                                  tfbs == 0 ~ 0),
                 n = x) %>%
          dplyr::select(celltype, locus, gene, n, tfbs) %>%
          dplyr::ungroup() 
        
      }, mc.cores = n.cores)
      tfbs.obj <- as.data.frame(do.call(rbind, tfbs.obj))
      tfbs.obj <- tfbs.obj %>% dplyr::group_by(n) %>% dplyr::summarise(tfbs = sum(tfbs)/length(tfbs))
      obj.fun.res[["tfbs"]] <- tfbs.obj
    } ## end of tfbs
    
    if("tfbs.marker" %in% obj.fun) {
      tfbs.marker.obj <- mclapply(X = as.list(unique(df$n)), FUN = function(x) {
        ## Get info on proposal
        proposal.set <- df %>% dplyr::filter(n == x) %>% dplyr::select(-n)
        ## Match proposal elements with TF expressed data
        proposal.set <- ct_exp %>% mutate(tf.hit = 1) %>% 
          left_join(x = proposal.set, y = ., by = c("gene"="gene_name", "celltype"="type")) %>%
          mutate(tf.hit = replace_na(tf.hit, replace = 0))
        
        ## If there are any TF expressed in correct cell type then add in motifbreakR info
        if(any(proposal.set$tf.hit == 1)) {
          ## What tf.gene are a hit
          tf.genes <- proposal.set %>% dplyr::filter(tf.hit == 1) %>% pull(gene)
          ## What SNPs have motif broken?
          for(tf.gene in tf.genes) {
            tfbs.disrupt <- mb.tfbs.disrupt %>% dplyr::filter(geneSymbol %in% tf.gene)
            ## Now add a point for all lead SNPs with a MB disruption for that tf.gene
            for(tfbs in unique(tfbs.disrupt$LEAD_SNP)) {
              proposal.set[which(proposal.set$locus == tfbs), "tf.hit"] <- 1
            }
          }
        } ## end of IF condition 
        proposal.set <- proposal.set %>% mutate(n = x) %>% ## add back generation info
          dplyr::select(celltype, locus, gene, n, tf.hit) %>%
          dplyr::group_by(n) %>%
          dplyr::summarise(tfbs = sum(tf.hit)/length(tf.hit))
        
      }, mc.cores = n.cores)
      tfbs.marker.obj <- as.data.frame(do.call(rbind, tfbs.marker.obj))
    }
  }
  
  
  #####################################
  ## Fitness mean objective function ##
  
  for(n in 1:length(obj.fun)) {
    if(n == 1) {
      fitness.obj <- obj.fun.res[[obj.fun[n]]]
    } else {
      add_test <- obj.fun.res[[obj.fun[n]]]
      fitness.obj <- left_join(x = fitness.obj, y = add_test, by = "n")
    }
  }
  fitness.obj[is.na(fitness.obj)] <- 0
  fitness.obj$fitness <- apply(X = fitness.obj %>% dplyr::select(-n), MARGIN = 1, FUN = function(x) {
    mean(x)
  })
  
}

################################
## Get objective function hit ##
################################

get.of.hits <- function(data, obj.fun) { ## Data is df of locus, gene and cell type; obj.fun is a vector of OF to use
  obj.fun.hit <- list()
  ########################################
  ## Combined gene marker - CEVM and FM ##
  if("marker.gene" %in% obj.fun) {
    marker.gene <- combined.gene.marker %>% mutate(marker.gene = 1)
    
    marker.gene.obj <- data %>% left_join(x = ., y = marker.gene, by = c("gene", "celltype")) %>%
      mutate(marker.gene = replace_na(marker.gene, replace = 0))
    obj.fun.hit[["marker.gene"]] <- as.data.frame(marker.gene.obj)
  } ## End of combined gene marker
  
  ###########
  ## Magma ##
  if("magma" %in% obj.fun) {
    magma.obj <- data %>% left_join(x = ., y = magma, by = "gene") %>%
      mutate(magma = replace_na(magma, replace = 0))
    obj.fun.hit[["magma"]] <- as.data.frame(magma.obj)
  } ## End of magma
  
  ##########################
  ## Highly mutated genes ##
  if("cancer.gene" %in% obj.fun) {
    hi.mut.obj <- data %>% left_join(x = ., y = breast.cancer.gene, by = "gene") %>%
      mutate(cancer.gene = replace_na(cancer.gene, replace = 0))
    obj.fun.hit[["cancer.gene"]] <- as.data.frame(hi.mut.obj)
  } ## End of hi mut
  
  ############################
  ## PPI Objective Function ##
  if(any(c("ppi", "intra.ppi", "inter.ppi", "marker.ppi", "go.cc.ppi") %in% obj.fun)) {
    match.set <- ppi %>% left_join(x = ., y = data, by = c("gene1"="gene"),
                                   relationship = "many-to-many") %>%
      dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
      left_join(x = ., y = data, by = c("gene2"="gene"),
                relationship = "many-to-many") %>%
      dplyr::rename("celltype2"="celltype", "locus2"="locus") %>%
      na.omit()
    if("intra.ppi" %in% obj.fun) {
      intra.df <-  match.set %>% mutate(intra.commun = case_when(celltype1 == celltype2 ~ 1,
                                                                 celltype1 != celltype2 ~ 0)) %>%
        ## Label each cases as 1 or 0 for having same celltype or not, then filter for 1.
        dplyr:::filter(intra.commun == 1) %>%
        group_by(grp = paste(pmax(gene1, gene2), ## Keep only intra combo (no celltalk combo)
                             pmin(gene1, gene2), sep = "_")) %>%
        dplyr::filter(!grp %in% celltalk_gene$grp) %>%
        ungroup() %>%
        dplyr::select(-grp) %>%
        distinct()
      if(nrow(intra.df) > 0) {
        intra.df <- data.frame(gene = unique(c(intra.df$gene1, intra.df$gene2)),
                               intra.ppi = 1)
        ##
        intra.res <- data %>% left_join(x = ., y = intra.df, by = "gene") %>%
          mutate(intra.ppi = replace_na(intra.ppi, replace = 0))
      } else {
        intra.res <- data %>% mutate(intra.ppi = 0)
      }
      obj.fun.hit[["intra.ppi"]] <- as.data.frame(intra.res)
    } ## End of intra.ppi
    
    ##############
    ## PPI Only ##
    if("ppi" %in% obj.fun) {
      ppi.genes <- data.frame(gene = unique(c(match.set$gene1, match.set$gene2)),
                              ppi = 1)
      ppi.res <- data %>% left_join(x = ., y = ppi.genes, by = "gene") %>%
        mutate(ppi = replace_na(ppi, replace = 0))
      obj.fun.hit[["ppi"]] <- as.data.frame(ppi.res)
    } ## End of PPI 
    
    #####################
    ## marker ppi only ##
    if("marker.ppi" %in% obj.fun) {
      all.markers <- combined.gene.marker %>% mutate(marker.gene = 1) %>% distinct()
      
      marker.obj <- match.set %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
        left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
        na.omit()
      
      if(nrow(marker.obj) == 0) {
        marker.res <- data %>% mutate(marker.ppi = 0)
      } else {
        df1 <- marker.obj %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
        df2 <- marker.obj %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
        marker.ppi.set <- rbind(df1, df2) %>% mutate(marker.ppi = 1) %>% distinct()
        
        marker.res <- data %>% left_join(x = ., y = marker.ppi.set, by = c("gene", "celltype")) %>%
          mutate(marker.ppi = replace_na(marker.ppi, replace = 0))
      }
      obj.fun.hit[["marker.ppi"]] <- as.data.frame(marker.res)
    } ## End of marker ppi
    
    ###############
    ## Inter ppi ##
    if("inter.ppi" %in% obj.fun) {
      all.markers <- combined.gene.marker %>% mutate(marker.gene = 1) %>% distinct()
      marker.obj <- match.set %>% left_join(x = ., y = all.markers, by = c("gene1"="gene", "celltype1"="celltype")) %>%
        left_join(x = ., y = all.markers, by = c("gene2"="gene", "celltype2"="celltype")) %>%
        na.omit() %>% ## Keep positive hits.
        dplyr::select(-marker.gene.x, -marker.gene.y) ## Get rid of columns no longer needed
      if(nrow(marker.obj) == 0) {
        inter.ppi.obj <- data %>% mutate(inter.ppi = 0)
      } else {
        inter.obj <- marker.obj %>% left_join(x = ., 
                                              y = celltalkdb[,c("source_genesymbol", "inter.ppi")], 
                                              by = c("gene1"="source_genesymbol"),
                                              relationship = "many-to-many") %>%
          distinct() %>%
          dplyr::rename("source_gene1" = "inter.ppi") %>%
          left_join(x = ., 
                    y = celltalkdb[,c("source_genesymbol", "inter.ppi")],
                    by = c("gene2"="source_genesymbol"),
                    relationship = "many-to-many") %>%
          distinct() %>% 
          dplyr::rename("source_gene2" = "inter.ppi") %>%
          .[which(rowSums(.[,c("source_gene1", "source_gene2")], na.rm = TRUE) == 1),] %>%
          
          ## Now trying to find targets the same way we did for source genes
          left_join(x = ., 
                    y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                    by = c("gene1"="target_genesymbol"),
                    relationship = "many-to-many") %>%
          distinct() %>%
          dplyr::rename("target_gene1" = "inter.ppi") %>%
          left_join(x = ., 
                    y = celltalkdb[,c("target_genesymbol", "inter.ppi")],
                    by = c("gene2"="target_genesymbol"),
                    relationship = "many-to-many") %>%
          distinct() %>%
          dplyr::rename("target_gene2" = "inter.ppi") %>%
          .[which(rowSums(.[,c("target_gene1", "target_gene2")], na.rm = TRUE) == 1),] %>%
          ## Removing target and source identity column
          dplyr::select(-source_gene1, -source_gene2, -target_gene1, -target_gene2)
        if(nrow(inter.obj) == 0) {
          inter.ppi.obj <- data %>% mutate(inter.ppi = 0)
        } else {
          df1 <- inter.obj %>% dplyr::select(gene1, celltype1) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
          df2 <- inter.obj %>% dplyr::select(gene2, celltype2) %>% 'colnames<-' (c("gene", "celltype")) %>% distinct()
          inter.ppi.set <- rbind(df1, df2) %>% mutate(inter.ppi = 1) %>% distinct()
          
          inter.ppi.obj <- data %>% left_join(x = ., y = inter.ppi.set, by = c("gene", "celltype")) %>%
            mutate(inter.ppi = replace_na(inter.ppi, replace = 0))
        }
      }
      obj.fun.hit[["inter.ppi"]] <- as.data.frame(inter.ppi.obj)
    } ## End of inter.ppi
  } ## End of all PPI function
  
  #################################
  ## Promoter Objective Function ##
  if("promoter" %in% obj.fun) {
    promoter.obj <- data %>% left_join(x = ., y = promoter.df, by = c("locus", "gene")) %>%
      mutate(promoter = replace_na(promoter, replace = 0)) 
    
    obj.fun.hit[["promoter"]] <- promoter.obj
  } ## End of promoter
  
  #############################
  ## ATAC Objective Function ##
  ## Marker atac peak 
  if("marker.atac" %in% obj.fun) {
    m.atac.obj <- data %>% left_join(x = ., y = atac.marker, by = c("locus", "celltype")) %>%
      mutate(marker.atac = replace_na(marker.atac, replace = 0))
    
    obj.fun.hit[["marker.atac"]] <- m.atac.obj
  } ## End of marker atac
  ## Common atac peak
  if("common.atac" %in% obj.fun) {
    c.atac.obj <- data %>% left_join(x = ., y = common.atac.hit, by = c("locus", "celltype")) %>%
      mutate(common.atac = replace_na(common.atac, replace = 0)) 
    obj.fun.hit[["common.atac"]] <- c.atac.obj
  } ## End of common atac
  
  ########################
  ## lncrna interaction ##
  if("lncrna" %in% obj.fun) {
    ## Match proposal elements with ppi
    match.set <- rna.protein.map %>% left_join(x = ., y = data, by = c("gene1"="gene"),
                                               relationship = "many-to-many") %>%
      dplyr::rename("celltype1"="celltype", "locus1"="locus") %>%
      left_join(x = ., y = data, by = c("gene2"="gene"),
                relationship = "many-to-many") %>%
      dplyr::rename("celltype2"="celltype", "locus2"="locus") %>% 
      na.omit()
    
    if(nrow(match.set)==0) { ## If match.set is 0. meaning no lncrna and protein interactions
      lncrna.res <- data %>% mutate(lncrna = 0)
    } else { ## If there is a match then do the following
      lncrna.genes <- data.frame(gene = unique(c(match.set$gene1, match.set$gene2)),
                                 lncrna = 1)
      lncrna.res <- data %>% left_join(x = ., y = lncrna.genes, by = "gene") %>%
        mutate(lncrna = replace_na(lncrna, replace = 0))
    }
    
    obj.fun.hit[["lncrna"]] <- lncrna.res
  } ## end of lncrna
  
  ## Combine the obj.fun.hit together 
  for(n in 1:length(obj.fun)) {
    if(n == 1) {
      obj.hits <- obj.fun.hit[[obj.fun[n]]]
    } else {
      add_test <- obj.fun.hit[[obj.fun[n]]]
      obj.hits <- left_join(x = obj.hits, y = add_test, by = c("locus", "gene", "celltype"))
    }
  }
  
  ## Return final res
  return(obj.hits)
} 

#########################################
## Get consensus for gene and celltype ##
#########################################

## Post GA analysis
## Requires 5 column: "celltype" "locus"    "gene"     "n"        "gen"   
consensus.proposal <- function(top.proposal, loci, nearbygenes) {
  genes <- subset(top.proposal, select = c(locus, gene))
  geneframe <- data.frame(locus=character(length(loci)), gene = character(length(loci)), Freq = integer(length(loci)))
  j <- 1
  for (i in unique(genes$locus)) {
    genetab <- as.data.frame(table(subset(genes, locus==i)))
    topgene <- which(genetab$Freq == max(genetab$Freq))
    if(length(topgene)>1) {
      topgene <- which(genetab$Freq == max(genetab$Freq)) %>% .[1]
    }
    geneframe[j,1] <- as.character(genetab[topgene,1])
    geneframe[j,2] <- as.character(genetab[topgene,2])
    geneframe[j,3] <- genetab[topgene,3]
    j <- j+1
  }
  
  cells <- subset(top.proposal, select = c(locus, celltype))
  cellframe <- data.frame(locus=character(length(loci)), celltype = character(length(loci)), Freq = integer(length(loci)))
  j <- 1
  for (i in unique(cells$locus)) {
    celltab <- as.data.frame(table(subset(cells, locus == i)))
    topcell <- which(celltab$Freq == max(celltab$Freq))
    if(length(topcell)>1) {
      topcell <- which(celltab$Freq == max(celltab$Freq)) %>% .[1]
    }
    cellframe[j,1] <- as.character(celltab[topcell,1])
    cellframe[j,2] <- as.character(celltab[topcell,2])
    cellframe[j,3] <- celltab[topcell,3]
    j <- j+1
  }
  
  result <- merge(geneframe, cellframe, by = "locus") %>%
    left_join(x = ., y = distinct(nearbygenes[, c("gene_name", "gene_type")]), by = c("gene"="gene_name"))
}

#####################################
## Get max hit from scaling gp run ##

get.scale.hit <- function(data, gen, n_loci) {
  ## Pull result from multi.gp
  plotdata.res <- list()
  
  for(n in 1:11) {
    ## Loop through each OF and process
    set1 <- multi.gp[[n]]
    set1 <- as.data.frame(do.call(rbind, set1[1:200]))
    set1$of <- names(set1)[2]
    names(set1) <- c("n", "of_score", "fitness", "gen", "of")
    ## Store res in list
    plotdata.res[[n]] <- set1
  }
  ## List to dataframe
  plotdata.res <- as.data.frame(do.call(rbind, plotdata.res))
  
  ## Get max hit
  max.hit <- plotdata.res %>% dplyr::filter(gen == gen) %>%
    dplyr::group_by(of) %>%
    dplyr::summarise(max = max(of_score)) %>%
    dplyr::ungroup() %>%
    mutate(max_loci = max * n_loci)
  
  return(max.hit)
}

######################################
## Get consensus gene and cell type ##

get.consensus.run <- function(data, loci, nearbygenes) {
  genes <- subset(data, select = c(locus, gene))
  geneframe <- data.frame(locus=character(length(loci)), gene = character(length(loci)), Freq = integer(length(loci)))
  
  j <- 1
  for (i in unique(genes$locus)) {
    genetab <- as.data.frame(table(subset(genes, locus==i)))
    topgene <- which(genetab$Freq == max(genetab$Freq))
    if(length(topgene)>1) {
      topgene <- which(genetab$Freq == max(genetab$Freq)) %>% .[1]
    }
    geneframe[j,1] <- as.character(genetab[topgene,1])
    geneframe[j,2] <- as.character(genetab[topgene,2])
    geneframe[j,3] <- genetab[topgene,3]
    j <- j+1
  }
  
  cells <- subset(data, select = c(locus, celltype))
  cellframe <- data.frame(locus=character(length(loci)), celltype = character(length(loci)), Freq = integer(length(loci)))
  
  j <- 1
  for (i in unique(cells$locus)) {
    celltab <- as.data.frame(table(subset(cells, locus == i)))
    topcell <- which(celltab$Freq == max(celltab$Freq))
    if(length(topcell)>1) {
      topcell <- which(celltab$Freq == max(celltab$Freq)) %>% .[1]
    }
    cellframe[j,1] <- as.character(celltab[topcell,1])
    cellframe[j,2] <- as.character(celltab[topcell,2])
    cellframe[j,3] <- celltab[topcell,3]
    j <- j+1
  }
  
  result <- merge(geneframe, cellframe, by = "locus") %>%
    left_join(x = ., y = distinct(nearbygenes[, c("gene_name", "gene_type")]), by = c("gene"="gene_name"))
}

###############################
## get objective function OF ##
## function to pull out no omics obj fun score each generation 
getOFmean <- function(data, type, omic) {
  df <- as.data.frame(do.call(rbind, data[1:100])) %>%
    dplyr::group_by(gen) %>%
    dplyr::summarise(marker.gene = mean(marker.gene),
                     magma = mean(magma),
                     cancer.gene = mean(cancer.gene),
                     ppi = mean(ppi),
                     inter.ppi = mean(inter.ppi),
                     intra.ppi = mean(intra.ppi),
                     lncrna = mean(lncrna),
                     promoter = mean(promoter),
                     marker.atac = mean(marker.atac),
                     common.atac = mean(common.atac),
                     fitness = mean(fitness)) %>%
    mutate(type = type,
           omic = omic)
  return(df)
}

