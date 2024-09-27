#######################
#######################
## Loading libraries ##
#######################
#######################

message("Loading libraries")

suppressMessages(c(
  library(tidyverse),
  library(data.table),
  library(parallel),
  library(GenomicRanges)
))

##########################################
##########################################
## Load functions for genetic algorithm ##
##########################################
##########################################

message("Loading functions")

source("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/code/gp_fun_update.R")

###############
###############
## Load data ##
###############
###############

## Objective function
obj.fun <- c("marker.gene", "magma", "cancer.gene", "ppi", "marker.ppi", "inter.ppi", "intra.ppi",
             "lncrna", "promoter","marker.atac", "common.atac")

load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/vs_controls_051524.rda")
## Get control group id
control_group <- unique(vs_controls$LABEL); rm(vs_controls)


for(control in control_group) { ## Will loop through each control set and perform GA. Each control set has diff SNP location. Needs diff data nearbygenes, promoter, etc.
  
  ###############
  ## Load data ##
  ###############
  
  ## SNP data
  load(paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/control/",control,"_snp.rda"))
  
  data <- data %>% mutate(CHR_ID = paste0("chr", CHR_ID),
           SNP_ID_CURRENT = paste0(CHR_ID, ":", CHR_POS)) %>%
    distinct()
  ## Promoter data
  load(paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/control/",control,"_ld_promoter.df.rda"))
  promoter.df <- promoter.df %>% mutate(locus = LEAD_SNP) %>% ## This was changed with new ctrl
    dplyr::select(-LEAD_SNP) %>%
    distinct()
  ## ATAC marker data
  load(paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/control/",control,"_ld_atac.marker.rda"))
  atac.marker <- atac.marker %>% dplyr::rename("marker.atac"="atac") %>%
    mutate(locus = LEAD_SNP) %>%
    dplyr::select(-LEAD_SNP) %>%
    distinct()
  ## ATAC common data
  load(paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/control/",control,"_ld_atac.common.rda"))
  common.atac.hit <- common.atac.hit %>% mutate(locus = LEAD_SNP) %>%
    dplyr::select(locus, GroupReplicate) %>% 
    dplyr::rename("celltype"="GroupReplicate") %>%
    mutate(common.atac = 1) %>%
    distinct()
  ## Nearbygene
  load(paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/control/",control,"_nearbygenes.rda"))
  nearbygenes <- nearbygenes_10xgenomics; rm(nearbygenes_10xgenomics)
  
  
  ## Other data not control specific 
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast_sct_marker.list_pairwise.rda")
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/marker.list.bg.rda"); rm(marker.list); gc()
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast_cevm.rda")
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast.cancer.gene.rda")
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/fusion.gene.rda")
  ## Merge fusion.gene with breast.cancer.gene
  fusion.gene <- fusion.gene %>% dplyr::rename("cancer.gene"="fusion") ## Rename colname for rbind
  breast.cancer.gene <- rbind(breast.cancer.gene, fusion.gene) ## rbind together
  breast.cancer.gene <- distinct(breast.cancer.gene); rm(fusion.gene) ## Keep unique observation
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/celltalkdb.rda")
  celltalkdb$inter.ppi <- 1 ## Add inter.ppi column
  
  ## Clean nearbygenes to get only map snp to gene
  loci2gene <- nearbygenes %>% dplyr::select(snp, gene_name) %>% 
    remove_rownames() %>%
    dplyr::rename("cytoband"="snp") %>%
    distinct()
  ## unique locus name
  loci <- unique(loci2gene$cytoband) 
  
  ############
  ## PPI df ##
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ppi.map.rda"); rm(map.genes)

  ppi.threshold <- 0.0 ## Confidence score for ppi interaction by stringdb. Medium confid 0.4; High confid 0.7; highest confid 0.9
  
  ## Update the link data to have gene symbol
  str.db.link <- str.db.link %>% dplyr::select(protein1, protein2, experiments, 
                                               experiments_transferred, combined_score) %>%
    left_join(x = ., y = str.db.info, by = c("protein1"="string_id")) %>%
    dplyr::rename("gene1"="gene_symbol") %>%
    left_join(x = ., y = str.db.info, by = c("protein2"="string_id")) %>%
    dplyr::rename("gene2"="gene_symbol") %>%
    .[,c("protein1", "protein2", "gene1", "gene2", "experiments", "experiments_transferred", "combined_score")]
  
  ## Select PPI based on nearby gene set and ppi threshold
  ppi <- str.db.link %>% dplyr::filter(gene1 %in% nearbygenes$gene_name, 
                                       gene2 %in% nearbygenes$gene_name) %>% ## Filter for specific genes found in gwas 200kb
    group_by(grp = paste(pmax(gene1, gene2), pmin(gene1, gene2), sep = "_")) %>%
    dplyr::filter(row_number() == 1) %>%
    ungroup() %>%
    dplyr::select(-grp) %>%
    distinct() %>%
    mutate(experiments_score = compute_prior_away(score = experiments),
           experiments_transferred_score = compute_prior_away(score = experiments_transferred)) %>% 
    mutate(experiments_both_prior_corrected = 1 - (1 - experiments_score) * (1 - experiments_transferred_score)) %>% 
    mutate(final_score = (experiments_both_prior_corrected * (1 - 0.041)) + 0.041) %>%
    filter(final_score > ppi.threshold) %>% ## Filter PPI based on ppi.thres
    dplyr::select(gene1, gene2); rm(str.db.info, str.db.link); gc()
  
  ## Prepare intra.cellular ppi
  ## Remove gene combo found in celltalkdb
  celltalk_gene <- celltalkdb %>% dplyr::select(source_genesymbol, target_genesymbol) %>%
    dplyr::rename("gene1"="source_genesymbol", "gene2"="target_genesymbol") %>%
    group_by(grp = paste(pmax(gene1, gene2), 
                         pmin(gene1, gene2), sep = "_"))
  
  intra.ppi.gene <- ppi %>% 
    group_by(grp = paste(pmax(gene1, gene2), 
                         pmin(gene1, gene2), sep = "_")) %>%
    dplyr::filter(!grp %in% celltalk_gene$grp) %>%
    dplyr::filter(row_number() == 1) %>%
    ungroup() %>%
    dplyr::select(-grp) %>%
    distinct()
  
  ######################
  ## Load Magma genes ##
  magma <- fread("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/UKB_Breast_cancer.gene_score.txt") %>%
    'colnames<-' ("n") %>%
    separate(n, sep = ":", c("gene", "Mscore")) %>% 
    mutate(magma = 1) %>% dplyr::select(-Mscore)
  
  ##################################
  ## Combine findmarkers and cevm ##
  combined.gene.marker <- rbind(cevm.genes, marker.hits) %>% distinct()
  
  #############################
  ## Load lncRNA interaction ##
  load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/rna.protein.map.rda")
  rna.protein.map <- rna.protein.map %>% 'colnames<-' (c("gene1", "gene2")) %>%
    group_by(grp = paste(pmax(gene1, gene2), pmin(gene1, gene2), sep = "_")) %>%
    dplyr::filter(row_number() == 1) %>%
    ungroup() %>%
    dplyr::select(-grp) %>%
    distinct()
  
  ###########################
  ## Generate 1k proposals ##
  mp_list <- list() 
  
  for(n in 1:10) { ## Generate 10 different controls using the loci
    celltype <- as.character(unique(sct_rna$CellType))
    n_proposal <- 1e3
    multi_proposal <- mclapply(X = as.list(1:n_proposal), FUN = function(p) {
      generate_proposal() %>% mutate(n = p)
    },
    mc.cores = 10)
    
    multi_proposal <- as.data.frame(do.call(rbind, multi_proposal))
    mp_list[[n]] <- multi_proposal
  }
  
  rm(multi_proposal)
  
  ##############################################################################################
  message("Running GA Model")
  
  n.run <- data.frame(run = 1:10) %>% t() %>% as.data.frame() %>% as.list() ## How many controls do we have
  ## Will run through process the 10 controls for each control set. Code uses at most 100 cores. Adjust to your core size
  mclapply(X = n.run, FUN = function(x) {
    num.run <- as.numeric(unlist(x))
    
    multi_proposal <- mp_list[[num.run]] ## Select the control replicates
    
    ##################
    ## Run GA model ##
    multi.gp <- gwas_gp(init.proposal = multi_proposal, generations = 200, obj.fun = obj.fun, khan.method = TRUE, selection.method = "progressive", parent.size = 1e2, proposal.size = 1e3, n.cores = 10, mut.rate = 0.01)
    save(multi.gp, file = paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/data/progressive_", control, "_run_", num.run , "_ppithres", ppi.threshold, "_092324.rda"))
  }, mc.cores = 10)
  date()

  ##############################################################################################
  
  } ## GA for loop end







