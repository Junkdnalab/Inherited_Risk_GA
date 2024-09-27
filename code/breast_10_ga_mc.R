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
  library(Seurat)
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


message("Loading Data")

load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast_sct_marker.list_pairwise.rda")
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast_nearbygenes_200kb_notelmer_10xgenomic_with_LD.rda") ## candidates near snps
nearbygenes <- nearbygenes_10xgenomics; rm(nearbygenes_10xgenomics)
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/gwas_snps.rda")  ## gwas snp
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/marker.list.bg.rda"); rm(marker.list); gc()
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_promoter.df.rda")
## Clean up promoter data for OF
promoter.df <- promoter.df %>% mutate(locus = LEAD_SNP) %>%
  dplyr::select(-LEAD_SNP) %>%
  distinct()
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast_cevm.rda")
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_atac.marker.rda")
## Clean up markerATAC for OF
atac.marker <- atac.marker %>% dplyr::rename("marker.atac"="atac") %>%
  mutate(locus = LEAD_SNP) %>%
  dplyr::select(-LEAD_SNP) %>%
  distinct()
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/breast.cancer.gene.rda")
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/fusion.gene.rda")
## Merge fusion.gene with breast.cancer.gene
fusion.gene <- fusion.gene %>% dplyr::rename("cancer.gene"="fusion") ## Rename colname for rbind
breast.cancer.gene <- rbind(breast.cancer.gene, fusion.gene) ## rbind together
breast.cancer.gene <- distinct(breast.cancer.gene); rm(fusion.gene) ## Keep unique observation
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/celltalkdb.rda")
celltalkdb$inter.ppi <- 1 ## Add inter.ppi column for OF
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_common_atac_peak.rda")
## Clean up commonATAC for OF
common.atac.hit <- common.atac.hit %>% dplyr::select(locus, GroupReplicate, LEAD_SNP) %>%
  dplyr::rename("celltype"="GroupReplicate") %>%
  mutate(common.atac = 1) %>%
  distinct() %>%
  mutate(locus = LEAD_SNP) %>%
  dplyr::select(-LEAD_SNP) %>%
  distinct()

## Get map of loci to gene
loci2gene <- nearbygenes %>% dplyr::select(snp, gene_name) %>% 
  remove_rownames() %>%
  dplyr::rename("cytoband"="snp") %>%
  distinct()
## Get unique loci
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

## Select PPI for genes found in nearbygene set and ppi.threshold
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
  filter(final_score > ppi.threshold) %>% ## ppi.threshold filtering
  dplyr::select(gene1, gene2);# rm(str.db.info, str.db.link); gc()

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

## Combine findmarkers and cevm
combined.gene.marker <- rbind(cevm.genes, marker.hits) %>% distinct()

## Load lncRNA interaction
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/rna.protein.map.rda")
rna.protein.map <- rna.protein.map %>% 'colnames<-' (c("gene1", "gene2")) %>%
  group_by(grp = paste(pmax(gene1, gene2), pmin(gene1, gene2), sep = "_")) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  dplyr::select(-grp) %>%
  distinct()

# save(atac.marker, breast.cancer.gene, breast.htftarget, combined.gene.marker, celltalkdb, celltalk_gene, common.atac.hit, ct_exp, exp.gene,
#      intra.ppi.gene, loci2gene, magma, marker.hits, mb.remap.df, mb2p.df, mb2p.tf, nearbygenes, ppi, promoter.df, rna.protein.map, max.hit,
#      file = "/drive-pool/data/peter_data/repo/ga_gwas_manu/data/breast_essential_092424.rda")

###########################
## Generate 1k proposals ## Creates 10 sets of 1,000 proposals
message("Generating proposals")
mp_list <- list() 

for(n in 1:10) {
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

#########################
#########################
## Genetic programming ##
#########################
#########################

message("Running GA Model")

## Objective function
obj.fun <- c("marker.gene", "magma", "cancer.gene", "ppi", "marker.ppi", "inter.ppi", "intra.ppi",
             "lncrna", "promoter","marker.atac", "common.atac")

n.run <- data.frame(run = 1:10) %>% t() %>% as.data.frame() %>% as.list() ## HOw many runs are we doing? 10

## Code to run GA model after loading all files. Code uses at most 100 cores. Adjust to your core size
mclapply(X = n.run, FUN = function(x) {
  num.run <- as.numeric(unlist(x))
  
  multi_proposal <- mp_list[[num.run]]

  multi.gp <- gwas_gp(init.proposal = multi_proposal, generations = 200, obj.fun = obj.fun, khan.method = TRUE, selection.method = "progressive", parent.size = 1e2, proposal.size = 1e3, n.cores = 10, mut.rate = 0.01)
  save(multi.gp, file = paste0("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/data/progressive_breast_run_", num.run , "_ppithres", ppi.threshold, "_092324.rda"))
}, mc.cores = 10)

