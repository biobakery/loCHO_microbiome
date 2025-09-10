#script that contains functions to load the various cleaned data types for analysis

#Function that loads the taxonomic data for LowCHO data

library(readxl)
library(mgsub)
library(dplyr)
library(stringr)

load_taxa_data <- function(level="SGB", keep_unknown=TRUE){
  
  m4_data_file <- "Data/Taxonomic_profiles/metaphlan_4_0_6_Oct22_fixed_profiles.tsv"
  
  bug4 <- read.csv(file = m4_data_file, sep = '\t', header = T, check.names = F)
  
  #remove F-50179105-003S_S46
  row.names(bug4) = bug4$`# taxonomy`
  bug4$`# taxonomy` = NULL
  
  if(level=="SGB"){
    bug = keep_SGB(bug4)
  }else{
    bug = keep_otherLevel(bug4, level=level)
  }
  
  bug = data.frame(t(bug), check.names = F)
  bug4 = data.frame(t(bug4), check.names = F)
  
  if(keep_unknown){
    bug$UNKNOWN = bug4$UNCLASSIFIED
  }

  return(bug)
}

load_meta_data <- function(m4_data=bug){
  
  metadata_file <- "Data/Metadata/Hills Jackson LoCHO MS_Metagenomics and Metadata.xlsx - Hills Jackson LoCHO MS_all data.csv"
  metadata <- read.csv(file = metadata_file, header = T, check.names = F)

  row.names(metadata) = metadata$`Metagenomics Sample Number`
  metadata = metadata[match(row.names(bug), row.names(metadata)), ]
  
  
  metadata$`Diet Name` = mgsub(metadata$`Diet Name`, c("LoCHO_PROT", "LoCHO_FAT", "HiCHO"), c("High Protein", "High Fat", "High Carb"))
  metadata$`Diet Order Names` = mgsub(metadata$`Diet Order Names`, c("LoCHO_PROT - LoCHO_FAT", "LoCHO_FAT - LoCHO_PROT"), c("High Protein --> High Fat", "High Fat --> High Protein"))
  
  #load in read count file
  ##load in read count table
  read_counts <- read.table("Data/Taxonomic_profiles/kneaddata_read_count_table.tsv",
                            sep="\t", header=T)
  
  read_counts$Sample <- gsub("-003.*", "-003S", read_counts$Sample)
  colnames(read_counts)[1] <- "Metagenomics Sample Number"
  read_counts <- read_counts %>% filter(`Metagenomics Sample Number` %in% metadata$`Metagenomics Sample Number`)
  read_counts <- read_counts %>% mutate(total_counts=final.pair1 + final.pair2 + final.orphan1 + final.orphan2)
  
  metadata <- metadata %>% left_join(read_counts)
  metadata <- data.frame(metadata, check.names=F)
  rownames(metadata) <- metadata$Metagenomics.Sample.Number
  return(metadata)
}


load_EC_data <- function(metadata=metadata, strat=FALSE, remove_ungroup=FALSE){
  ec_file <- "Data/Funcational_profiles/merged_ecs.tsv"
  ecs <- read.csv(file=ec_file,
    sep="\t", header=T, row.names=1,
    check.names = F)
  
  names(ecs) <- gsub("_.*", "", names(ecs))
  
  if(strat){
    filt_ecs <- keep_unstrat(ecs, rev=TRUE)
  }else{
    filt_ecs <- keep_unstrat(ecs)
  }
  if(remove_ungroup){
    filt_ecs <- filt_ecs[!grepl("UNGROUP*", rownames(filt_ecs)),]
  }
  #make sure its the same order as the metaedata
  filt_ecs <- filt_ecs[,rownames(metadata)]
  #flip it so ECs are the columns
  filt_ecs <- data.frame(t(filt_ecs), check.names = F)
  return(filt_ecs)
}


load_geneFam_data <- function(metadata=metadata, strat=FALSE, remove_unmapped=FALSE){
  gene_file <- "Data/Funcational_profiles/merged_genefamilies.tsv"
  
  genes <- read.csv(file=gene_file,
    sep="\t", header=T, row.names=1,
    check.names = F)
  
  names(genes) <- gsub("_.*", "", names(genes))
  
  if(strat){
    filt_genes <- keep_unstrat(genes, rev = TRUE)
  }else{
    filt_genes <- keep_unstrat(genes)
  }
  
  if(remove_unmapped){
    filt_genes <- filt_genes[!grepl("READS_UNMAPPED", rownames(filt_genes)),]
  }

  filt_genes <- filt_genes[,rownames(metadata)]
  filt_genes <- data.frame(t(filt_genes), check.names = F)
  return(filt_genes)
  
}


load_pathway_data <- function(metadata=metadata, strat=FALSE, remove_unknown=FALSE){
  pathway_file <- "Data/Funcational_profiles/merged_pathways.tsv"
  pathways <- read.csv(file=pathway_file, 
    sep="\t", header=T, row.names=1,
    check.names = F)
  
  names(pathways) <- gsub("_.*", "", names(pathways))
  
  if(strat){
    filt_path <- keep_unstrat(pathways, rev = TRUE)
  }else{
    filt_path <- keep_unstrat(pathways)
  }
  
  if(remove_unknown){
    filt_path <- filt_path[!grepl("UNMAPPED", rownames(filt_path)),]
    filt_path <- filt_path[!grepl("UNINTEGRATED", rownames(filt_path)),]
  }
  
  filt_path <- filt_path[,rownames(metadata)]
  filt_path <- data.frame(t(filt_path), check.names = F)
  return(filt_path)
}

load_pfam_data <- function(metadata=metadata, strat=FALSE, remove_unknown=FALSE){
  pfam_file <- "Data/Funcational_profiles/merged_pfams.tsv"
  
  pfams <- read.csv(file=pfam_file, sep="\t", header=T, row.names=1, check.names=F)
  
  names(pfams) <- gsub("_.*", "", names(pfams))
  
  if(strat){
    filt_pfams <- keep_unstrat(pfams, rev=T)
  }else{
    filt_pfams <- keep_unstrat(pfams)
  }
  
  if(remove_unknown){
    filt_pfams <- filt_pfams[!grepl("UNGROUPED", rownames(filt_pfams)),]
  }
  
  filt_pfams <- filt_pfams[,rownames(metadata)]
  filt_pfams <- data.frame(t(filt_pfams), check.names = F)
  return(filt_pfams)
}


#Load fecalMBX data directly from metabolon datasheet
load_fecalMbx_data <- function(metadata){
  mbx_file <- "Data/Metabolomic_profiles//HILL-0105-20MDTA CLIENT DATA TABLE (FECES).XLSX"
  
  fecal_mbx <- read_excel(mbx_file, sheet=4)  
  
  colnames(fecal_mbx) <- gsub("-002S", "-003S", colnames(fecal_mbx))
  
  fecal_mbx <- data.frame(t(fecal_mbx))
  colnames(fecal_mbx) <- unlist(fecal_mbx[1,])
  
  fecal_mbx <- fecal_mbx[-1,]
  
  fecal_mbx <- fecal_mbx[rownames(metadata),]
  
  fecal_mbx <- fecal_mbx %>% mutate_if(is.character, as.numeric)
  
  return(fecal_mbx)
  
}

Impute_MBX_NAs <- function(x){
  ### takes a row
  min_value <- min(x, na.rm = T)
  x[is.na(x)] <- min_value
  return(x)
}


load_rawfecalMbx_data <- function(metadata){
  mbx_file <- "Data/Metabolomic_profiles/HILL-0105-20MDTA CLIENT DATA TABLE (FECES).XLSX"
  
  fecal_mbx <- read_excel(mbx_file, sheet=2)
  fecal_mbx <- data.frame(fecal_mbx[-c(1:13),],check.names = F, check.rows = F)
  rownames(fecal_mbx) <- fecal_mbx[,2]
  fecal_mbx <- fecal_mbx[,-c(1:13)]
  fecal_mbx <- fecal_mbx[-1,]
  
  
  colnames(fecal_mbx) <- gsub("-002S", "-003S", colnames(fecal_mbx))
  fecal_mbx <- data.frame(t(fecal_mbx),check.names = F)
  row_names <- rownames(fecal_mbx)
  #convert to numeric
  fecal_mbx <- apply(fecal_mbx, 2, as.numeric)
  
  ## if something is not detected we will replace it by the minimum/2
  fecal_mbx <- apply(fecal_mbx, 2, Impute_MBX_NAs)
  rownames(fecal_mbx) <- row_names
  fecal_mbx <- fecal_mbx[rownames(metadata),]
  
  return(fecal_mbx)
}

#Load fecalMBX data directly from metabolon datasheet
load_serumMbx_data <- function(metadata){
  mbx_file <- "Data/Metabolomic_profiles/HILL-0105-20MDTA CLIENT DATA TABLE (SERUM).XLSX"
  
  serum_mbx <- read_excel(mbx_file, sheet=4)  
  
  colnames(serum_mbx) <- gsub("-002S", "-003S", colnames(serum_mbx))
  
  serum_mbx <- data.frame(t(serum_mbx))
  colnames(serum_mbx) <- unlist(serum_mbx[1,])
  
  serum_mbx <- serum_mbx[-1,]
  
  serum_mbx <- serum_mbx[rownames(metadata),]
  
  serum_mbx <- serum_mbx %>% mutate_if(is.character, as.numeric)
  
  return(serum_mbx)
  
}

load_rawserumMbx_data <- function(metadata){
  mbx_file <- "Data/Metabolomic_profiles/HILL-0105-20MDTA CLIENT DATA TABLE (SERUM).XLSX"
  
  serum_mbx <- read_excel(mbx_file, sheet=2)  
  
  serum_sample_ids <- data.frame(serum_mbx[12:14,-c(1:12)], check.names = F)
  rownames(serum_sample_ids) <- serum_sample_ids[,1]
  
  
  mbx_file_fecal <- "Data/Metabolomic_profiles/HILL-0105-20MDTA CLIENT DATA TABLE (FECES).XLSX"
  
  fecal_mbx <- read_excel(mbx_file_fecal, sheet=2)
  fecal_sample_ids <- data.frame(fecal_mbx[11:13,-c(1:12)], check.names = F)
  rownames(fecal_sample_ids) <- serum_sample_ids[,1]
  
  
  matching_order <- apply(serum_sample_ids[1:3, ], 2, function(col) {
    # Match based on the content in the first 3 rows
    apply(fecal_sample_ids[1:3, ], 2, function(reference_col) identical(col, reference_col))
  })
  col_indices <- max.col(matching_order, ties.method = "first")
  serum_sample_ids <- serum_sample_ids[, col_indices]
  identical(unname(serum_sample_ids), unname(fecal_sample_ids))
  
  mapping_fecal_to_serum <- colnames(fecal_sample_ids)
  names(mapping_fecal_to_serum) <- colnames(serum_sample_ids)
  ##okay so now we have the matching we can map them...
  
  
  serum_mbx <- data.frame(serum_mbx[-c(1:14),],check.names = F, check.rows = F)
  
  
  
  rownames(serum_mbx) <- serum_mbx[,2]
  serum_mbx <- serum_mbx[,-c(1:13)]
  serum_mbx <- serum_mbx[-1,]
  
  sample_match_indexs <- match(colnames(serum_mbx), names(mapping_fecal_to_serum))
  colnames(serum_mbx) <- mapping_fecal_to_serum[sample_match_indexs]
  
  row_names <- rownames(serum_mbx)
  #convert to numeric
  serum_mbx <- apply(serum_mbx, 2, as.numeric)
  
  ## if something is not detected we will replace it by the minimum/2
  serum_mbx <- apply(serum_mbx, 2, Impute_MBX_NAs)
  rownames(serum_mbx) <- row_names
  colnames(serum_mbx) <- gsub("-002S", "-003S", colnames(serum_mbx))
  serum_mbx <- serum_mbx[,rownames(metadata)]
  serum_mbx <- data.frame(t(serum_mbx),check.names = F)
  return(serum_mbx)
}


Annotate_MBX <- function(mbx){
  
  #load annotation file
  annotation_file <- "Data/Metabolomic_profiles/HMDB_chemical_taxonomy.csv"
  name_files <- "Data/Metabolomic_profiles/Keto Microbiome Publication- Metabolite IDs.xlsx"
  
  mbx_annos <- read.csv(annotation_file, header=T, check.names=F)
  name_files <- read_excel(name_files)
  
  match(name_files$HMDB, mbx_annos$HMDB.ID)
  
  melt_mbx <- melt(as.matrix(mbx))
  
}
