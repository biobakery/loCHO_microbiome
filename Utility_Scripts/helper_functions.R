### LowCHO utility function

library(stringr)

#filter a feature table by prevelance and abundance
is_common <- function(otu_df,cutoff=0.1, abundance=0.0001){
  num_row = dim(otu_df)[1]
  num_col = dim(otu_df)[2]
  otu_barcode = otu_df > abundance
  common_index = apply(otu_barcode,1,mean) > cutoff
  return(common_index)
}

#Take a M4 formatted table and filter it to the SGB level
keep_SGB <-function(dat_taxa){
  
  temp <- dat_taxa[grepl("t__", rownames(dat_taxa)),]
  SGB_name <- rownames(temp) 
  #get suffix
  SGB_name_suffix <- gsub(".*t__", "", SGB_name)
  #get prefix
  SGB_name_prefix <- sapply(SGB_name, 
                            function(x) str_extract(x,
                                                    "^(.*?)(?=(\\|g__GGB|\\|o__OFGB|\\|c__CFGB|\\|p__PFGB|\\|t__))"))
  #fix prefix to the last known assignment
  SGB_name_prefix <- gsub('.*[kpcofgst]__', "", SGB_name_prefix)
  #combined
  SGB_comb <- paste(SGB_name_prefix, SGB_name_suffix, sep="_")
  rownames(temp) <- SGB_comb
  return(temp)
}

keep_otherLevel <- function(dat_taxa, level){
  
  dat_taxa <- data.frame(dat_taxa, check.names=F)
  dat_taxa <- data.frame(t(dat_taxa), check.names = F)
  names(dat_taxa) = gsub(".*\\|", "", names(dat_taxa))
  dat_taxa <- data.frame(t(dat_taxa), check.names = F)
  
  if(level=="Phylum"){
    ret_frame = dat_taxa[grep("p__.*", row.names(dat_taxa), value = T), ]
  }else if(level=="Class"){
    ret_frame = dat_taxa[grep("c__.*", row.names(dat_taxa), value = T), ]
  }else if(level=="Order"){
    ret_frame = dat_taxa[grep("o__.*", row.names(dat_taxa), value = T), ]
  }else if(level=="Family"){
    ret_frame = dat_taxa[grep("f__.*", row.names(dat_taxa), value = T), ]
  }else if(level=="Genus"){
    ret_frame = dat_taxa[grep("g__.*", row.names(dat_taxa), value = T), ]
  }else if(level=="Species"){
    ret_frame = dat_taxa[grep("s__.*", row.names(dat_taxa), value = T), ]
  }
  return(ret_frame)
}

keep_Phylum <- function(dat_taxa){
  dat_taxa <- data.frame(dat_taxa, check.names=F)
  dat_taxa <- data.frame(t(dat_taxa), check.names = F)
  names(dat_taxa) = gsub(".*\\|", "", names(dat_taxa))
  dat_taxa <- data.frame(t(dat_taxa), check.names = F)
  phylum = dat_taxa[grep("p__.*", row.names(dat_taxa), value = T), ]
}

#Keeps unstratified pathway data
keep_unstrat <-function(dat_path, rev=FALSE){
  #if rev then keep strat...
  if(rev){
    temp = dat_path[grepl("\\|",rownames(dat_path)),]
  }else{
    temp = dat_path[!grepl("\\|",rownames(dat_path)),]
  }
  return(temp)
}


#manual download missed annotations...
# function if HMDB Accession
ChemTax_HMDB <- function(h){
  message(h)
  hmdb_url <- paste0("https://hmdb.ca/metabolites/",h,".xml")
  if(length(grep("Location", curlGetHeaders(hmdb_url))) > 0){
    h_new <- stringr::str_extract(grep("Location", curlGetHeaders(hmdb_url), value = TRUE), "HMDB[0-9]*")
    hmdb_url <- paste0("https://hmdb.ca/metabolites/",h_new,".xml")
  }
  Sys.sleep(sample(5,1)*0.1)
  if(RCurl::url.exists(hmdb_url)){
    hmdb_page <- xml2::read_xml(hmdb_url) 
    node_class <- xml2::xml_text(xml2::xml_find_all(hmdb_page, "//class"))
    node_subclass <- xml2::xml_text(xml2::xml_find_all(hmdb_page, "//sub_class"))
    cbind(h, node_subclass, node_class)
  }
}



