# Emerald Bender 2026
# prepares tree and data files for analysis in BayesTraits


# define function 'bt_treeprep'
# tree.file accepts a string naming tree file in external directory OR a phylo object from env.
# data accepts a data table from env with binomials as row names, either in df or matrix format
# here, data has already had add_taxonomy modify the row names, IE the row names are up to date on the input data object
# shortBranch.thresh set by default to branches shorter than 50k years collapsed
# drop.tips allows manual specification of a vector of otu names to drop from tree
# synonym.table allows inclusion of 2-column synonym table for replacing names in phylogeny and data
#       with original name in column 1, new name in column 2
#
# output is a treedata object (geiger) and a dataframe summarizing why dropped data points were removed
# function removes duplicated data rows, reduces excessively short branches to length 0,
#   collapses zero length branches to polytomies, and matches names in phylogeny to OTL database
# this assumes that data row names have already been cleaned with rotl
#
# duplicates.script / taxonomy.script let you point at your own copies of these shared
# helper scripts; by default they're located relative to this file (two levels up, in
# the Scripts/ directory shared by all wrapper projects), so no machine-specific path
# needs to be hardcoded here.

# locate this file's own directory even when sourced (not just run interactively)
.bt_treeprep_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) NA_character_
)
.bt_helpers_default_dir <- if (!is.na(.bt_treeprep_dir)) dirname(dirname(.bt_treeprep_dir)) else NA_character_

bt_treeprep<- function(tree.file, data, synonym.table = NULL, shortBranch.thresh = 0.05,
                       drop.tips = NULL, cont_colnum= NULL, filestring = NULL,
                       write_predictors = TRUE, scale = TRUE, 
                       duplicates.script = file.path(.bt_helpers_default_dir, "getDuplicates_cladewise.R"),
                       taxonomy.script = file.path(.bt_helpers_default_dir, "add_taxonomy.R")) {

  library(dplyr)
  library(geiger)
  library(rotl)
  library(janitor)
  library(phytools)
  library(tibble)
  source(duplicates.script)
  source(taxonomy.script)

  dsts<- na.omit(data) # remove any data rows with na values
  data_dropped<- data %>% mutate(
    reason_dropped = NA
  )        
  rownames(data_dropped)<- gsub(" ", "_", rownames(data_dropped))
  # reading in tree file
  
  if ( is.phylo(tree.file) == T){ # if phylo object in environment
    tre<- tree.file  # is.phylo placement at beginning avoids logical T/F error with grepl on phylo object 
  } else {
  
  
  if (grepl("nwk", tree.file) && is.character(tree.file) == T){ # if path to newick tree
  tre<- read.newick(tree.file)
  }
  
  if (grepl("nex", tree.file)&& is.character(tree.file) == T){ # if path to nexus tree
    tre<- read.nexus(tree.file)
    }
    

  if ( is.null(tre)){
  stop("tree.file argument must be the name of a nexus or newick file, or a phylo object in the r environment")  
  }
  }
  
  tre$tip.label<- gsub( "_", " ", tre$tip.label)
  if(!is.null(drop.tips)){
    tre<- drop.tip(tre, drop.tips)
  }
  
  
# read in synonym files from directory or pull from environment
if (!is.null(synonym.table)) {
if (is.character(synonym.table)) {  
  
  if (grepl("txt", synonym.table)) {
    synonym.table.env<- read.table(synonym.table, quote = "")
  }
  
  if (grepl("txt", synonym.table)) {
    synonym.table.env<- read.csv(synonym.table, quote = "")
  }
}
  
  if (is.matrix(synonym.table) | is.data.frame(synonym.table)){ 
    synonym.table.env<- synonym.table
    
  }
  
    synonym.table.env[,1]<- gsub("_", " ", synonym.table.env[,1])
    synonym.table.env[,2]<- gsub("_", " ", synonym.table.env[,2])

    
  for (i in 1:nrow(synonym.table.env)) {
    tre$tip.label[tre$tip.label == synonym.table.env[1,i]] <- synonym.table.env[2,i]
  }
  }
  
  
  
  # changing format of phylo names to agree with rotl
  tre_rotlnames<- add_taxonomy(tre)
  
  #write.table(tre_rotlnames, file = paste0(filestring, "tre_rotlnames.txt"))  
  #tre_rotlnames<- read.table(paste0(filestring, "tre_rotlnames.txt"))
  print("done ROTL names phy")

  # keep original unmodified tree
  tre$node.label<- NULL
  tre.named.dup<- multi2di(tre)
  
  # replace tip names with rotl matches
  tre.named.dup<- tre_rotlnames
  dup.tips<- as.numeric(which(duplicated(tre.named.dup$tip.label)))
  if(length(dup.tips)>0){
  tre.named<- drop.tip(tre.named.dup, dup.tips) # remove any newly synonymized tips
  } else {tre.named<- tre.named.dup}
  
  # collapsing 0 length internal branches to polytomies
  tre.named$edge.length[tre.named$edge.length < shortBranch.thresh]<- 0
  tre.named.edit<- di2multi(tre.named)
  tre.named.edit$node.label<- NULL
  
  # repeat find short branches to drop short terminal branches
  # first, find overlap between data and tree
  # so that I'm not dropping tips unnecessarily
  # to begin, remove duplicates in data
  print("Tree renamed, converted to hard polytomies")
  # switch name formatting
  
  if ('binomial' %in% colnames(data)) {
    temp.colnames<- colnames(data)
    data_tax<- add_taxonomy(data)
    rownames(data)<-data_tax$binomial_updated
    data<- data[,temp.colnames]
    
  } else  {
    temp.colnames<- colnames(data)
    data<- rownames_to_column(data)
    colnames(data)<- c("binomial", temp.colnames)
    data_tax<- add_taxonomy(data)
    rownames(data)<-data_tax$binomial_updated
    data<- data[,temp.colnames]
  }
  
  rownames(data)<- gsub(" ", "_", rownames(data))
  tre.named.edit$tip.label<- gsub(" ", "_", tre.named.edit$tip.label)
  
  print(class(tre.named.edit))
  print(class(data))
  
  dat.rownames<- rownames(data)
  data<- as.matrix(data)
  rownames(data)<- dat.rownames
  
  treedat_prepDups<- treedata(tre.named.edit, data=data)
  temp_dropped<- rownames(data)[!rownames(data) %in% treedat_prepDups$phy$tip.label]
  data_dropped$reason_dropped[rownames(data_dropped) %in% temp_dropped]<- "No match between tip and data after rotl"
  all_duplicates<- vector()
  
  if(!is.null(cont_colnum)){  # Duplicated continues values removed for every continuous column. 
                              # Column numbers specified in function call. Does not run if no columns specified
  
    for (i in 1:length(cont_colnum)) { # repeated for each column
      temp_duplicates<- getDuplicates_cladewise(treedat_prepDups$phy, treedat_prepDups$data[,i])
      temp_duptax<- na.omit(unlist(temp_duplicates[,!(names(temp_duplicates) %in% c("trait", "clade_size"))]))
      all_duplicates<- append(all_duplicates, temp_duptax )
    }


  # finally, remove all tips with duplicated continuous values
  if(length(all_duplicates)>0){
  tre_noDups<- drop.tip(treedat_prepDups$phy, unique(all_duplicates))
  print(c("Dropped:", unique(all_duplicates)))
  } else {tre_noDups<- treedat_prepDups$phy}
  } else{
    
    tre_noDups<- treedat_prepDups$phy
  }
  data_dropped$reason_dropped[rownames(data_dropped) %in% unique(all_duplicates)]<- "Duplicated datapoint"
  
  print("Duplicated data points removed")
  dat.rownames<- rownames(data)
  data<- as.matrix(data)
  rownames(data)<- dat.rownames
  
  tredat<- treedata(tre_noDups, data)
  
  # finding erroneously long branches in the timetree.org phylogenies.... this happens often
  treeHeights<- data.frame(nodeHeights(tredat$phy))
  heightsVec<-sort(table(treeHeights$X2),decreasing=TRUE)[1:10]
  height<- as.numeric(names(heightsVec[1]))  # this finds most common node height in phylogeny, which will be tip height in ultrametric (or mostly extant) phylos

  tooLong<-which(treeHeights$X2> (height+0.01))
  
  if(length(tooLong)>= 1){
    tooLong_diff<-treeHeights$X2[tooLong]-height
    tredat$phy$edge.length[tooLong]<- tredat$phy$edge.length[tooLong]-tooLong_diff
    
    # check that trim was successful
    nodeHeights(tredat$phy)[tooLong,2]
    max(nodeHeights(tredat$phy))
  }
  
  tredat<- treedata(tredat$phy, data)
  print("Trimmed erroneously long branches")
  # and now make vector of tips with very short terminal branches
  # shortBranch.thresh determines minimum length of branch
  short_terminals<- tredat$phy$edge[which(tredat$phy$edge.length < shortBranch.thresh),]
  while( length(unlist(short_terminals)) > 0) {
    drop_tips<- vector()
  for (i in (1:length(unique(short_terminals[,1])))) {
    temp_node<- unique(short_terminals[,1])[i]
    temp_tips<- short_terminals[c(short_terminals[,1]==temp_node),2]

    if(length(temp_tips) >1){
    select_tip<- sample(temp_tips, (length(temp_tips)-1))  # drop all but one member of clades with rapid diversification
    } else {select_tip<- temp_tips}
    drop_tips<- append(drop_tips, select_tip)

  }
  if (length(drop_tips)>0) {
 data_dropped$reason_dropped[rownames(data_dropped) %in% tredat$phy$tip.label[drop_tips]]<- "Short terminal branch"
    
 tredat$phy<- drop.tip(tredat$phy, drop_tips)
  }
 short_terminals<- tredat$phy$edge[which(tredat$phy$edge.length < shortBranch.thresh),]
  }
  
  print(paste(length(short_terminals), "short terminals")) # if all is working right, this will be a zero length vector
  
  
  # confirm name formatting
  rownames(data)<- gsub(" ", "_", rownames(data))
  tredat$phy$tip.label<- gsub(" ", "_", tredat$phy$tip.label)
  
  print("Dropped tips with very short terminal branch lengths")
  
  species<- rownames(data)
  
  if(scale == T) {
  data <- data.frame(data) %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(., scale = FALSE))))
  }
  
  data<- as.matrix(data)
  rownames(data)<- species
  tredat.noshorts<- treedata(tredat$phy, data)
  
  write.nexus(tredat.noshorts$phy, file = paste0(filestring, "timetree.nex"))
  write.table(tredat.noshorts$data, file = paste0(filestring, "dat.txt"), col.names = F, quote = F, sep = "\t")
  write.table(tredat.noshorts$data[,1], file = paste0(filestring, "response_dat.txt"), col.names = F, quote = F, sep = "\t")

  
  if(write_predictors == TRUE){
    for (i in 2:ncol(tredat.noshorts$data)) {
      write.table(tredat.noshorts$data[,i], file = paste0(filestring, colnames(data)[i], "predictor_dat.txt"), col.names = F, quote = F, sep = "\t")
    }
  }
  return(list(tredat.noshorts, data_dropped))
  
}



