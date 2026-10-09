# plot betas and node scalars on trees
# node count vs fabric transformations plot
# trace plots
# data on tree
# posterior distributions of parameters

library(ggplot2)
library(ggtree)
library(phytools)
library(dplyr)
library(colorspace)
library(janitor)
library(bayestestR)
library(patchwork)

if (!interactive()) pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)

if (length(args)>0) {
  master_file <- args[1]
  varrates_file<- args[2]
  tre_file<- args[3]
  
  } else { # finding necessary tree plotting files
    
  all_files<- list.files()
  master_file<- all_files[grep("(?=.*master)(?=.*.csv)",all_files, perl = TRUE)][1]
  df_all<- read.csv(master_file, check.names = F)
  names(df_all) <- gsub("[()\\* ]", "", names(df_all)) 
  
  varrates_file_prep<- colnames(df_all)[grep("MeanBetaBLNZ-", colnames(df_all))]
  varrates_file<- unlist(strsplit(varrates_file_prep[1], "MeanBetaBLNZ-", fixed = TRUE))[2]
  
  log_file_prep<- sapply(varrates_file_prep, function (i) gsub(".VarRates.txt.csv", ".Log.txt", i))
  log_file_prep2<- sapply(log_file_prep, function(i) gsub ("MeanBetaBLNZ-", "", i))
  log_file<- log_file_prep2[1]
  
  tree_file_prep<- readLines(log_file, n = 50)
  tree_file_prep2<- tree_file_prep[grep("Tree File Name", tree_file_prep)]
  tre_file_prep3<- trimws(unlist(strsplit(tree_file_prep2, "Tree File Name:" ))[2])
  tre_file<- sub(".*/", "", tre_file_prep3)

}



  ################## Load tree & data
if (!exists("df_all")) {
  df_all<- read.csv(master_file) # used to plot significant node tranformations
}
  df<- read.csv(varrates_file) # used to match BT node numbers to R phylo object
  tree<- read.nexus(tre_file) # used for plotting and for aligning 
  
  all_logfiles<- list()
  for (i in 1:length(log_file_prep2)) {
    tbl<- readLines(log_file_prep2[i], n = 100)
    start_line <- grep("Iteration.*Lh|Lh.*Iteration", tbl)
    all_logfiles[[i]] <- read.delim(log_file_prep2[i], header = TRUE, skip = start_line-1)
    
  }
  names(all_logfiles)<- log_file_prep2
  print("all files read in")
  
  combined_logfiles<- bind_rows(all_logfiles, .id= "source")
  
  colnames(df_all)<- gsub("...", "_", colnames(df_all), fixed = T)
  colnames(df_all)<- gsub("..", "_", colnames(df_all), fixed = T)
  colnames(df_all)<- gsub(".", "_", colnames(df_all), fixed = T)
  colnames(df_all)<- sub("_$", "", colnames(df_all))
  
  
  nruns<- length(grep("MeanBetaBLNZ-", colnames(df_all)))
  meanbeta_cols<- grep("MeanBetaBLNZ-", colnames(df_all))
  meanscalar_cols<- grep("MeanNon1Scalar-", colnames(df_all))
  df_all<- clean_names(df_all)
  df<- clean_names(df)
  treeheight <- max(node.depth.edgelength(tree))
  
  outName <- gsub(".csv", "", master_file)
  homedir<- getwd()
  dir.create(paste0(outName, "_plots"))
  setwd(normalizePath(paste0(homedir, paste0("/", outName, "_plots", "/"))))
  
  
  ################## Align fabric post-processor output with Ape tree format
  # Sort post-processor by branch length, then by taxon name
  df = df %>% arrange(`branch_length`, `taxa_1`)
  df = df[-1, ] # Delete 1st row
  
  # Process tree to df
  tree_df = tidytree::as_tibble(ape::makeNodeLabel(tree, method = "number"))
  colnames(tree_df) = c( "parent_node", "desc_node", "branch_len", "label")
  tree_df = tree_df %>% arrange(`branch_len`, `desc_node`)
  tree_df = tree_df[-nrow(tree_df), ] # Delete last row
  
  # Align Data
  df_aligned = data.frame(
    df1_selected = df[, c("id", "branch_length", "mean_beta_bl_nz", "mean_non_1_scalar","mean_beta_bl" )],
    df2_selected = tree_df[, c("desc_node", "branch_len", "label")]
  )
  
  
  # select significant nodes
  nodes_BTid<- df_all[df_all$no_sig_nodes == nruns & df_all$taxa_3!="", 
                      ]$id
  test<- df_all[df_all$no_sig_nodes == nruns,]
  nodes_BTid_tbl<- df_all[df_all$no_sig_nodes == nruns,] #& df_all$taxa_2!="", ]$ID
  nodes_BTid_tbl<- nodes_BTid_tbl[nodes_BTid_tbl$taxa_3 != "",]
  
  nodes_num<- df_aligned$df2_selected.desc_node[df_aligned$df1_selected.id %in% nodes_BTid]
  
  
  # select significant betas
  betas_BTid<- df_all[df_all$no_sig_betas == nruns & df_all$taxa_2!="", 
  ]$id
  betas_BTid_tbl<- df_all[df_all$no_sig_betas == nruns,] #& df_all$taxa_2!="", ]$ID
  betas_BTid_tbl<- betas_BTid_tbl[betas_BTid_tbl$taxa_2 != "",]
  betas_num<- df_aligned$df2_selected.desc_node[df_aligned$df1_selected.id %in% betas_BTid]
  betas_vals<- rowMeans(betas_BTid_tbl[,c(meanbeta_cols[1]: meanbeta_cols[length(meanbeta_cols)])])
  
  
  # Cleanup data for plotting betas
  data_to_map = data.frame(df_aligned$df2_selected.desc_node, df_aligned$df1_selected.mean_beta_bl_nz, df_aligned$df1_selected.mean_non_1_scalar, df_aligned$df1_selected.mean_beta_bl)
  
  colnames(data_to_map) = c( "node", "betas", "nscalars", "betas_Nonsignif")
  data_to_map$betas = round(data_to_map$betas, 1) # Round betas
  data_to_map$nscalars = round(data_to_map$nscalars, 1) # Round node scalars
  
  # get mean beta values
  beta_vals<- df_all[df_all$id %in% betas_BTid,][,c(meanbeta_cols[1]: meanbeta_cols[length(meanbeta_cols)])]
  beta_vals[,(nruns+1)]<- rowMeans(beta_vals)
  rownames(beta_vals)<- betas_BTid
  
  # get mean node scalar values
  node_vals<- df_all[df_all$id %in% nodes_BTid,][,c(meanscalar_cols[1]: meanscalar_cols[length(meanscalar_cols)])] 
  node_vals[,(nruns+1)]<- rowMeans(node_vals)
  rownames(node_vals)<- nodes_BTid
  
  ########## Plot tree with betas and scalars

  # setting plot dimensions for outfile
  dim_h<- (length(tree$tip.label)*0.07)

  # Create the plot and view
  # https://pmassicotte.github.io/paletteer_gallery/
  p = ggtree(tree, aes(color=nscalars),  layout = "fan", size=0.5) %<+% data_to_map +
    #paletteer::scale_colour_paletteer_c("grDevices::Spectral")
    #theme(text=element_text(family = "Nimbus Sans"))+
    scale_colour_continuous_diverging(palette ="Blue-Red", mid=1.2, l1=0.1, c1=100, p1=0.2, p2=0.5, rev=F
    )
  
  #scale_colour_gradient2(low="#009593FF", mid= "#EDE5C6FF",high = "#D35C79FF", midpoint = 4)
  # scale_colour_paletteer_c("grDevices::Zissou 1")
  # scale_color_paletteer_c("grDevices::Spectral", direction=-1)
  # scale_colour_paletteer_c("oompaBase::greenscale", direction=1) # Matrix trees
  # scale_color_paletteer_c("oompaBase::redscale") # Red matrix
  # scale_color_paletteer_c("pals::coolwarm")
  # scale_color_paletteer_c("harrypotter::ravenclaw") 
  # scale_color_paletteer_c("pals::kovesi.rainbow_bgyr_35_85_c72")
  p +
    #geom_tiplab(offset=0.01, color="black", fontface=3, size=3) +
    geom_point2(aes(subset= node %in% nodes_num ), fill= "firebrick", color= "black", shape= 21, size=3, stroke=0.2)+
    #geom_text2(aes(subset = node %in% top10_betas$node, label = betas), 
    #           hjust = -0.1, color = "black", size = 4) +
    ggplot2::xlim(0, (treeheight*1.5) )+
    geom_tiplab(size=3, color= "gray40")
  
  
  ggsave(filename = paste0(outName, "_scalarTrePlot.svg"),  width=dim_h, height=dim_h, units="in", 
        dpi=300, limitsize = F)
  
  # betas plotting
  p_betas = ggtree(tree, aes(color=betas), layout = "fan", size=0.5) %<+% data_to_map +
    #paletteer::scale_colour_paletteer_c("grDevices::Spectral")
    #theme(text=element_text(family = "Nimbus Sans")) +
    paletteer::scale_colour_paletteer_c("grDevices::PiYG", direction=-1) 
  
  #scale_colour_continuous_diverging(palette ="Blue-Red", mid=1.2, l1=0.1, c1=100, p1=0.2, p2=0.5, rev=F)
  p_betas +
    #geom_tiplab(offset=0.01, color="black", fontface=3, size=3) +
    geom_point2(aes(subset= node %in% betas_num), fill = "firebrick", color= "black", shape= 21, size=3, stroke=0.2)+
    #geom_text2(aes(subset = node %in% top10_betas$node, label = betas), 
    #           hjust = -0.1, color = "black", size = 4) +
    ggplot2::xlim(0, (treeheight*1.5))+
    geom_tiplab(size=3, color= "gray50")
  
  ggsave(filename = paste0(outName, "_betaSignificantTrePlot.svg"), width=dim_h, height=dim_h, units="in", 
          limitsize = F)
  
  
  
  
  # betas plotting
  p_betas_v2 = ggtree(tree, aes(color=betas_Nonsignif),  layout = "fan", size=0.5) %<+% data_to_map +
    #paletteer::scale_colour_paletteer_c("grDevices::Spectral")
    #theme(text=element_text(family = "Nimbus Sans")) +
  
    paletteer::scale_colour_paletteer_c("ggthemes::Temperature Diverging",
                                        direction=-1, 
                                        limits = c(unlist(ci(data_to_map$betas_Nonsignif))[c(2,3)]),
                                        oob = scales::squish
                                        ) 
  
  #scale_colour_continuous_diverging(palette ="Blue-Red", mid=1.2, l1=0.1, c1=100, p1=0.2, p2=0.5, rev=F)
  p_betas_v2 +
    #geom_tiplab(offset=0.01, color="black", fontface=3, size=3) +
    geom_point2(aes(subset= node %in% betas_num), fill = "firebrick", color= "black", shape= 21, size=3, stroke=0.2)+
    #geom_text2(aes(subset = node %in% top10_betas$node, label = betas), 
    #           hjust = -0.1, color = "black", size = 4) +
    ggplot2::xlim(0,  (treeheight*1.5))+
    geom_tiplab(size=3, 
                #size= 0.3*dim_h
                color= "gray40")
  
  
  ggsave(filename =  paste0(outName, "_betaNonSignifTrePlot.svg"), width=dim_h, height=dim_h, units="in", 
limitsize = F)
  

# trace plots for logfiles

  alpha_trace= 0.5
  # likelihood
lh_plot<-ggplot(combined_logfiles, aes(x= Iteration, y=Lh, group = source))+
  geom_line(size = 0.3, aes(color = source, fill = source), alpha = alpha_trace)+
  theme_light()+
  #paletteer::scale_color_paletteer_d("ggprism::warm_pastels")+

theme(legend.position = "none", plot.title = element_text(size = 10))+
  labs(title="Trace Plot of Posterior Likelihood, all Fabric Chains")  
  lh_plot
  
  
var_plot<- ggplot(combined_logfiles, aes(x= Iteration, y=Var, group = source))+
    geom_line(size = 0.2, aes(color = source, fill = source), alpha = alpha_trace)+
    theme_light()+
  #paletteer::scale_color_paletteer_d("ggprism::warm_pastels")+
  
    theme(legend.position = "none", plot.title = element_text(size = 10))+
  labs(title="Trace Plot of Variance, all Fabric Chains")  
  
alpha_plot<- ggplot(combined_logfiles, aes(x= Iteration, y=Alpha, group = source))+
    geom_line(size = 0.2, aes(color = source, fill = source), alpha = alpha_trace)+
    theme_light()+
  #paletteer::scale_color_paletteer_d("ggprism::warm_pastels")+
  
    theme(legend.position = "none", plot.title = element_text(size = 10))+
    labs(title="Trace Plot of Alpha, all Fabric Chains")  
  
  beta_cols<- grep("^Beta|^beta", colnames(combined_logfiles))
  beta_names<- colnames(combined_logfiles)[beta_cols]
  beta_plots<- list()
  
  for (i in beta_cols) {
    
  iter<- which(colnames(combined_logfiles)[i] %in% beta_names)
  temp_name<- beta_names[iter]
  
  beta_plots[[i]]<- ggplot(combined_logfiles, aes(x= Iteration, y= !!sym(temp_name), group = source))+
    geom_line(size = 0.2, aes(color = source, fill = source), alpha = alpha_trace)+
    theme_light()+
    #paletteer::scale_color_paletteer_d("ggprism::warm_pastels")+
    
    theme(legend.position = "none", plot.title = element_text(size = 10))+
    labs(title="Trace Plot of Beta, all Fabric Chains")  
    
  }
  
  lh_plot+var_plot+alpha_plot+unlist(beta_plots)
  ggsave( filename = paste0(outName, "_",  "AllFabricTrace.png"), width = 11, height = 8)
  

    


