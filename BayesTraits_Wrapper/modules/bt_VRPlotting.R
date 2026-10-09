# Chris Organ 2026
# Calculates variable rates from BayesTraits analysis and maps onto a time tree
# Can handle tree sampling during the VR
# Plot the rate gradient on the time tree (credit to Emi Bender)
#
# Input: VR tree posterior and original (time calibrated) tree or trees
# Output: Plot of time tree with rates mapped on

library(phytools)
library(phangorn) # provides maxCladeCred()
library(ggtree) # BiocManager::install("ggtree")
library(ggplot2)
library(paletteer)
library(tidytree)
library(dplyr)
library(ggnewscale)
library(deeptime)
library(geiger)


# PROCESS BRANCH LENGTHS ----------------------------------------------------
# Can handle for a distribution of VR topologies or a single topology
# Uses prop.part (standard ape function) to get Clade IDs instantly
get_clade_data = function(phy) {

  # Tips: ID is the label, Node is 1..Ntip
  tips = tibble(
    node = 1:Ntip(phy),
    id   = phy$tip.label,
    len  = phy$edge.length[match(1:Ntip(phy), phy$edge[,2])]
  )
  
  # Internal: ID is sorted descendant tips, Node is (Ntip+1)..End
  # prop.part returns a list where index i corresponds to node (Ntip + i)
  parts = prop.part(phy)
  internal_ids = sapply(parts, function(x) paste(sort(phy$tip.label[x]), collapse="|"))
  
  internal = tibble(
    node = (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy)),
    id   = internal_ids,
    len  = phy$edge.length[match((Ntip(phy) + 1):(Ntip(phy) + Nnode(phy)), phy$edge[,2])]
  )
  
  return(bind_rows(tips, internal))
}


# LOAD & PROCESS DATA ----------------------------------------------------
# NOTE: Trees saved from BayesTrees may cause R to crash...
# timetree  : the time-calibrated tree (phylo/multiPhylo object) the VR analysis was run on
#                  if this is a tree set, the mcc is used
# trees_vr  : the posterior of VarRates trees (multiPhylo object) output by BayesTraits
#             (read in with e.g. ape::read.nexus() before calling this function)
# todrop    : optional vector of tip labels to drop from timetree (e.g. taxa cut from
#             the VR analysis due to missing data) before mapping rates onto it
# Returns a list: mapped_data (rates per clade) and mcc_tree (the MCC tree the
# rates were mapped onto — pass this into VR_plot(), not the original timetree).
mapped_rates = function(timetree, trees_vr, todrop) {

  # Drop tips so time tree matches VR tree
  if (!is.null(todrop)) {
    timetree = ape::drop.tip(timetree, todrop)
  }

  # Aggregate Posterior (VR Trees)
  # Calculate mean length for every unique clade across all trees
  vr_means = bind_rows(lapply(trees_vr, get_clade_data)) %>%
    group_by(id) %>%
    summarize(mean_VRBL = mean(len, na.rm = TRUE))

  # MCC Tree
  timetree_mcc = maxCladeCred(timetree, tree=TRUE, rooted=TRUE)

  # Map to MCC (Time Tree)
  mapped_data = get_clade_data(timetree_mcc) %>%
    rename(TBL = len) %>%
    left_join(vr_means, by = "id") %>%
    mutate(rates = round(mean_VRBL / TBL, 2)) %>%
    arrange(desc(rates)) # Sort for top 10

  return(list(mapped_data = mapped_data, mcc_tree = timetree_mcc))
}


# PLOT VR HEATMAP----------------------------------------------------
# Palette gallery: https://pmassicotte.github.io/paletteer_gallery/
# timetree    : the MCC tree the rates were mapped onto — use mapped_rates()'s
#               $mcc_tree return value, not the original (possibly multi-tree) timetree
# mapped_data : the $mapped_data return value from mapped_rates()
# outdir      : directory to save the plot into (created if missing)
# outName     : base filename (job name) used for the saved plot
VR_plot = function(timetree, mapped_data, color_palette){

  p_rates = ggtree(timetree, aes(color = rates), size = 0.35) %<+% mapped_data +
    scale_color_paletteer_c(color_palette, na.value = "grey80", direction = 1) +
    geom_tiplab(offset = 0.002, color = "black", fontface = 3, size = 3) +
    geom_text2(aes(subset = node %in% head(mapped_data$node, 20), label = rates),
              hjust = -0.2, color = "red", size = 5) +
    xlim(0, max(node.depth.edgelength(timetree)) * 1.2)

  # Save
  fig_save_path = file.path(path_results, paste0(job_name, "_results/", job_name, "_rates.svg"))
  ggsave(fig_save_path, p_rates, width=8, height=6)

}



# p_geo = ggtree(scaled_tree, aes(color=rates), size=0.75) %<+% target_data +
#   scale_color_paletteer_c("pals::kovesi.rainbow_bgyr_35_85_c72") +
#   coord_geo(xlim=c(-80, 5), ylim=c(-2, Ntip(scaled_tree)+2), neg=TRUE, abbrv=FALSE, dat=list("epochs", "periods"), pos=list("bottom", "bottom")) +
#   scale_x_continuous(breaks = seq(-80, 0, 20), labels=abs(seq(-80, 0, 20))) +
#   scale_y_continuous(guide = NULL) +
#   theme_classic() +
#   guides(r = guide_axis_stack(guide_geo("epochs", neg=TRUE, rot=-90, size="auto", height=unit(1, "line")), guide_axis(), spacing=unit(0, "line"))) +
#   #geom_tiplab(offset=0.01, color="black", fontface=3, size=2.5) +
#   geom_text(aes(label = node_ages), hjust = -0.3, size = 3) # Adjust `hjust` for label placement
# p_geo = revts(p_geo)
# print(p_geo)

# invisible(list(rates_data = target_data, p_rates = p_rates, p_geo = p_geo))

  