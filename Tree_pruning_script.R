##This script was modified from a script produced by WL Testo and JA Pelosi at the University of Florida

##set working directory

setwd("~/Dropbox/Documents/Ciscoes/FTOL-analysis")

##load geiger and dependencies

library("geiger")
library("dplyr")
library("phytools")
library("ape")

##load tree
tree<-read.tree("full_time_tree_Rabosky.nw")
tree<-ladderize(tree,right=F)

##check tree
plot(tree, cex=0.3)

##read in taxon list
data<-read.table("list_of_taxa_to_keep.txt", header = F)

##check data
data

###check names in list against tree
name_list<-name.check(tree,data)

##make list of names not in data file
checked_names<-name_list$tree_not_data

##convert to a vector
names_to_drop<-as.vector(checked_names)

##trim tree to just taxa represented in list
trimmed_tree<-drop.tip(tree,names_to_drop)

##visualize trimmed tree
plot(trimmed_tree, cex=0.35)

write.tree(trimmed_tree)
