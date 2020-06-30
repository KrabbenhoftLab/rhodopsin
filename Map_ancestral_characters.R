#Script for visualizing the output of an ancestral state reconstruction by mapping ancestral characters onto a phylogenetic tree
#KM Eaton, June 2020, University at Buffalo

#ABOUT THIS SCRIPT:
#This script takes a tree with branch lengths in Newick format, a tree with named ancestors in newick format (output from RAxML)
#And two csv files - one containing info on the ancestral states and one containing info on the states of the tips of the tree. 
#These csvs should have two columns each - one with the species/ancestor name, and one with the state (this column was called trait for me).

#USER INPUT:
#Since this is a different format than what we usually use with these scripts, I'll be defining user inputs as we go along.
#Any time you see a file name in quotes (it'll show up green if you're using Rstudio, which you should be), that'll be something the user has to define. 
#I'll keep it annotated as clearly as I can to prevent any confusion!

#STEP 1: Set your working directory. Make sure all the files you'll need are in this directory, including your two tree files and your two csvs.
setwd("~/Dropbox/Documents/Ciscoes/FTOL-analysis")

#STEP 2: Load in the following libraries. If you don't have them installed yet, you'll have to do: install.packages("libraryname")
require(tibble)
require(ggtree)
require(treeio)
require(ape)
require(dplyr)
require(tidytree)
require(phytools)

#STEP 3: Merge your two tree files. This step will ONLY work if your two tree files have exactly the same topology.
#If they don't have the same topology, that's weird and you did something wrong with RAxML.
#The first tree we'll load in is the one you input into RAxML. This has branch lengths, but the ancestors aren't labelled. (Replace "trimmed-time-tree.nw" with your own file.)
ultrametric <- read.tree("trimmed-time-tree.nw")
#The second tree we'll load in is the one you got in the output of RAxML. This has labelled ancestors, but no branch lengths. (Replace "RAxML_node_labelled_rooted_tree.tre" with your own file.)
labeled <- read.tree("RAxML_node_labelled_rooted_tree_6_23_20.tre")
#This next step makes a new list in the ultrametric tree; the node labels in the labeled tree fill that list.
ultrametric$node.label <- labeled$node.label
#This last line is optional, but it's good to do. It saves the "ultrametric" tree to a file, where now the nodes of the tree have labels according to how they were assigned in RAxML. 
write.tree(file = "branch_lengths_and_nodes_tree.tre", ultrametric)
#Now the "ultrametric" object is ready for us to use in our analyses.

#STEP 4: Merge your two csv files.
#I personally think it's much easier to do this step completely in a text editor like sublime text, so that's what I'm going to describe here, but if you know how to do it in R, be my guest.
#You should have two csvs: one containing the ancestral node names and associated trait, and one containing the tip names and associated trait.
#Copy the entire content of the tip file into the ancestral node file. Just add it right to the end. This should basically double the length of your file (in terms of line numbers).
#Add a line at the start of the file to give the columns names. You want the column with the ancestral node names & species names to be called "label" and the one with the trait to be called "trait". (Exclude the quotes from the actual names, just trying to make it clear here.) Make sure to include a comma as a separator.
#The order of the columns should be label,trait.
#Finally, do a search and replace in the trait column. You want to assign your traits numerical values. I only had two traits, so I just made mine 0 and 1, but you can go higher if you need to. Just make sure you write it down somewhere what the numbers correspond to.
#Search for ,traitvalue0 (I searched for ,F when I was looking to replace phenylalanine), and replace it with ,0.
#Do the same with the rest of your traits.
#For Katie: 0 = Phe, 1 = Tyr.

#STEP 5: Load this new, big csv into R, and incorporate the trait data into your tree file. 
#This loads your csv file into R. Replace "RAxML_ancestral_chars_numeric.csv" with the name of the big csv file you generated in STEP 4.
node_and_tip_traits <- read.csv("RAxML_ancestral_chars_numeric.csv")

#This changes the format of your tree file we generated way back in STEP 3. 
tibble_tree <- as_tibble(ultrametric)
#The following line is optional, just run it if you want to see what the data looks like. It's not that exciting. 
tibble_tree

#This adds your trait data from the csv file to the tree file. It merges them by the column "label", which contains the species/node names.
joined_tibble <- full_join(tibble_tree, node_and_tip_traits, by = 'label')

#This just converts between file formats for the tree again.
tree_data_object <- as.treedata(joined_tibble)
#This line is optional, just lets you see what you just created.
tree_data_object

#STEP 6: Plot a super cool tree with your traits on it! 
#First, you'll need to load the tree into ggtree. Use aes to tell ggtree that you want it to color the branches and tips of the tree by the "trait" column. 
#You can change the layout of the tree if you want, but I think it looks best circular.
#You can also change the size of the tree. 0.15 is pretty small, but I had a giant tree with lots of taxa (>2000), so if you have less, you can definitely increase the size.
#geom_rootedge just lets you see the root of the tree 
#geom_tiplab lets you see the tip labels (i.e., the species names), and change the size of them. Again, mine are pretty small since the tree has so many taxa.
#you have to have geom_tiplab() invoked to be able to see the species names at all
#geom_treescale lets you include a scale for branch lengths
#scale_color_manual lets you set colors for each trait. I set it equal to "blue" (for 0) and "red" (for 1), but you can do whatever you'd like.
#And I hate the legend that R puts in so I usually get rid of it with theme(legend.position = 'none')
#Finally, when you use aes in the initial ggtree, it colors the tip labels the same as the tips of the tree. If you wanted all the tip labels to be the same color (like the default black), you can remove aes from the initial calling of ggtree and instead do + geom_tree(aes(color = trait))
ggtree(tree_data_object, aes(color = as.factor(trait)), layout = "circular", size = 0.15) + 
  geom_rootedge() + 
  geom_tiplab(size = 0.2) +
  geom_treescale() + 
  scale_color_manual(values = c("blue","red")) + 
  theme(legend.position = "none")

#OUTPUT: A super cool colored tree file! You can save it as a pdf from the plot window, or export it as a svg file using the following:
ggsave("FTOL_RAxML_traits.svg", dpi = 300)
