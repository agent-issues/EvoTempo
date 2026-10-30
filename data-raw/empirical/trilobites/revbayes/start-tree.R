# Starting tree for the trilobite RevBayes chains: the median-height
# consensus tree of Paterson et al. 2019's BEAST2 FBD analysis with sampled
# ancestors (BEAST2_contrees/CambrianTrilos_106_DisGamStr_SABD_time_Median_con.tre
# in ../raw/Dataset_4_BEAST_and_BEAST2.zip).
#   * heights there are above the youngest tip (486 Ma); they are converted to
#     Myr before the end of the sampling window (485.4 Ma; prepare-data.R);
#   * tips are set to their (fixed) ages;
#   * zero-length sampled-ancestor branches become 0.01-Myr side branches and
#     each internal node is made at least 0.01 Myr older than its oldest child
#     (the chain's collapse/expand move restores sampled ancestors);
#   * the ingroup constraint is checked.
#
# Run from this directory after prepare-data.R: Rscript start-tree.R  (needs ape)

library(ape)
member <- "SI_file4_BEAST_and_BEAST2/BEAST2_contrees/CambrianTrilos_106_DisGamStr_SABD_time_Median_con.tre"
f <- unzip("../raw/Dataset_4_BEAST_and_BEAST2.zip", files = member, exdir = tempdir())
tr <- read.nexus(f)
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
tips <- read.delim("data/taxa.tsv")
stopifnot(setequal(tr$tip.label, tips$taxon))
tr <- multi2di(tr)
nTip <- Ntip(tr)
depth <- node.depth.edgelength(tr)
youngest <- which.min(tips$min_age[match(tr$tip.label, tips$taxon)])
age <- depth[youngest] - depth + tips$min_age[match(tr$tip.label[youngest], tips$taxon)]
age[seq_len(nTip)] <- tips$min_age[match(tr$tip.label, tips$taxon)]
for (i in postorder(tr)) {
  parent <- tr$edge[i, 1]
  child <- tr$edge[i, 2]
  age[parent] <- max(age[parent], age[child] + 0.01)
}
tr$edge.length <- age[tr$edge[, 1]] - age[tr$edge[, 2]]
stopifnot(all(tr$edge.length > 0), is.binary(tr), is.rooted(tr))
ingroup <- readLines("data/clade_ingroup.txt")
stopifnot(is.monophyletic(tr, ingroup))
write.tree(tr, "data/start_tree.nwk", digits = 10)
cat(sprintf("start tree: root %.2f Ma, ingroup %.2f Ma (absolute)\n",
            age[nTip + 1] + 485.4, age[getMRCA(tr, ingroup)] + 485.4))
