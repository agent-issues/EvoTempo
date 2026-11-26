# Starting tree for the RevBayes chains, from the starting tree of the BEAST2
# analysis (8+dna_R.xml, TreeParser newick).
#
# Chains started from random FBD trees stalled with a crown age of 3-7 Ma and
# the crown fossils (Spheniscus urbinai, S. megaramphus, S. muizoni,
# Pygoscelis grandis, Madrynornis mirandus) outside the crown (in 8 chains,
# after ~4 000 iterations, the only exception was P. grandis inside a 4.7 Ma
# crown, at a sampled age of 3.1 Ma, in one chain): a fossil cannot
# enter a crown younger than itself, and the crown cannot grow older past the
# fossils sitting on its stem, so the chains did not leave that region (see
# README).  The XML's starting tree has these fossils inside a crown of
# 11.3 Ma.  Its tip heights do not all fall within the fossil age intervals
# (BEAST2 resets them), so:
#   * each fossil age outside its interval is moved inside it, 1 % of the
#     interval width from the nearer bound;
#   * each internal node is made at least 0.01 Myr older than its oldest child
#     (sampled ancestors in the XML tree become short side branches; the chain's
#     collapse/expand move restores them);
#   * extant tips are set to age 0 exactly.
#
# Run from this directory after prepare-data.R: Rscript start-tree.R  (needs ape)

library(ape)
xmlFile <- unzip("../raw/xmls.zip", files = "xmls/8+dna_R.xml", exdir = tempdir())
xml <- paste(readLines(xmlFile, warn = FALSE), collapse = "\n")
newick <- sub(".*<input name='newick'>\\s*([^<]*?)\\s*</input>.*", "\\1", xml)
tr <- read.tree(text = paste0(sub(";?$", "", newick), ";"))

tips <- read.delim("data/taxa.tsv")
stopifnot(setequal(tr$tip.label, tips$taxon))
nTip <- Ntip(tr)
depth <- node.depth.edgelength(tr)
extant <- tips$taxon[tips$max_age == 0]
age <- max(depth[match(extant, tr$tip.label)]) - depth

lo <- tips$min_age[match(tr$tip.label, tips$taxon)]
hi <- tips$max_age[match(tr$tip.label, tips$taxon)]
margin <- 0.01 * (hi - lo)
tipAge0 <- age[seq_len(nTip)]
tipAge <- pmin(pmax(tipAge0, lo + margin), hi - margin)
tipAge[hi == 0] <- 0
age[seq_len(nTip)] <- tipAge

# edges in post-order: a child's age is final before its parent is fixed
for (i in postorder(tr)) {
  parent <- tr$edge[i, 1]
  child <- tr$edge[i, 2]
  age[parent] <- max(age[parent], age[child] + 0.01)
}
tr$edge.length <- age[tr$edge[, 1]] - age[tr$edge[, 2]]
stopifnot(all(tr$edge.length > 0))

crown <- getMRCA(tr, extant)
cat(sprintf("crown age %.3f Ma, root age %.3f Ma; fossils moved into interval: %d\n",
            age[crown], age[nTip + 1],
            sum(hi > 0 & abs(tipAge - tipAge0) > 1e-6)))
write.tree(tr, "data/start_tree.nwk", digits = 10)
