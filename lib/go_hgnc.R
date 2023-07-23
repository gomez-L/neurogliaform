
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Script performing Gene Set Enrichment Analyses from GO gene sets and HGNC gene sets
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

source("src/julien/lib/go_hgnc_lib.R")

go <- function(universe,genelist) {
  
  go <- read.obo(gzfile("data/go-basic.obo.gz"))
  
  gene.sets <- readRDS("data/data/mouse_go_gene_sets.rds")
  
  # remove global GO terms that are associated to almost all genes
  gene.sets$sets <- gene.sets$sets[!gene.sets$sets %in% c("molecular_function [GO:0003674]","cellular_component [GO:0005575]","biological_process [GO:0008150]")]
  gene.sets <- gene.sets[lengths(gene.sets$sets)>0,]

  universe <- readLines("out/universe.txt")
  genes <- intersect(genelist,universe)
  
  # Compute an enrichment
  X <- test.gene.set.enrichment(genes,gene.sets,gene.universe = universe)
  X$go_id <- sub(".*\\[(GO:[0-9]+)\\]","\\1",X$GeneSet)
  X$go_nm <- V(go)[X$go_id]$namespace

}


hgnc <- function(universe,genelist) {
  
  gene.sets <- hgnc.gene.sets()
  
  universe <- readLines("out/universe.txt")
  genes <- intersect(genelist,universe)
  
  X <- test.gene.set.enrichment(genes,gene.sets,gene.universe = universe)

}


