

library(igraph)
library(Matrix)
library(IRanges)
library(ggplot2)

# read OBO format and return an igraph with all "is_a" relations
read.obo <- function(obo.file) {
  pat <- "^([^:]+) *: *(.*)$"
  obo <- readLines(obo.file)
  obo <- obo[grepl("^\\[",obo) | grepl(pat,obo)]
  obo <- DataFrame(
    section = cumsum(grepl("^\\[",obo)),
    key = sub(pat,"\\1",obo),
    value = sub(pat,"\\2",obo)
  )
  obo <- subset(obo,section %in% section[key=="[Term]"])
  obo <- subset(obo,key %in% c("is_a","name","id","is_obsolete","namespace"))
  obo$section <- factor(obo$section)
  obo <- DataFrame(
    goid = with(subset(obo,key=="id"),drop(CharacterList(split(value,section)))),
    term = with(subset(obo,key=="name"),drop(CharacterList(split(value,section)))),
    is_obsolete = with(subset(obo,key=="is_obsolete"),drop(CharacterList(split(value,section)))) %in% "true",
    namespace = with(subset(obo,key=="namespace"),drop(CharacterList(split(value,section)))),
    is_a = with(subset(obo,key=="is_a"),CharacterList(split(value,section)))
  )
  obo$is_a <- sub(" *!.*$","",obo$is_a)
  
  is_a <- stack(setNames(obo$is_a,obo$goid))
  is_a$type <- "is_a"
  graph_from_data_frame(is_a,vertices = as.data.frame(obo[1:4]))
}


# Test GOterm enrichement for the given list of genes
test.gene.set.enrichment <- function(genes,gene.sets,gene.universe=NULL) {
  # remove genes that that are not annotated
  gene.sets <- gene.sets[lengths(gene.sets$sets)>0,]
  
  # keep only annotation of the genes in the universe
  if (is.null(gene.universe)) {
    gene.universe <- unlist(gene.sets$gene.names)
  }
  gene.sets <- gene.sets[any(gene.sets$gene.names %in% gene.universe),]

  # determine number of gene with an annotation
  num_genes <- sum(any(gene.sets$gene.names %in% genes))
  
  # compute number of gene into each set
  N <- table(unlist(gene.sets$sets))
  n <- table(unlist(gene.sets$sets[any(gene.sets$gene.names %in% genes)]))
  
  # extract the list of gene for each set
  U <- stack(gene.sets$sets)
  U$gene <- unstrsplit(gene.sets$gene.names[gene.sets$gene.names %in% genes],"|")[as.integer(U$name)]
  U$value <- factor(U$value)
  U <- U[U$gene!="",]
  U <- splitAsList(U$gene,U$value)
  
  # Create the output object
  gsea <- DataFrame(GeneSet=names(N),N=as.vector(N),n=as.vector(n[names(N)]),genes=U[names(N)])  
  gsea$n[is.na(gsea$n)] <- 0
  gsea$pval <- 1 - phyper(gsea$n-1L,gsea$N,nrow(gene.sets)-gsea$N,num_genes)
  gsea$qval <- p.adjust(gsea$pval,"fdr")
  gsea <- gsea[order(gsea$pval,gsea$N),]
  gsea
}



hgnc.gene.sets <- function() {
  g <- local({
    fam <- read.table("data/family.csv.gz",sep=",",header=TRUE,stringsAsFactors = FALSE)
    fam$fam_name <- fam$name
    fam$name <- NULL
    h <- read.table("data/hierarchy_closure.csv.gz",sep=",",header=TRUE,stringsAsFactors = FALSE)
    h <- subset(h,child_fam_id!=parent_fam_id)
    graph_from_data_frame(h[c(2,1,3)],vertices=fam)  
  })
  
  
  A <- local({
    A <- ancestors(g)
    ij <- Matrix::which(A,arr.ind=TRUE)
    D <- DataFrame(node = rownames(A)[ij[,1]], ancestor = colnames(A)[ij[,2]])
    splitAsList(D$ancestor,D$node)
  })
  
  x <- read.table("data/hgnc_complete_set_Enriched.txt",sep="\t",comment="",quote='"',stringsAsFactors = FALSE,header=TRUE)
  x <- DataFrame(x)
  x <- x[x$status=="Approved",]
  x$gene_family_id <- CharacterList(strsplit(x$gene_family_id,"|",fixed=TRUE))
  x$mgd_id <- CharacterList(strsplit(x$mgd_id,"|",fixed=TRUE))
  x$mgd_symbol <- local({
    z <- read.table("data/mgi2ensembl.tsv.gz",sep="\t",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)  
    extractList(z$`Gene name`,match(x$mgd_id,z$`MGI ID`))
  })
  x <- x[lengths(x$mgd_id)>0 & lengths(x$mgd_symbol),]
  
  x$gene_family_ancestors_id <- local({
    U <- stack(x$gene_family_id)
    U$ancestors <- A[U$value]
    U$value <- NULL
    U <- splitAsList(unname(unlist(U$ancestors)),U$name[togroup(U$ancestors@partitioning)])
    unique(U)
  })
  
  x$gene_family_ancestors_name <- relist(V(g)[unlist(x$gene_family_ancestors_id)]$fam_name,x$gene_family_ancestors_id )
  
  
  GL <- merge(stack(x$mgd_symbol),stack(x$gene_family_ancestors_name),by="name")
  GL <- splitAsList(GL$value.y,GL$value.x)
  GL <- DataFrame(gene.names = as(names(GL),"CharacterList"),sets = GL)
  GL
}















