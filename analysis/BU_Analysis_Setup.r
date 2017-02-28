##############################################################################
### SETUP
##############################################################################

# set the repoDir to be the 'ipfcore' repo directory
repoDir <- "/protected/projects/pulmseq/Genentech_IPF_Analysis/ipfcore/"
# move to the repoDir so everything can be loaded and saved appropriately
setwd(repoDir)

##############################################################################
### LOAD LIBRARIES
##############################################################################
require(DESeq2)
require(corrplot)
require(org.Hs.eg.db)
require(ggplot2)
require(enrichR)

# For use with enrichR
dbs <- c("Human_Gene_Atlas", "Biocarta_2016",
         "KEGG_2016", "WikiPathways_2016",
         "Reactome_2016", "GO_Biological_Process_2015",
         "GO_Cellular_Component_2015", "GO_Molecular_Function_2015")

##############################################################################
### LOAD DATA
##############################################################################
vsd <- readRDS("data/Genentech_IPF_vsd.RDS")
dds <- readRDS("data/Genentech_IPF_dds2.RDS")

##############################################################################
### WRITE FUNCTIONS
##############################################################################
# For use in combining the sequencing and phenotype data
remove0s <- function(str){
  temp <- str
  str_length <- nchar(temp)
  last_char <- substr(temp, str_length, str_length)
  if(last_char!="0" & last_char!=" ") {
    return(temp)
  } else{
    temp <- substr(temp, 1, str_length-1)
    return(remove0s(temp))
  }
}

deseqSummary <- function(result, p.val=0.05){
  #print(paste("Number of non converging:", summary(is.na(result$padj))[["TRUE"]]))
  res.ens <- rownames(result)[which(result$padj < p.val)]
  res.ens.up <- rownames(result)[which(result$padj < p.val & result$log2FoldChange > 0)]
  res.ens.dn <- rownames(result)[which(result$padj < p.val & result$log2FoldChange < 0)]
  print(paste("Number of significant genes:", length(res.ens)))
  if(length(res.ens)==0){
    return(NULL)
  }
  #print(paste("Number at uncorrected p:", length(rownames(result)[which(result$pvalue < p.val)])))
  sigGenesEnsembl <- list("ensembl.all"=res.ens, "ensembl.up"=res.ens.up, "ensembl.dn"=res.ens.dn)
  sigGenesSymbols <- list("symbols.all"=ens2Symbol(sigGenesEnsembl$ensembl.all),
                          "symbols.up"=ens2Symbol(sigGenesEnsembl$ensembl.up),
                          "symbols.dn"=ens2Symbol(sigGenesEnsembl$ensembl.dn))
  
  return(list("ensembl"=sigGenesEnsembl, "symbols"=sigGenesSymbols))
}


saveGMX <- function(genes, filename, set.name, set.description="NA"){
  filename <- paste(Sys.Date(), "_", filename, ".gmx", sep="")
  write.table(c(set.name, set.description, genes), quote=F, row.names=F, col.names=F, file = filename)
  
}

ens2Symbol <- function(ensembl){
  require(org.Hs.eg.db)
  x <- org.Hs.egENSEMBL
  xx <- as.list(org.Hs.egENSEMBL2EG)
  
  # now go entrez to gene symbol
  gene_symbols_db <- org.Hs.egSYMBOL
  mapped_symbols <- mappedkeys(gene_symbols_db)
  mapped_symbols <- as.list(gene_symbols_db[mapped_symbols])
  
  inds <- match(ensembl, names(xx))
  inds <- inds[!is.na(inds)]
  
  entrez <- xx[inds] 
  entrez <- unlist(entrez)
  
  symbols <- mapped_symbols[match(entrez, names(mapped_symbols))]
  symbols <- unlist(symbols)
  return(as.character(symbols))
  
}
############################################################################
#### PLAYGROUND FOR ANALYSIS FUNCTION
############################################################################
deseqAnalysis <- function(dds, mod.full, mod.reduced, varOfInt, fdr=0.05, filename="temp", genesets.save=T,
                          enrichr.run=T, enrichr.fdr=0.25,
                          enrichr.dbs=c("Human_Gene_Atlas", "Biocarta_2016",
                                        "KEGG_2016", "WikiPathways_2016",
                                        "Reactome_2016", "GO_Biological_Process_2015",
                                        "GO_Cellular_Component_2015", "GO_Molecular_Function_2015")){
  
  filename1 <- paste(Sys.Date(), "_", filename, ".RDS", sep="")
  
  #  model <- formula(paste("~", paste(covariates, collapse=" + "), sep=""))
  #  m1 <- model.matrix(model, data=colData(dds))
  #  m1 <- m1[ ,-which(apply(unname(m1), 2, sum)==0)]
  #  print(model)
  
  #  model <- formula(paste("~", paste(c(covariates, varOfInt), collapse=" + "), sep=""))
  #  m2 <- model.matrix(model, data=colData(dds))
  #  m2 <- m2[ ,-which(apply(unname(m2), 2, sum)==0)]  
  #  print(model)
  
  #  deseqRes <- DESeq(dds, full=m2, reduced=m1, test="LRT")
  deseqRes <- DESeq(dds, full=mod.full, reduced=mod.reduced, test="LRT")
  
  saveRDS(deseqRes, file = filename1)
  
  res <- results(deseqRes, name=resultsNames(dds)[length(resultsNames(dds))])
  
  sigGenes <- deseqSummary(res, fdr)
  if(!is.null(sigGenes)){
    if(genesets.save){
      saveGMX(genes = sigGenes$symbols$symbols.up,
              filename = paste(filename, varOfInt, "up", sep="_"),
              set.name = paste(filename, varOfInt, "up", sep="_"),
              set.description = paste("Genes up with ", varOfInt, "; below fdr ", as.character(fdr)))
      
      saveGMX(genes = sigGenes$symbols$symbols.dn,
              filename = paste(filename, varOfInt, "down", sep="_"),
              set.name = paste(filename, varOfInt, "up", sep="_"),
              set.description = paste("Genes down with ", varOfInt, "; below fdr ", as.character(fdr)))
      
    }
    
    if(enrichr.run){
      sigGenes.enrich <- enrichFullGeneList(up.genes = sigGenes$symbols$symbols.up,
                                            dn.genes = sigGenes$symbols$symbols.dn,
                                            databases = enrichr.dbs, 
                                            fdr.cutoff = enrichr.fdr)
      
      sigGenes.enrich <- sigGenes.enrich[order(sigGenes.enrich$qval), ]
      
      filename2 <- paste(Sys.Date(), filename, "enrichment.txt", sep="_")
      write.table(sigGenes.enrich, file = filename2, quote=F, row.names=F, col.names=T)
      
    }
  }
}

deseqFinish <- function(dds, varOfInt, fdr=0.05, filename="temp", genesets.save=T,
                        enrichr.run=T, enrichr.fdr=0.25,
                        enrichr.dbs=c("Human_Gene_Atlas", "Biocarta_2016",
                                      "KEGG_2016", "WikiPathways_2016",
                                      "Reactome_2016", "GO_Biological_Process_2015",
                                      "GO_Cellular_Component_2015", "GO_Molecular_Function_2015")){
  filename1 <- paste(Sys.Date(), "_", filename, ".RDS", sep="")
  saveRDS(dds, file = filename1)
  
  res <- results(dds, name=resultsNames(dds)[length(resultsNames(dds))])
  
  sigGenes <- deseqSummary(res, fdr)
  
  if(!is.null(sigGenes)){
    if(genesets.save){
      saveGMX(genes = sigGenes$symbols$symbols.up,
              filename = paste(filename, varOfInt, "up", sep="_"),
              set.name = paste(filename, varOfInt, "up", sep="_"),
              set.description = paste("Genes up with ", varOfInt, "; below fdr ", as.character(fdr)))
      
      saveGMX(genes = sigGenes$symbols$symbols.dn,
              filename = paste(filename, varOfInt, "down", sep="_"),
              set.name = paste(filename, varOfInt, "up", sep="_"),
              set.description = paste("Genes down with ", varOfInt, "; below fdr ", as.character(fdr)))
      
    }
    
    if(enrichr.run){
      sigGenes.enrich <- enrichFullGeneList(up.genes = sigGenes$symbols$symbols.up,
                                            dn.genes = sigGenes$symbols$symbols.dn,
                                            databases = enrichr.dbs, 
                                            fdr.cutoff = enrichr.fdr)
      
      sigGenes.enrich <- sigGenes.enrich[order(sigGenes.enrich$qval), ]
      
      filename2 <- paste(Sys.Date(), filename, "enrichment.txt", sep="_")
      write.table(sigGenes.enrich, file = filename2, quote=F, row.names=F, col.names=T)
      
    }
  }
}

