Mgene_EXPLORER <- function(cpg, methy, express, win)
{
  # annot <- readRDS("/projects/b1042/HouLab/Jingzhe/geneExpress/annot.RDS")
  
  library(biomaRt)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  annot <- subset(annot, Name %in% cpg)
  annot$chr <- sub("^chr", "", annot$chr)
  
  # distance
  annot$start_pos <- annot$pos - win
  annot$end_pos <- annot$pos + win
  
  # filter methylation data by input cpg
  M <- as.data.frame(methy)
  M <- M[rownames(M) %in% cpg, ]
  M <- as.matrix(M)
  
  # load biomart for matching ensembl annotation
  print("Start loading ensembl dataset hsapiens_gene_ensembl")
  # ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
  print("Ensembl dataset hsapiens_gene_ensembl loaded successfully")
  
  # clean ensembl id (need to think about non unique condition)
  ensg_ids = c(sub("\\..*", "", rownames(express)))
  
  # match ensembl id with gene symbol for gene expression data
  geneExp <- as.data.frame(express)
  geneExp$ensembl_gene_id <- ensg_ids
  
  # only keep the first occurrence of each ensembl id after remove dot and number after dot
  geneExp <- geneExp[!duplicated(geneExp$ensembl_gene_id), ]
  
  # get annotation (ensembl gene id, gene symbol, start position, end position)
  results <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position", "strand", "chromosome_name"), 
                   filters = 'ensembl_gene_id', values = ensg_ids, mart = ensembl)
  
  # merge expression data with annotation data
  geneExp <- merge(geneExp, results, all = FALSE)
  
  # keep the ensembl id if no annotated gene symbol
  geneExp$external_gene_name[geneExp$external_gene_name == ""] <- geneExp$ensembl_gene_id[geneExp$external_gene_name == ""]
  
  # creat output data
  output <- list()
  
  # create methylation data and expression data for each cpg
  for (i in 1:length(cpg)) {
    
    # make sure all genes are in same chromosome with this cpg
    expdat <- geneExp[geneExp$chromosome_name == annot$chr[annot$Name == cpg[i]], ]
    
    # filter all neighboring gene expression data related with cpg by distance (old selection)
    # expdat <- expdat[ expdat$start_position >= annot$start_pos[annot$Name == cpg[i]] &  expdat$end_position <= annot$end_pos[annot$Name == cpg[i]], ]
    
    # calculate distance between cpg to gene, if gene strand is positive, distance is from cpg pos to gene start pos, otherwise, to end pos
    expdat$distance <- ifelse(expdat$strand == 1,
                              abs(annot$pos[annot$Name == cpg[i]] - expdat$start_position),
                              abs(annot$pos[annot$Name == cpg[i]] - expdat$end_position))
    
    # filter all neighboring gene expression data related with cpg by distance (new selection)
    expdat <- expdat[expdat$distance <= win,]
    
    # if no genes are within the distance, skip this CpG and move to the next CpG
    if (nrow(expdat) == 0) {
      next
    }
    
    # Remove unnecessary variables
    expdat <- subset(expdat, select = -c(start_position, end_position, strand, chromosome_name))
    
    # temp dataframe for selecting gene with maxium correlation
    temp_res <- data.frame(correlation = numeric(0), p_value = numeric(0), gene = character(0), ensemblID = character(0), distance = numeric(0))
    
    # clean expression data for calculating
    subexpdat <- subset(expdat, select = -c(ensembl_gene_id, external_gene_name, distance))
    subexpdat <- as.matrix(subexpdat)
    
    # calculate correlation between cpg and each gene in expression data
    for (j in 1:nrow(expdat)) {
      # store all correlation and gene symbol for this cpg
      test_result <- cor.test(M[rownames(M) == cpg[i],], subexpdat[j,], method = "pearson", use="complete.obs")
      temp_res[j, "correlation"] <- test_result$estimate
      temp_res[j, "p_value"] <- test_result$p.value
      temp_res[j, "gene"] <- expdat$external_gene_name[j]
      temp_res[j, "ensemblID"] <- expdat$ensembl_gene_id[j]
      temp_res[j, "distance"] <- expdat$distance[j]
    }
    
    # Filter results with absolute correlation > 0.2
    # temp_res <- temp_res[abs(temp_res$correlation) > 0.1, ]
    
    # remove missing that caused by zero expression
    temp_res <- na.omit(temp_res)
    
    # sort by correlation 
    temp_res$abs_correlation <- abs(temp_res$correlation)
    temp_res_sorted <- temp_res[order(-temp_res$abs_correlation), ]
    temp_res_sorted <- temp_res_sorted[, -ncol(temp_res_sorted)]
    
    # if the correlation of this cpg is not empty data, give this data a row with that cpg name
    if (nrow(temp_res_sorted) != 0) {
      temp_res_sorted$cpg <- cpg[i]
    }
    
    # Add results to output list
    output[[cpg[i]]] <- temp_res_sorted
    
  }
  
  # rbind all correlation data
  output <- do.call(rbind, output)
  rownames(output) <- NULL
  return(output)
}