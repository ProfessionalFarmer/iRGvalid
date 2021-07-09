getGeneExpressionDf <- function(project, candidate, referenece.genes){
  
  # project = "TCGA-BRCA"
  # target.gene = "HLA-A"
  # referenece.genes = c("CIAO1", "CNBP",  "HEY1",  "UBC")

  load(paste0(project,".TPM.paired.T.N.rdata"))
  
  if(  length(setdiff(referenece.genes, row.names(project.df)))!=0   ){
    na.gene = setdiff(referenece.genes, row.names(project.df) )
    showNotification(paste0(na.gene," may have low expression") )
    
    referenece.genes = setdiff(referenece.genes, na.gene)
  }
    
  tmp.genes <- c(candidate, referenece.genes)
  
  
  candidate.df <- data.frame( t(project.df[tmp.genes,]), check.names = FALSE )
  paired.num <- ncol(project.df)/2
  candidate.df$Group <- rep(c("N", "T"), c(paired.num, paired.num))
  candidate.df
  
}


normalize <- function( gene.exp.df, target.gene, referenece.genes){
  res = list()
  
  
  # 最后一列为Group，对 target.gene normalize之前计算correlation
  library(corrplot)
  res$correlation.before <- round(cor(gene.exp.df[,-c(ncol(gene.exp.df))], method = c("pearson")  ), 3) # spearman pearson


  ###############################分别对target.gene进行
  target.gene.exp <- unlist( gene.exp.df[,target.gene] )
  post.target.gene.exp <- target.gene.exp - unlist( rowMeans(gene.exp.df[,referenece.genes]) )
  
  res$normalize.correlation.value <- round(cor(post.target.gene.exp, target.gene.exp, method = "pearson"),3)
  
  
  # normalize 前后的数据
  res$target.df <- data.frame(Sample = rownames(gene.exp.df), check.names = F,
                          `Before normalization` = target.gene.exp,
                          `After normalization` = post.target.gene.exp)
  
  res$scatter = ggpubr::ggscatter(res$target.df, x= "Before normalization", y = "After normalization")
  
  
  ########################################## 以上是correlation的分析，以下是组合分析
  library(foreach)
  
  # 从1个到所有gene拿来normalize
  all.res <- foreach::foreach(panel.size=1:length(referenece.genes), .combine = rbind) %do% {
    
    size.combination <- gtools::combinations(length(referenece.genes), panel.size, referenece.genes )
    
    
    # 遍历特定size下的所有组合
    com.res <- foreach(row.ind=1:nrow(size.combination), .combine = rbind ) %do% {
      combination.genes <- size.combination[row.ind,]
      # 转成string 用于保存
      combination.genes.string <- paste(combination.genes,sep = '', collapse = ", ")
      
      combination.genes.exp <- data.frame(gene.exp.df[,combination.genes])
      
      # normalize前后的相关性
      post.target.gene.exp <- target.gene.exp - rowMeans(combination.genes.exp)
      correlation <- round(cor(post.target.gene.exp, target.gene.exp, method = "pearson"),3)
      
      
      c(panel.size, correlation, combination.genes.string)
    }
    com.res
  }
  
  
  all.res <- as.data.frame(all.res, stringsAsFactors = FALSE)
  all.res <- all.res[ order(all.res$V2,decreasing = T), ]
  
  # all.res$Gene <- paste(c(target.gene, project), sep="", collapse = "-")
  colnames(all.res) =  c("Panel size", "Correlation Rt value", "Panel")
  
  res$all.combination = all.res
  
  
  res

}


parseReferenceGene <- function(reference.gene){

  #reference.gene = c("A\n B,\nC, C\nD")
  reference.gene = stringr::str_replace_all(reference.gene, "\n",",")
  reference.gene = stringr::str_replace_all(reference.gene, ",,",",")
  reference.gene = stringr::str_remove_all(reference.gene, " ")
  reference.gene = unique(unlist(strsplit(reference.gene, ",|\\n")))
  
  reference.gene
  
}




