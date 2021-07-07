getGeneExpressionDf <- function(project, candidate, referenece.genes){
  
  # project = "TCGA-BRCA"
  # target.gene = "HLA-A"
  # referenece.genes = c("CIAO1", "CNBP",  "HEY1",  "UBC")
  
  load(paste0(project,".TPM.paired.T.N.rdata"))
  referenece.genes = unlist( strsplit(referenece.genes,"\n") )
  
  tmp.genes <- c(candidate, referenece.genes)
  
  
  candidate.df <- data.frame( t(project.df[tmp.genes,]), check.names = FALSE )
  paired.num <- ncol(project.df)/2
  candidate.df$Group <- rep(c("N", "T"), c(paired.num, paired.num))
  candidate.df
  
}


normalize <- function( gene.exp.df, referenece.genes, target.gene ){
  res = list()
  
  # target.gene的原始数据
  target.gene.exp <- unlist( gene.exp.df[,target.gene] )
  combination.genes.exp <- data.frame(gene.exp.df[,referenece.genes])
  
  # normalize
  post.target.gene.exp <- target.gene.exp - rowMeans(combination.genes.exp)
  
  # normalize前后的相关性
  correlation <- round(cor(unlist(post.target.gene.exp), target.gene.exp, method = "pearson"),3)
  
  target.df <- data.frame(Sample = rownames(gene.exp.df), check.names = F,
                          `Before normalization` = target.gene.exp,
                          `After normalization` = post.target.gene.exp)

  
  res$correlation = correlation
  res$target.df = target.df
  
  res
  
}


normalize.all.combination <- function( gene.exp.df, target.gene, referenece.genes ){
  
  library(foreach)

  total.panel.mem <- referenece.genes

  # target.gene的原始数据
  target.gene.exp <- unlist( gene.exp.df[,target.gene] )
  
  # 从1个到所有gene拿来normalize
  all.res <- foreach::foreach(panel.size=1:length(total.panel.mem), .combine = rbind) %do% {
    
    size.combination <- gtools::combinations(length(total.panel.mem), panel.size, total.panel.mem )

    # 遍历特定size下的所有组合
    com.res <- foreach(row.ind=1:nrow(size.combination), .combine = rbind ) %do% {
      
      combination.genes <- size.combination[row.ind,]
      
      # 转成string 用于保存
      combination.genes.string <- paste(combination.genes,sep = '', collapse = ", ")
      
      combination.genes.exp <- data.frame(gene.exp.df[,combination.genes])
      
      # normalize前后的相关性
      post.target.gene.exp <- target.gene.exp - rowMeans(combination.genes.exp)
      correlation <- round(cor(post.target.gene.exp, target.gene.exp, method = "pearson"),3)
      
      c(panel.size, correlation, combination.genes.string )
    }
    com.res
  }

  all.res <- as.data.frame(all.res, stringsAsFactors = FALSE)
  all.res$V2 <- as.numeric(all.res$V2)
  all.res <- all.res[ order(all.res$V2,decreasing = T), ]
  colnames(all.res) <- c("Panel size", "Correlation", "Panel")
  
  all.res
  
}
  
  
