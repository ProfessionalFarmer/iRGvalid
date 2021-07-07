library(dplyr)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 清洗上一年的数据，共shiny使用
load("/data/home2/Zhongxu/work/from35/housek/20200922-all.combination.rdata")

all.projects.comb = reshape2::dcast(all.projects.comb, V1+V3~Gene, value.var="V2")
colnames(all.projects.comb)[1:2] <- c("Panel size","Panel")
all.projects.comb[,3:ncol(all.projects.comb)] = loonR::convertDfToNumeric(all.projects.comb[,3:ncol(all.projects.comb)] )

save(all.projects.comb, file="/data/home2/Zhongxu/work/iRGvalid/data/all.combination.rdata")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 过滤paried sample中不表达的，节省空间，并且转换成基因名
for(nf in list.files(path = "/data/home2/Zhongxu/work/from35/housek/TCGA", 
                     full.names = FALSE, pattern = "TPM.paired.N.rdata") ){
  
  cat(nf,"\n")
  project = stringr::str_remove_all(nf,".RNA-TPM.paired.N.rdata")
  cat(project,"\n")
  
  load(paste0("/data/home2/Zhongxu/work/from35/housek/TCGA/",project,".RNA-TPM.paired.N.rdata"))
  load(paste0("/data/home2/Zhongxu/work/from35/housek/TCGA/",project,".RNA-TPM.paired.T.rdata"))
  
  project.df <- merge(normal.df, tumor.df, by=0, all=TRUE)
  row.names(project.df) <- project.df$Row.names
  project.df = project.df[,-c(1)]
  
  project.df = CancerSubtypes::FSbyVar(project.df, cut.type = "cutoff", value = "0.1")
  cat(dim(project.df))
  
  id.mapping <- loonR::id_mapping(row.names(project.df), key = "ENSEMBL", column = c("SYMBOL")  )
  
  
  # 如果没有map到或者重复则删除
  id.mapping <- id.mapping[ !duplicated(id.mapping$SYMBOL), ] 
  id.mapping <- id.mapping[ !is.na(id.mapping$SYMBOL), ] 
  
  project.df <- project.df[id.mapping$Ref, ]
  rownames(project.df) = id.mapping$SYMBOL
  
  save(project.df, file=paste0("/data/home2/Zhongxu/work/iRGvalid/data/", project, ".TPM.paired.T.N.rdata")  )
  
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


getGeneExpressionDf <- function(project, candidate, referenece.genes){
  
  project = "TCGA-BRCA"
  candidate = "HLA-A"
  referenece.genes = "CIAO1\nCNBP\nHEY1\nUBC"
  
  load(paste0("Data/", project,".TPM.paired.T.N.rdata"))
  referenece.genes = unlist( strsplit(referenece.genes,"\n") )
  
  tmp.genes <- c(candidate, referenece.genes)
  
  
  candidate.df <- data.frame( t(project.df[tmp.genes,]), check.names = FALSE )
  paired.num <- ncol(project.df)/2
  candidate.df$Group <- rep(c("N", "T"), c(paired.num, paired.num))
  candidate.df
  
}





