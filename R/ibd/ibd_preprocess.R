ibd_preprocess <- function(metagenomics_df){
  ## df2 as the metagenomics data
  df2 <- metagenomics_df
  rownames(df2) = df2[,1]
  df2 <- df2[,-1]

  ## filter out features with more zeros
  library(data.table)
  DT<-data.table(df2)
  DT$Total_0s<-rowSums(DT==0)
  zeros <- DT$Total_0s / (length(colnames(DT))-1)
  pathway <- rownames(df2)
  #plot(density(zeros))
  df8 <- DT[zeros <= 0.8,]
  name8 <- pathway[zeros <= 0.8]

  ## add pseudo 1
  df8 <- as.data.frame(df8+1)
  rownames(df8) <- name8
  df8 <- df8[,-1639]

  ## filter-out pathways with low RSD
  library(dplyr)
  library(matrixStats)

  columns <- colnames(df8)
  df8 <- df8 %>%
    mutate(Mean= rowMeans(.[columns]), stdev=rowSds(as.matrix(.[columns])))
  df8$rsd <- df8$stdev/df8$Mean
  boxplot(df8$rsd)
  plot(density(df8$rsd))
  quantile_0.05 = quantile(df8$rsd,0.05)
  df8 <- df8[df8$rsd >= quantile_0.05,]
  df8 <- df8[,!(colnames(df8) %in% c("stdev", "Mean", "rsd"))]

  ## central Log tranformation
  library(compositions)
  df5 <- clr(df8)

  return(list(output = df8))
}
