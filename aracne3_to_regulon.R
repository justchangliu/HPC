aracne3_to_regulon <- function(net.file, anno=NULL, MI.thres=0, regul.size=50){
  net = read.delim(net.file)
  if(!is.null(anno)){
    if(ncol(anno)!=2 | !is.element(class(anno),c("matrix","data.frame","Matrix"))){
      stop("anno should contain two columns: 1-original symbol, 2-new symbol")
    }
    colnames(anno) = c("old","new")
    ## Convert gene symbols:
    net$regulator.values = anno$new[match(net$regulator.values, anno$old)]
    net$target.values = anno$new[match(net$target.values, anno$old)]
  }
  ## network filtering
  net = net[net$mi.values>MI.thres,]
  ## Total MR set
  mr = unique(net$regulator.values)
  regul = lapply(mr, function(mri){
    ## Regulon for MR_i
    tmp.net.data <- net[which(net$regulator.values == mri),]
    ## Sort interactions in the decreasing order of count and MI
    tmp.net.data <- tmp.net.data[order(tmp.net.data$count.values,
                                       tmp.net.data$mi.values, decreasing = TRUE),]
    # print(dim(tmp.net.data))
    ## Top 50 interactions
    tmp.net.data <- tmp.net.data[seq(from = 1, to = min(c(regul.size, nrow(tmp.net.data))), by = 1),]
    ## Regulatory mode = spearman correlation score
    tmp.net.data$am.values <- tmp.net.data$scc.values
    ## Regulatory weight = scaled MI
    tmp.net.data$aw.values <- (tmp.net.data$mi.values/max(tmp.net.data$mi.values))
    tmp.regul <- list(am = as.numeric(tmp.net.data$am.values),
                      aw = (as.numeric(tmp.net.data$aw.values)))
    names(tmp.regul$aw) <- as.character(tmp.net.data$target.values)
    names(tmp.regul$am) <- as.character(tmp.net.data$target.values)
    return(tmp.regul)
  })
  names(regul)=mr
  return(regul)
}
