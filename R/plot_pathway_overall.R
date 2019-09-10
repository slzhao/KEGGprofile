##' plot_pathway_overall
##' 
##' The function will plot the correlation distributions for each enriched pathway (result from find_enriched_pathway
##'  function), and then Wilcoxon tests would be used to estimate the significance of values distribution between 
##'  genes in each pathway and all other genes.
##' 
##' 
##' 
##' @inheritParams plot_profile
##' @param gene_values A data.frame or matrix with gene ID as rownames. Each column represent gene value in one condition.
##' @param pathwayNumInFigure a integer specifying max number of (top) pathways in each direction in the boxplot.
##' @param rankByVar a character string specifying variable (sample) name by which to rank the boxplot.
##' @importFrom reshape2 melt
##' @export
##' @return p values for Wilcoxon tests in each pathway
##' @examples data(pro_pho_expr)
##' data(pho_sites_count)
##' gene_values<-pro_pho_expr[row.names(pho_sites_count)[which(pho_sites_count>=10)],]
##' plot_pathway_overall(gene_values=gene_values[,1:3])
plot_pathway_overall<-function(gene_values,species="hsa",pathwayNumInFigure=5,rankByVar=colnames(gene_values)[1]) {
  kegg_enriched_pathway=find_enriched_pathway(row.names(gene_values),species=species,returned_pvalue=1,returned_adjpvalue=1,returned_genenumber=5)
  names(kegg_enriched_pathway[[2]])=make.names(kegg_enriched_pathway[[1]]$Pathway_Name)
  geneValuesInPathway<-lapply(kegg_enriched_pathway[[2]],function(x) gene_values[intersect(x,row.names(gene_values)),,drop=FALSE])

  ##Genes in Pathway vs Genes out of pathway test. Not using them at this time  
#  geneValuesDiffPInPathway=sapply(geneValuesInPathway,function(x) {genesNotInPathway=setdiff(row.names(gene_values),row.names(x));pValue=rep(NA,ncol(x));for (i in 1:ncol(x)) {pValue[i]=wilcox.test(x[,i],gene_values[genesNotInPathway,i])$p.value};return(pValue)})
#  geneValuesDiffPInPathway=t(geneValuesDiffPInPathway)
#  colnames(geneValuesDiffPInPathway)=colnames(gene_values)
  
  dataForPlot=NULL
  for (i in 1:length(geneValuesInPathway)) {
    dataForPlot=rbind(dataForPlot,cbind(geneValuesInPathway[[i]],Pathway=names(geneValuesInPathway)[i]))
  }
  dataForPlot=reshape2::melt(dataForPlot, id.vars="Pathway",variable.name = "Sample",value.name = "value")
  
  #order by median value of first variable (column)
  temp=dataForPlot[which(dataForPlot[,"Sample"]==rankByVar),]
  pathwayOrder=tapply(temp$value,temp$Pathway,median,na.rm=T)
  pathwayOrder=names(pathwayOrder)[order(pathwayOrder)]
  if (length(pathwayOrder)>pathwayNumInFigure*2) {
    pathwayOrder=unique(c(head(pathwayOrder,pathwayNumInFigure),tail(pathwayOrder,pathwayNumInFigure)))
  }
  
  dataForPlot=dataForPlot[which(dataForPlot$Pathway %in% pathwayOrder),]
  dataForPlot$Pathway=factor(dataForPlot$Pathway,levels = pathwayOrder)
  
  p=ggplot(dataForPlot,aes_string(x="Pathway",y="value"))+geom_boxplot(aes_string(colour="Sample"))
  
  return(p+ coord_flip()+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}

