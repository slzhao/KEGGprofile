\name{plot_pathway}
\alias{plot_pathway}
\title{plot_pathway}
\usage{
  plot_pathway(gene_expr, line_col, groups,
    pathway_id = "00010", specis = "hsa", pathway_min = 5,
    database_dir = getwd(), ...)
}
\arguments{
  \item{...}{any other Arguments for function plot_profile}

  \item{gene_expr}{the matrix for gene expression,
  row.names should be NCBI gene ID, such as 67040, 93683}

  \item{line_col}{line color for expression in different
  samples in the pathway map, valid when type='lines'}

  \item{groups}{a character used to indicate expression
  values from different samples}

  \item{specis}{the specis id in KEGG database, 'hsa' means
  human, 'mmu' means mouse, 'rno' means rat, etc}

  \item{pathway_min}{The pathways with number of annotated
  genes less than pathway_min would be ignored}

  \item{database_dir}{the directory where the XML files and
  png files are located}

  \item{pathway_id}{the KEGG pathway id, such as '00010'}
}
\description{
  A wrapper for function download_KEGGfile, parse_XMLfile
  and plot_profile
}
\details{
  This wrapper function is developed to make the
  visualization process more easier. Firstly the existence
  of XML file and png file would be checked, if not, the
  download_KEGGfile function would be used to download the
  files. Then the parse_XMLfile function would be used to
  parse the XML file. At last the plot_profile function
  would be used to generate the pathway map.
}
\examples{
data(pro_pho_expr)
data(pho_sites_count)
#type='lines'
col<-col_by_value(pho_sites_count,col=colorRampPalette(c('white','khaki2'))(4),breaks=c(0,1,4,10,Inf))
temp<-plot_pathway(pro_pho_expr,bg_col=col,line_col=c("brown1","seagreen3"),groups=c(rep("Proteome ",6),rep("Phosphoproteome ",6)),magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110",max_dist=5)
#type='bg'
pho_expr<-pro_pho_expr[,7:12]
temp<-apply(pho_expr,1,function(x) length(which(is.na(x))))
pho_expr<-pho_expr[which(temp==0),]
col<-col_by_value(pho_expr,col=colorRampPalette(c('green','black','red'))(1024),range=c(-6,6))
temp<-plot_pathway(pho_expr,type="bg",bg_col=col,text_col="white",magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110")
}
\seealso{
  \code{\link{download_KEGGfile}},
  \code{\link{parse_XMLfile}}, \code{\link{plot_profile}}
}

