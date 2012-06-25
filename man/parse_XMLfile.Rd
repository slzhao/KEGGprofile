\name{parse_XMLfile}
\alias{parse_XMLfile}
\title{parse_XMLfile}
\usage{
  parse_XMLfile(pathway_id, specis, database_dir = getwd())
}
\arguments{
  \item{database_dir}{the directory where the XML files and
  png files are located}

  \item{pathway_id}{the KEGG pathway id, such as '00010'}

  \item{specis}{the specis id in KEGG database, 'hsa' means
  human, 'mmu' means mouse, 'rno' means rat, etc}
}
\value{
  a matrix containing genes in this pathway, and their
  names, locations etc, which could be used in the function
  plot_profile as param KEGG_database
}
\description{
  The function parses KEGG XML (KGML) files
}
\details{
  This function will parse the KEGG XML (KGML) file. Then a
  matrix with genes in this pathway and related infomations
  will be returned. This matrix can be used for plot the
  expression profiles on the pathway figure.
}
\examples{
XML2database<-parse_XMLfile(pathway_id="04110",specis="hsa",database_dir=system.file("extdata",package="KEGGprofile"))
}

