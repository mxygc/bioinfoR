\name{run.david}
\alias{run.david}
\title{submit gene lists to DAVID and save results to local files}
\description{1. term enrichment, 2. term clustering, 3. gene clustering, 4. gene annotation}
\usage{
    run.david <- function(
                          gene.ids, 
                          dir.out = NULL,
                          perl = "perl",
                          auth.email = "vishalrp@uci.edu", 
                          tasks = c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot"), 
                          species = "Homo sapiens", 
                          id.types = c("AFFYMETRIX_3PRIME_IVT_ID"), 
                          list.names = NULL, 
                          categories = c("OMIM_DISEASE", "COG_ONTOLOGY", "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE", "GOTERM_BP_FAT", "GOTERM_CC_FAT", "GOTERM_MF_FAT", "BBID", "BIOCARTA", "KEGG_PATHWAY", "INTERPRO", "PIR_SUPERFAMILY", "SMART"),
                          te.count = 2, te.ease = 0.1,
                          tc.similarity.gene.overlap = 3,
                          tc.initial.group.member = 3,
                          tc.final.group.member = 3,
                          tc.multi.linkage = 0.5,
                          tc.similarity = 0.5,
                          gc.similarity.term.overlap = 4,
                          gc.initial.group.member = 4,
                          gc.final.group.member = 4,
                          gc.multi.linkage = 0.5,
                          gc.similarity = 0.35,
                          verbose = FALSE,
                          envir = parent.frame())
}
\arguments{
  \item{gene.ids}{A character vector (single ID list) or a list of characters (multiple ID lists); when it is a unnamed list, names can be specified by \code{list.names}.}
  \item{dir.out}{Where the output is stored. Default the current workspace.}
  \item{perl}{Where perl is installed. Default \code{perl} defined by your $PATH is used}
  \item{auth.email}{An autheticated email address for DAVID web service. Default \code{vishalrp@uci.edu} is used but you'd better register and use your own, to avoid overusing the default account.}
  \item{tasks}{A vector of characters from the 4 tasks \code{TermEnrich}, \code{TermClust}, \code{GeneClust}, \code{GeneAnnot}.}
  \item{species}{Species of your IDs (e.g. \code{Rattus norvegicus}, \code{Mus musculus}, etc. Default: \code{Homo sapiens}.}
  \item{id.types}{A single or a vector (multiple ID lists) of characters specifying the ID types. (e.g. \code{AFFYMETRIX_3PRIME_IVT_ID}, \code{REFSEQ_MRNA}, \code{ENSEMBL_GENE_ID}, \code{ENTREZ_GENE_ID}, etc. Default: \code{AFFYMETRIX_3PRIME_IVT_ID}.}  
  \item{list.names}{The name(s) of input list(s). If provided, the length of list.names should match the length of the input list. Default: NULL}
  \item{categories}{A vector of characters specifying the pathways/databases for functional analyese. Default: \code{c("OMIM_DISEASE", "COG_ONTOLOGY", "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE", "GOTERM_BP_FAT", "GOTERM_CC_FAT", "GOTERM_MF_FAT", "BBID", "BIOCARTA", "KEGG_PATHWAY", "INTERPRO", "PIR_SUPERFAMILY", "SMART")}.}
  \item{te.count}{A positive integer specifying Term Enrichment \code{Count Threshold} (see \url{http://david.abcc.ncifcrf.gov/content.jsp?file=functional_annotation.html#E3}). Default: 2}
  \item{te.ease}{A positive float specifying Enrichment \code{EASE Threshold} (see \url{http://david.abcc.ncifcrf.gov/content.jsp?file=functional_annotation.html#E3}). Default: 0.1}
  \item{tc.similarity.gene.overlap}{A positive float specifying Term Clustering \code{Similarity Gene Overlap}. Default: 3} 
  \item{tc.initial.group.member}{A positive integer specifying Term Clustering \code{Initial Group Membership}. Default: 3} 
  \item{tc.final.group.member}{A positive integer specifying Term Clustering \code{Final Group Membership}. Default: 3}
  \item{tc.multi.linkage}{A positive float specifying Term Clustering \code{Multiple Linkage Threshold}. Default: 0.5} 
  \item{tc.similarity}{A positive float specifying Term Clustering \code{Similarity Threshold}. Default: 0.5}
  \item{gc.similarity.gene.overlap}{A positive float specifying Gene Clustering \code{Similarity Term Overlap}. Default: 4} 
  \item{gc.initial.group.member}{A positive integer specifying Geme Clustering \code{Initial Group Membership}. Default: 4} 
  \item{gc.final.group.member}{A positive integer specifying Gene Clustering \code{Final Group Membership}. Default: 4}
  \item{gc.multi.linkage}{A positive float specifying Gene Clustering \code{Multiple Linkage Threshold}. Default: 0.5} 
  \item{gc.similarity}{A positive float specifying Gene Clustering \code{Similarity Threshold}. Default: 0.35}



}
\details{ }
\value{ No value returned}
\note{ 
    1. A maximum of 200 lists can be submitted per day per IP; 
    2. Each list should not exceed 3000 genes; 
    3. ID type \code{OFFICIAL_GENE_SYMBOL} not supported by DAVID due its ambiguity (see \url{http://david.abcc.ncifcrf.gov/forum/viewtopic.php?f=14&t=885&sid=d860a473f348c9bcd57e669ddf308fee}).
}
\author{ The author's name of the Perl core scripts not found on DAVID. :( }
\references{
  DAVID Web Service Perl Client \url{http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html
  \url{http://david.abcc.ncifcrf.gov/webservice/sample_clients/PerlClient-1.1.zip}}
 }
\examples{
    data("ids")
    run.david(gene.ids = ids, dir.out = "results", tasks = c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot"), species = "Homo sapiens", id.types = c("AFFYMETRIX_3PRIME_IVT_ID", "REFSEQ_MRNA", "ENSEMBL_GENE_ID", "ENTREZ_GENE_ID"))
}
