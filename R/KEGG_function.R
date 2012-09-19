###############################################################################
# Data: 2012-09-19
# Author: Shilin Zhao (zhaoshilin@gmail.com)
###############################################################################

#' @importMethodsFrom AnnotationDbi as.list
#' @name pho_sites_count
#' @title number of phosphorylation sites quantified for each gene
#' @description This data set is a data.frame with number of phosphorylation sites quantified for each gene in the analysis.
#' @docType data
#' @usage pho_sites_count
#' @source Olsen, J.V., et al. (2010) Quantitative phosphoproteomics reveals 
#' widespread full phosphorylation site occupancy during mitosis, Sci Signal, 3, ra3.
NULL
#' @name pro_pho_expr
#' @title expression profiles in proteome and phosphoproteome
#' @description This data set is from a previously published data of proteome and phosphoproteome analysis in different cell phase.
#' The column 1-6 are proteome data and column 7-12 are phosphoproteome data in this data.frame. The 6 time points are G1, G1/S, Early S, Late S, G2, Mitosis.
#' @docType data
#' @usage pro_pho_expr
#' @source Olsen, J.V., et al. (2010) Quantitative phosphoproteomics reveals 
#' widespread full phosphorylation site occupancy during mitosis, Sci Signal, 3, ra3.
NULL
##' download_KEGGfile
##' 
##' The function download XML files and png files from KEGG website to local disk
##' 
##' If pathway_id is set as 'all', all KEGG pathway ids in KEGG.db package will be used and downloaded from KEGG website
##' 
##' @param pathway_id the KEGG pathway id, such as '00010'
##' @param specis the specis id in KEGG database, 'hsa' means human, 'mmu' means mouse, 'rno' means rat, etc
##' @param target_dir the local directory where the downloaded files are saved
##' @export
##' @examples download_KEGGfile(pathway_id="00010",specis='hsa')
download_KEGGfile<-function(pathway_id="00010",specis='hsa',target_dir=getwd()) {
	if (pathway_id=='all') { #download files for all pathway baesd on KEGG.db package
		require(KEGG.db)
		keggpathway2gene <- as.list(KEGG.db::KEGGPATHID2EXTID)
		pathway_id<-names(keggpathway2gene)[grep(specis,names(keggpathway2gene))]
	} else {
		pathway_id<-paste(specis,pathway_id,sep="")
	}
	for (x in 1:length(pathway_id)) {
		print (paste("Downloading files: ",x,"/",length(pathway_id),sep=""))
		download.file(paste("http://www.genome.jp/kegg-bin/download?entry=",pathway_id[x],'&format=kgml',sep=""),paste(target_dir,"/",pathway_id[x],".xml",sep=""))
		download.file(paste("http://www.genome.jp/kegg/pathway/",specis,"/",pathway_id[x],'.png',sep=""),paste(target_dir,"/",pathway_id[x],".png",sep=""),mode="wb")
	}
}

##' parse_XMLfile
##' 
##' The function parses KEGG XML (KGML) files
##' 
##' This function will parse the KEGG XML (KGML) file. Then a matrix with genes in this pathway and related infomations will be returned. This matrix can be used for plot the expression profiles on the pathway figure.
##' 
##' @inheritParams download_KEGGfile
##' @param database_dir the directory where the XML files and png files are located
##' @export
##' @return a matrix containing genes in this pathway, and their names, locations etc, which could be used in the function plot_profile as param KEGG_database
##' @examples XML2database<-parse_XMLfile(pathway_id="04110",specis="hsa",database_dir=system.file("extdata",package="KEGGprofile"))
parse_XMLfile<-function(pathway_id,specis,database_dir=getwd()) {
	#get pathway gene and their location in the pic
	#except three global maps
	if (pathway_id=="01100" | pathway_id=="01110" | pathway_id=="01120") {
		print(paste("Skip global maps:",pathway_id,sep=""))
		return(NULL)
	}
	require(XML)
	inter1<-xmlTreeParse(paste(database_dir,"/",specis,pathway_id,".xml",sep=""),useInternal=TRUE)
	inter2<-getNodeSet(inter1,"//entry")
	inter3<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "name"))
	inter4<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "type"))
	inter5<-sapply(inter2,  function(xxx) getNodeSet(xxx,".//graphics"))
	inter_graphic_type<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "type"))
	inter6<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "x"))
	inter7<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "y"))
	inter8<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "width"))
	inter9<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "height"))
	inter10<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "name"))
	result<-NULL
	for (i in 1:length(inter4)) {
		if (inter4[[i]]=="gene" & inter_graphic_type[[i]]!="line") {
			temp<-strsplit(inter3[[i]]," ")[[1]]
			name<-strsplit(inter10[[i]],",")[[1]][1]
			name<-gsub('\\.\\.\\.',"",name)
			for (n in 1:length(temp)) {
				result<-rbind(result,c(temp[n],inter6[[i]],inter7[[i]],inter8[[i]],inter9[[i]],name))
			}
		}
	}
	result[,1]<-gsub(paste(specis,":",sep=""),'',result[,1])
	return(result)
}

##' plot_profile
##' 
##' The function plot gene expression profiles on KEGG pathway maps
##' 
##' There are two visualization methods to represent gene expression profiles: 'background' and 'lines'.
##' The first one is applicable for analysis with only one sample or one type of data, which divides the gene polygon into several sub-polygons to represent different time points.
##' And each sub-polygon has a specific background color to represent expression changes in that time point. The second method plots lines with different colors in the gene polygon to represent different samples or different types of data.
##' The dynamic changes of lines mean the profiles of genes in different time points.
##' 
##' @param gene_expr the matrix for gene expression, row.names should be NCBI gene ID, such as 67040, 93683
##' @param pathway_name the specis id and KEGG pathway id, such as 'hsa00010'
##' @param KEGG_database the matrix returned by function parse_XMLfile, which contains genes in this pathway, and their names, locations etc
##' @param groups a character used to indicate expression values from different samples
##' @param bg_col background color for gene rectangles in the pathway map
##' @param line_col line color for expression in different samples in the pathway map, valid when type='lines'
##' @param text_col the colors for text in the pathway map. A color matrix generated by function \code{\link{col_by_value}} can be used here
##' @param border_col border color for gene rectangles in the pathway map. A color matrix generated by function \code{\link{col_by_value}} can be used here
##' @param text_cex cex for text in the pathway map. A color matrix generated by function \code{\link{col_by_value}} can be used here
##' @param magnify the coefficient used to magnify the gene rectangles
##' @param type the type of pathway map visulization, could be 'bg' or 'lines'. Default is 'bg'. See also 'Details'
##' @param pathway_min The pathways with number of annotated genes less than pathway_min would be ignored
##' @param genes_kept methods used for choosing genes when several genes corresponding to one location in pathway map. Default is 'foldchange', which kept the gene with largest fold changes. 'first' kept the first gene. 'random' chosed gene random. 'var' kept the gene with largest variation. 'abs' kept the gene with largest absolute value
##' @param max_dist The expression changes that represented by the distance from the bottom to the top of gene rectangle, valid when type='lines'. This param is used to ensure the dynamic changes of lines in different gene polygon represent equal variation. It would be calculated from the maximum changes of genes in this pathway by default. If max_dist=NA, then the lines would be plotted from top to bottom in each gene rectangle
##' @param lwd The line width when type='lines' 
##' @inheritParams parse_XMLfile
##' @export
##' @return a matrix containing genes maped in this pathway, and their names, expressions 
##' @examples XML2database<-parse_XMLfile(pathway_id="04110",specis="hsa",database_dir=system.file("extdata",package="KEGGprofile"))
##' data(pro_pho_expr)
##' temp<-plot_profile(pro_pho_expr,pathway_name="hsa04110",KEGG_database=XML2database,line_col=c("brown1","seagreen3"),groups=c(rep("Proteome ",6),rep("Phosphoproteome ",6)),magnify=1.2,database_dir=system.file("extdata",package="KEGGprofile"),max_dist=5)
plot_profile<-function(gene_expr,pathway_name,KEGG_database,groups,bg_col="white",text_col="black",line_col,border_col="grey",text_cex=0.25,magnify=1,type=c('lines','bg'),pathway_min=5,genes_kept=c('foldchange','first','random','var','abs'),specis='hsa',database_dir=getwd(),max_dist,lwd=1.2) {
	type <- if (missing(type)) 
				"lines" else match.arg(type)
	if (type == "lines" & ncol(gene_expr)<=1) {
		print("When type=='lines', You should have more than one time points")
	}
	genes_kept <- if (missing(genes_kept)) 
				"foldchange" else match.arg(genes_kept)
	if (missing(groups)) groups<-rep(1,ncol(gene_expr))
	groups<-factor(groups,levels=unique(groups))
	if (missing(line_col)) line_col<-rainbow(length(unique(groups)))
	if (is.matrix(bg_col) | is.data.frame(bg_col)) {} else {
		bg_col<-matrix(rep(bg_col,nrow(gene_expr)),ncol=1)
		row.names(bg_col)<-row.names(gene_expr)
	}
	if (is.matrix(border_col) | is.data.frame(border_col)) {} else {
		border_col<-matrix(rep(border_col,nrow(gene_expr)),ncol=1)
		row.names(border_col)<-row.names(gene_expr)
	}
	if (is.matrix(text_col) | is.data.frame(text_col)) {} else {
		text_col<-matrix(rep(text_col,nrow(gene_expr)),ncol=1)
		row.names(text_col)<-row.names(gene_expr)
	}
	return_expr<-NULL
	
	#the number of annotated genes
	genes<-intersect(KEGG_database[,1],row.names(gene_expr))
	if (length(genes)<pathway_min) {
		print (paste("The genes mapped in pathway ",pathway_name," were less than ",pathway_min,", skip this pathway.",sep=""))
		return()
	}
	
	#plot KEGG pic
	require(png)
	require(TeachingDemos)
	img  <-  readPNG(paste(database_dir,"/",pathway_name,".png",sep=""))
	width<-ncol(img)
	height<-nrow(img)
	err_x_location<-1
	err_y_location<-0
	png(paste(pathway_name,"_profile_",type,".png",sep=""),width=width*2,height=height*2,res=600)
	par(yaxs="i")
	par(xaxs="i")
	par(mar=c(0,0,0,0))
	plot(c(0,width),c(0,height),  type='n',xlab="",ylab="")
	rasterImage(img,  0,  0,  width,  height,interpolate=F)
	
	#plot gene profile in KEGG pic
	result_genes<-as.data.frame(KEGG_database[which(KEGG_database[,1] %in% genes),],stringsAsFactors=F)
	colnames(result_genes)<-c("genes","x","y","width","height","name")
	result_genes<-transform(result_genes, x = as.numeric(x), y = as.numeric(y),width = as.numeric(width),height = as.numeric(height))
	if (missing(max_dist) & type == "lines") {
		max_dist<-max(apply(matrix(gene_expr[result_genes[,1],]),1,function(x) range(x,na.rm=T)[2]-range(x,na.rm=T)[1]))
	}
	findUnique<-apply(result_genes[,2:3],1,function(x) paste(x,collapse=" "))
	temp<-split(as.data.frame(result_genes,stringsAsFactors=F),findUnique)
	for (xx in 1:length(temp)) {
		#for genes in same location, only keep one
		if (length(temp[[xx]][,1])>1) {
			if (genes_kept=="foldchange") {
				ChosedGene<-temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][,1],]),1,function(x) range(x,na.rm=T)[2]-range(x,na.rm=T)[1])),1]
			} else if (genes_kept=="first") {
				ChosedGene<-temp[[xx]][1,1]
			} else if (genes_kept=="random") {
				ChosedGene<-sample(temp[[xx]][,1],1)
			} else if (genes_kept=="var") {
				ChosedGene<-temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][,1],]),1,var)),1]
			} else if (genes_kept=="abs") {
				ChosedGene<-temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][,1],]),1,function(x) max(abs(x)))),1]
			}
			ChosedGene<-temp[[xx]][which(temp[[xx]][,1]==ChosedGene),]
		} else {ChosedGene<-temp[[xx]]}	
		i_max<-ncol(bg_col)
		for (i in 1:i_max) {
			plot_polygon(ChosedGene=ChosedGene,i=i,i_max=i_max,col=bg_col,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=magnify)
		}
		if (type == "lines") {
			ChosedGeneProfile<-as.matrix(gene_expr[as.character(ChosedGene[,1]),])
			ChosedGeneProfile<-sapply(split(ChosedGeneProfile,groups),function(x) x)		
			gene_dist<-gene_expr[as.character(ChosedGene[,1]),]
			gene_dist<-range(gene_dist,na.rm=T)[2]-range(gene_dist,na.rm=T)[1]
			if (is.na(max_dist) | gene_dist==0) {y_ratio<-1} else {
				y_ratio<-gene_dist/max_dist
				if (y_ratio>1) {y_ratio<-1}
			}
			old_par <- par(no.readonly = TRUE)
			subplot(matplot(ChosedGeneProfile,main="",type="l",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",col=line_col,lty=1,lwd=lwd),c(ChosedGene[,2]-ChosedGene[,4]/2*magnify+err_x_location,ChosedGene[,2]+ChosedGene[,4]*magnify/2+err_x_location),c(height-ChosedGene[,3]-ChosedGene[,5]*y_ratio/2*magnify-err_y_location,height-ChosedGene[,3]+ChosedGene[,5]*y_ratio/2*magnify-err_y_location))
			par(old_par)
		}
		plot_polygon(ChosedGene=ChosedGene,col=as.data.frame(NULL),i=1,i_max=1,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=magnify,border=border_col[as.character(ChosedGene[,1]),1],lty=1,lwd=0.5)

		text(ChosedGene[,2],height-ChosedGene[,3],labels=ChosedGene[,6],cex=text_cex,col=text_col[as.character(ChosedGene[,1]),1])
		return_expr<-rbind(return_expr,c(ChosedGene[,1],ChosedGene[,6],gene_expr[as.character(ChosedGene[,1]),]))
	}
	if (type == "lines") {
		legend("topright",legend=unique(groups),lwd=1,col=line_col,bty="n",cex=0.3)
		if (!(is.na(max_dist))) {
			polygon(c(width-66,width-20,width-20,width-66),c(height-60,height-60,height-77,height-77),border="grey")
			text(width-43,height-68.5,round(max_dist,2),cex=text_cex)
		}
	}
	dev.off()
	print (paste("The figure was generated in ",getwd(),"/",pathway_name,"_profile_",type,".png",sep=""))
	colnames(return_expr)<-c("Gene","Name",colnames(gene_expr))
	return(unique(return_expr))
}

##' plot_pathway
##' 
##' A wrapper for function download_KEGGfile, parse_XMLfile and plot_profile
##' 
##' This wrapper function is developed to make the visualization process more easier.
##' Firstly the existence of XML file and png file would be checked, if not, the download_KEGGfile function would be used to download the files.
##' Then the parse_XMLfile function would be used to parse the XML file.
##' At last the plot_profile function would be used to generate the pathway map.
##' 
##' @inheritParams plot_profile
##' @inheritParams parse_XMLfile
##' @param ... any other Arguments for function plot_profile
##' @export
##' @seealso \code{\link{download_KEGGfile}}, \code{\link{parse_XMLfile}}, \code{\link{plot_profile}}
##' @examples data(pro_pho_expr)
##' data(pho_sites_count)
##' #type='lines'
##' col<-col_by_value(pho_sites_count,col=colorRampPalette(c('white','khaki2'))(4),breaks=c(0,1,4,10,Inf))
##' temp<-plot_pathway(pro_pho_expr,bg_col=col,line_col=c("brown1","seagreen3"),groups=c(rep("Proteome ",6),rep("Phosphoproteome ",6)),magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110",max_dist=5)
##' #type='bg'
##' pho_expr<-pro_pho_expr[,7:12]
##' temp<-apply(pho_expr,1,function(x) length(which(is.na(x))))
##' pho_expr<-pho_expr[which(temp==0),]
##' col<-col_by_value(pho_expr,col=colorRampPalette(c('green','black','red'))(1024),range=c(-6,6))
##' temp<-plot_pathway(pho_expr,type="bg",bg_col=col,text_col="white",magnify=1.2,specis='hsa',database_dir=system.file("extdata",package="KEGGprofile"),pathway_id="04110")
plot_pathway<-function(gene_expr,line_col,groups,pathway_id="00010",specis="hsa",pathway_min=5,database_dir=getwd(),...) {
	if ((!file.exists(paste(database_dir,"/",specis,pathway_id,".xml",sep=""))) | (!file.exists(paste(database_dir,"/",specis,pathway_id,".png",sep="")))) {download_KEGGfile(pathway_id=pathway_id,specis=specis,target_dir=database_dir)}
	XML2data<-parse_XMLfile(pathway_id=pathway_id,specis=specis,database_dir=database_dir)
	if (is.null(XML2data)) {return()}
	return_expr<-plot_profile(gene_expr=gene_expr,KEGG_database=XML2data,groups=groups,line_col=line_col,pathway_name=paste(specis,pathway_id,sep=""),database_dir=database_dir,...)
	return(return_expr)
}

plot_polygon<-function(ChosedGene,i,i_max,col,height,err_x_location,err_y_location,magnify,border=F,...) {
	polygon(c(ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*(i-1)*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*i*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*i*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*(i-1)*magnify+err_x_location),
			c(height-ChosedGene[,3]-ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]-ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]+ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]+ChosedGene[,5]/2*magnify-err_y_location)
			,col=col[as.character(ChosedGene[,1]),i],border=border,...)
}

##' col_by_value
##' 
##' The function will transfer a numeric matrix into a matrix of colors, in which the colors represent the values of numeric matrix
##' 
##' A colorbar would also be ploted. The returned colors of the function can be used in function plot_profile. 
##' if breaks not equal to NA, col must have the same length with breaks-1. 
##' 
##' @param x a numeric matrix
##' @param col colors used to represent the values. (See also 'Details')
##' @param range values out of the range will be modified to in the range.
##' @param breaks a numeric vector of three or more cut points giving the number of intervals into which x is to be cut. See also 'Details'
##' @export
##' @return a matrix equal to x, but the values were instead by colors.
##' @examples data(pho_sites_count)
##' col<-col_by_value(pho_sites_count,col=colorRampPalette(c('white','khaki2'))(4),breaks=c(0,1,4,10,Inf))
col_by_value<-function(x,col,range=NA,breaks=NA) {
	if (is.na(range[1])) {} else {
		x[x<range[1]]<-range[1]
		x[x>range[2]]<-range[2]
	}
	if (is.na(breaks[1])) {
		ff <- seq(min(x),max(x), length=length(col))
		bg2<-apply(as.matrix(as.numeric(unlist(x))),1,function(x) rank(c(ff,x),ties.method ="min")[length(col)+1])
		dens <-matrix(bg2,nrow(x),ncol(x))
		result<-matrix(col[dens],nrow=nrow(x),ncol=ncol(x))
		row.names(result)<-row.names(x)
		image(x=1,y=as.matrix(ff),z=t(ff),col=col,xaxt="n",ylab="")
		box()
		return(result)
	} else {
		temp<-cut(as.numeric(unlist(x)),breaks=breaks,include.lowest=T)
		if (length(col)!=length(levels(temp))) {stop("length:col != length: cut result")}
		result<-matrix(col[as.numeric(temp)],nrow=nrow(x),ncol=ncol(x))
		row.names(result)<-row.names(x)
		par(mar=c(5,9,2,2))
		image(x=1,y=as.matrix(1:(length(breaks)-1)),z=t(1:(length(breaks)-1)),col=col,xaxt="n",yaxt="n",ylab="",xlab="")
		axis(2, at = 1:(length(breaks)-1),labels=levels(temp),las=1)
		box()
		return(result)
	}
}

##' find_enriched_pathway
##' 
##' The function will map the genes in KEGG pathway database, and then hypergegeometric tests would be used to estimate the significance of enrichment for each pathway
##' 
##' Only the pathways with p value <= returned_pvalue in hypergegeometric tests and number of annotated genes >= returned_genenumber would be taken as enriched and returned.
##' 
##' @param gene a numeric matrix
##' @param returned_pvalue the minimum p value for enriched pathways
##' @param returned_genenumber the minimum number of annotated genes for enriched pathways
##' @inheritParams download_KEGGfile
##' @export
##' @return a list with two parts \item{name stastic}{description a matirx containing the pathway IDs of enriched pathways, and their names, p values, number of annotated genes}\item{name detail}{description a list with the genes annotated for each pathway}
##' @examples data(pho_sites_count)
##' #the 300 genes with most phospholation sites quantified
##' genes<-names(rev(sort(pho_sites_count[,1]))[1:300])
##' pho_KEGGresult<-find_enriched_pathway(genes,specis='hsa')
find_enriched_pathway<-function(gene,specis="hsa",returned_pvalue=0.05,returned_genenumber=5) {
	require(KEGG.db)
	keggpathway2gene <- as.list(KEGGPATHID2EXTID)
	pathway2name<-as.list(KEGGPATHID2NAME)
	keggpathway2gene<-keggpathway2gene[names(keggpathway2gene)[grep(specis,names(keggpathway2gene))]]
	##map NEED gene to kegg result
	kegg_result<-lapply(keggpathway2gene,function(x) x[x %in% gene])
	##ratio
	kegg_result_length<-unlist(lapply(kegg_result,length))
	keggpathway2gene_length<-unlist(lapply(keggpathway2gene,length))
	ratio<-kegg_result_length/keggpathway2gene_length
	##pvalue
	pvalue<-NULL
	for (x in 1:length(kegg_result_length)) {
		pvalue[x]<-phyper(kegg_result_length[x], keggpathway2gene_length[x], length(unique(unlist(keggpathway2gene)))-kegg_result_length[x], length(unique(unlist(kegg_result))), lower.tail = F)
	}
	pathway_ID<-substr(names(kegg_result_length),4,8)
	names(kegg_result)<-pathway_ID
	result<-data.frame(Pathway_Name=as.character(pathway2name[pathway_ID]),Gene_Found=kegg_result_length,Gene_Pathway=keggpathway2gene_length,Percentage=round(ratio,2),pvalue,stringsAsFactors=F)
	row.names(result)<-pathway_ID
	temp<-which(result[,5]<=returned_pvalue & result[,2]>=returned_genenumber)
	result<-result[temp,]
	kegg_result<-kegg_result[temp]
	return(list(stastic=result,detail=kegg_result))
}
