preprocess.gsets <- function(gsetsFile=NULL, gsets=NULL, bkg, minGset=NULL,
                             do.unique=FALSE, verbose=TRUE)
{
    if ( !xor(is.null(gsetsFile),is.null(gsets)) )
        stop( "must specify either gsetsFile or gsets" )
    if ( any(bkg!=toupper(bkg)) )
        stop( "bkg!=toupper(bkg)" )
    if ( !is.null(minGset) && (minGset<=0 || minGset>=length(bkg)) )
        stop( "minGset must be in (1,",length(bkg),")")

    if ( !is.null(gsetsFile) ) {
        gsets <- getGeneSet(new("GeneSet",gsetsFile))
    }
    gsets <- lapply(gsets,toupper)
    gsets <- lapply(lapply(gsets,intersect,bkg),sort)
    if ( !is.null(minGset) ) {
        keep.idx <- sapply(gsets,length)>=minGset
        if ( sum(keep.idx)<length(gsets) )
            VERBOSE(verbose,"removed",length(gsets)-sum(keep.idx),"genesets (too few genes)\n")
        gsets <- gsets[keep.idx]
    }
    if ( do.unique )
        gsets <- unique(gsets)

    return(gsets)
}
file.ext <- function( filen, w.sep=FALSE )
{
  ext <- rev(unlist(strsplit( filen, "\\.")))[1]
  if ( w.sep )
    ext <- paste(".",ext,sep="")
  ext
}
list2table <- function( obj, fill=NA )
{
  if ( !is.list(obj) )
    stop( "input object must be a list" )

  mx <- matrix( fill, nrow=max(sapply(obj,length)), ncol=length(obj), dimnames=list(NULL,names(obj)) )
  for ( i in 1:length(obj) ) {
    if ( length(obj[[i]])<1 ) next
    mx[1:length(obj[[i]]),i] <- obj[[i]]
  }
  mx
}
fast.mean <- function( x, cls=NULL, na.rm=FALSE )
{
  if (is.vector(x))
    return( fast.mean(matrix(x,1,length(x)),cls=cls) )
  if (is.null(dim(x)))
    stop( "a 2D object expected" )
  if ( any(is.na(x)) && !na.rm )
    stop( "NA values not admissible unless na.rm==T" )

  In <- rep(1,ncol(x))

  if ( is.null(cls) )
  {
    nc <- (!is.na(x)) %*% In
    x[is.na(x)] <- 0
    return ( (x %*% In)/ncol(x) )
  }
  # ELSE ..
  #
  lev <- sort(unique(cls))
  In <- sapply( 1:length(lev), function (z) as.numeric(cls==lev[z]) )

  nc <- (!is.na(x)) %*% In
  x[is.na(x)] <- 0
  x.mn <- (x %*% In)/nc
  colnames(x.mn) <- if (is.null(levels(cls))) lev else levels(cls)

  return(x.mn)

  if ( FALSE ) {
  x.sd <- fast.sd( x, cls=cls )

  In <- rep(1,ncol(x.mn))

  # weighted mean
  #
  mn <- ( ( x.mn * t(nc/t(x.sd)) ) %*% In ) / ( t(nc/t(x.sd)) %*% In )
  sd <- 1/( t(nc/t(x.sd^2)) %*% In )
  }
}
collapse2genesymbol <- function( eset, colId, method=c('max','median') )
{
    method <- match.arg(method)

    eset <- eset[fData(eset)[,colId]!="",]
    eset <- eset[!is.na(fData(eset)[,colId]),]

    if ( method=="max" ) {
        ord <- order( apply(exprs(eset), 1, mad), decreasing=TRUE )
        eset <- eset[ord,]
        eset <- eset[match(unique(fData(eset)[,colId]),fData(eset)[,colId]),]
        featureNames(eset) <- fData(eset)[,colId]
        return(eset)
    }
    else if ( method=="median" )
        stop( "median not implemented yet" )

}
########################################################################
## TISSUE PROJECTION
########################################################################
tissue.projection <- function( dat, barcode, threshold, max.tissues=1, minGset=5, do.unique=FALSE, stub="dat",
                               do.heatmap=FALSE, annot=NULL, annotCol=NULL, gsetFile=NULL, verbose=TRUE )
{
    ## BEGING checks on inputs
    ##
    if ( class(dat)!="ExpressionSet" )
        stop( "dat expected to be an ExpressionSet: ", class(dat) )
    if ( class(barcode)!="ExpressionSet" )
        stop( "barcode expected to be an ExpressionSet: ", class(barcode) )
    if ( threshold<0 | threshold>=1 )
        stop( "threshold must be in [0;1) range: ", threshold )
    if ( max.tissues<1 | max.tissues>ncol(dat) )
        stop( "max.tissues must be between 1 and ", ncol(dat),": ", max.tissues )
    ##
    ## END checks

    out.stub <-
        file.path(PATH,paste0("results/somascan/",stub,"GSVA.maxTissue",max.tissues,".thresh",threshold))

    ## thresholding: 'discretize' the [0;1] expression matrix into a {0,1} matrix
    barcode01 <- barcode
    exprs(barcode01)[exprs(barcode)<=threshold] <- 0
    exprs(barcode01)[exprs(barcode)> threshold] <- 1

    ## count the number of tissues each transcript is expressed in
    nt <- rowSums(exprs(barcode01))
    VERBOSE(verbose, "#( genes expressed in <=", max.tissues,"): ", sum(nt>0 & nt<=max.tissues), "\n")

    ## keep only the transcripts expressed in at most max.tissues
    barcodeFiltered <- barcode01[nt<=max.tissues,]

    ## extract tissue-specific genesets
    fData(dat)$symbol <- toupper(fData(dat)$symbol)
    tissueGSets <- apply(exprs(barcodeFiltered),2,
                      function(X) toupper(unique(fData(barcodeFiltered)$gene_symbol[X==1])))
    tissueGSets <- preprocess.gsets(gsets=tissueGSets,bkg=fData(dat)$symbol,minGset=minGset,do.unique=do.unique )

    ## saving genesets to a table and writing to file
    gsetTable <- list2table(tissueGSets,fill="")
    gsetFile <- paste0( out.stub, ".genesets.xlsx" )
    write.xlsx(gsetTable,file=gsetFile)
    cat( paste0("GeneSet table saved to '",gsetFile,"'\n") )

    ## geneset projection
    datGSVA <- gsva(dat, tissueGSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)$es.obs

    ## generating heatmap, if requested
    ##
    if ( do.heatmap )
    {
        hcRow <- hcopt(dist(exprs(datGSVA)),method="ward.D")
        hcCol <- hcopt(dist(t(exprs(datGSVA))),method="ward.D")
        pheatmap(exprs(datGSVA),main=paste0("genes in at most ",max.tissues," tissue(s) - threshold=",threshold),
                 color=colGradient(c("blue","white","red"),length=15),
                 annotation_col = annot,
                 annotation_colors = annotCol,
                 cluster_rows=hcRow,
                 cluster_cols=hcCol,
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 fontsize_row=5,
                 labels_row=sapply(featureNames(datGSVA),substr,1,20),
                 scale = "none")
        file.name <- paste0( out.stub, ".xlsx" )
        write.xlsx( exprs(datGSVA)[rev(hcRow$order),rev(hcCol$order)],rowNames=TRUE, file=file.name )
    }
    return(datGSVA)
}
########################################################################
## GSETS DIFFANAL
########################################################################

gsets.diffanal <- function(gspDat.file, do.prune=FALSE, max.fdr=0.01,min.genes=0,max.genes=+Inf,outpath="somascan/",
                           do.heatmap=TRUE, fontsize_row=5, annot_col=NULL, annotColors=NULL,group="Centenarian",
                           heat.col=colGradient(c("purple","white","darkorange"),length=11),
                           border_color=NA,interactive=FALSE,do.write=TRUE,verbose=TRUE)
{
    ## load the data
    gspDat <- readRDS(file=gspDat.file)

    ## shorten pathway labels
    if ( do.prune ) {
        featureNames(gspDat) <- gsub("HALLMARK_","",featureNames(gspDat))
        featureNames(gspDat) <- gsub("REACTOME_","R_",featureNames(gspDat))
        featureNames(gspDat) <- gsub("KEGG_","K_",featureNames(gspDat))
        featureNames(gspDat) <- gsub("BIOCARTA_","B_",featureNames(gspDat))
    }
    ## perform differential analysis based on LM
    coeffIdx <- paste0("Cohort",group)
    gspDiff <- as.data.frame(
        t(sapply( 1:nrow(gspDat),
           function(i) {
               tmp <- lm(exprs(gspDat)[i,] ~ Cohort + Gender + Year.completed, data=pData(gspDat))
               summary(tmp)$coefficients[coeffIdx,c("Estimate","Pr(>|t|)")]
           }))
    )
    rownames(gspDiff) <- featureNames(gspDat)
    gspDiff <- gspDiff[order(gspDiff[,"Pr(>|t|)"]),]
    gspDiff$fdr <- p.adjust(gspDiff[,"Pr(>|t|)"], method="BH")
    gspDiff$group <- c("UP","DN")[as.numeric(gspDiff$Estimate<0)+1]

    if ( do.write )
    {
        ## save output in binary ..
        outfile <- paste0("workdir/",outpath,gsub(".RDS","",basename(gspDat.file)),".gspDiff.RDS")
        saveRDS(gspDiff,file=file.path(PATH,outfile))
        VERBOSE(verbose, paste0("diffanal binary file saved to '", outfile, "'\n") )

        ## .. and 'xlsx' format
        outfile <- paste0("results/",outpath,gsub(".RDS","",basename(gspDat.file)),".gspDiff.xlsx")
        write.xlsx(gspDiff,file=file.path(PATH,outfile),rowNames=TRUE)
        VERBOSE(verbose, paste0("diffanal excel file saved to '", outfile, "'\n") )
    }
    if (!do.heatmap) return(gspDiff)

    ## get significant markers first
    up <- gspDiff[gspDiff$Estimate>0 & gspDiff$fdr<=max.fdr ,]
    up <- up[order(up[,"Estimate"],decreasing=TRUE), ][1:min(max.genes,nrow(up)),,drop=FALSE]
    dn <- gspDiff[gspDiff$Estimate<0 & gspDiff$fdr<=max.fdr,]
    dn <- dn[order(dn[,"Estimate"],decreasing=FALSE),][1:min(max.genes,nrow(dn)),,drop=FALSE]
    tmp <- rbind(up,dn); tmp <- tmp[order(tmp[,"Estimate"],decreasing=TRUE),]

    gspDat1 <- gspDat[rownames(tmp),order(gspDat$Cohort),drop=FALSE]
    VERBOSE(verbose, nrow(gspDat1), "significant genesets at FDR <=",max.fdr,"\n")

    ## add row annotation
    annot_row <-
        data.frame(markers=c("down","up")[as.integer(gspDiff[featureNames(gspDat1),"Estimate"]>0)+1])
    rownames(annot_row) <- featureNames(gspDat1)
    annotColors$markers <- c(down="darkgreen",up="darkred")

    ## heatmap sorted by phenotype (columns) and markers esimates (columns)
    pheatmap.group(gspDat1, main="clustered within phenotype, sorted by markers",
         group="Cohort",
         color=heat.col,
         annotation_col = annot_col[sampleNames(gspDat1),,drop=FALSE],
         annotation_row = annot_row,
         annotation_colors = annotColors,
         cluster_rows=FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")
    if ( interactive ) { cat("next> "); scan() }

    ## heatmap sorted by phenotype (columns) and markers esimates (columns)
    pheatmap(exprs(gspDat1), main="sorted by phenotype and markers",
         color=heat.col,
         annotation_col = annot_col,
         annotation_row = annot_row,
         annotation_colors = annotColors,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")
    if ( interactive ) { cat("next> "); scan() }

    ## heatmap sorted by clustering
    hcRow <- hcopt(dist(exprs(gspDat1)),method="ward.D")
    hcCol <- hcopt(dist(t(exprs(gspDat1))),method="ward.D")
    pheatmap(exprs(gspDat1), main="sorted by clustering",
         color=heat.col,
         annotation_col = annot_col,
         annotation_row = annot_row,
         annotation_colors = annotColors,
         cluster_rows=hcRow,
         cluster_cols=hcCol,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")

    return(gspDiff)
}
gsets.longevity.diffanal <- function(gspDat, do.prune=FALSE, max.fdr=NULL,max.pval=NULL,outpath="somascan/",base.file,
                                     do.heatmap=TRUE, fontsize_row=12, annot_col=NULL, annotColors=NULL,
                                     heat.col=colGradient(c("purple","white","darkorange"),length=11),
                                     border_color=NA,interactive=FALSE,do.write=TRUE,verbose=TRUE)
{
    if (!xor(is.null(max.fdr),is.null(max.pval)))
        stop("can only select max.fdr or max.pval, not both")

    ## shorten pathway labels
    if ( do.prune ) {
        featureNames(gspDat) <- gsub("HALLMARK_","",featureNames(gspDat))
        featureNames(gspDat) <- gsub("REACTOME_","R_",featureNames(gspDat))
        featureNames(gspDat) <- gsub("KEGG_","K_",featureNames(gspDat))
        featureNames(gspDat) <- gsub("BIOCARTA_","B_",featureNames(gspDat))
    }
    ## perform differential analysis based on LM
    idx <- paste0("gspDat$FU.2group",levels(gspDat$FU.2group)[2])
    gspDiff <- as.data.frame(
        t(sapply( 1:nrow(gspDat),
           function(i) {
               fit <- lm(exprs(gspDat)[i,] ~ gspDat$Gender+gspDat$age.group+gspDat$FU.2group )
               summary(fit)$coefficients[idx,c("Estimate","t value","Pr(>|t|)")]
           }))
    )
    rownames(gspDiff) <- featureNames(gspDat)
    gspDiff <- gspDiff[order(gspDiff[,"Pr(>|t|)"]),]
    gspDiff$fdr <- p.adjust(gspDiff[,"Pr(>|t|)"], method="BH")
    gspDiff$group <- levels(gspDat$FU.2group)[as.numeric(gspDiff$Estimate>0)+1]

    if ( do.write )
    {
        ## save output in binary ..
        outfile <- paste0("workdir/",outpath,base.file,".longevity.gspDiffanal.RDS")
        saveRDS(gspDiff,file=file.path(PATH,outfile))
        VERBOSE(verbose, paste0("diffanal binary file saved to '", outfile, "'\n") )

        ## .. and 'xlsx' format
        outfile <- paste0("results/",outpath,base.file,".longevity.gspDiffanal.xlsx")
        write.xlsx(gspDiff,file=file.path(PATH,outfile),rowNames=TRUE)
        VERBOSE(verbose, paste0("diffanal excel file saved to '", outfile, "'\n") )
    }
    if (!do.heatmap) return(gspDiff)

    ## get significant markers first
    tmp <- {
        if (is.null(max.pval))
            gspDiff[gspDiff$fdr<=max.fdr,]
        else
            gspDiff[gspDiff[,"Pr(>|t|)"]<=max.pval,]
    }
    tmp <- tmp[order(tmp[,"Estimate"],decreasing=TRUE),]
    gspDat1 <- gspDat[rownames(tmp),order(gspDat$FU.2group)]
    VERBOSE(verbose, nrow(gspDat1), "significant genesets at q-/p-value <=",
            if(is.null(max.pval)) max.fdr else max.pval,"\n")

    ## add row annotation
    annot_row <-
        data.frame(markers=c("down","up")[as.integer(gspDiff[featureNames(gspDat1),"Estimate"]>0)+1])
    rownames(annot_row) <- featureNames(gspDat1)
    annotColors$markers <- c(down="darkgreen",up="darkred")

    ## heatmap clustered within phenotype (rows) and sorted by markers esimates (columns)
    pheatmap.group(gspDat1, main="clustered within phenotype, sorted by markers",
         group="FU.2group",
         color=heat.col,
         annotation_col = annot_col[sampleNames(gspDat1),,drop=FALSE],
         annotation_colors = annotColors,
         annotation_row = annot_row,
         cluster_rows=FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize=fontsize_row,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")

    if ( interactive ) { cat("next> "); scan() }

    ## heatmap sorted by phenotype (rows) and markers esimates (columns)
    pheatmap(exprs(gspDat1), main="sorted by phenotype and markers",
         color=heat.col,
         annotation_col = annot_col,
         annotation_colors = annotColors,
         annotation_row = annot_row,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize=fontsize_row,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")

    if ( interactive ) { cat("next> "); scan() }

    ## heatmap sorted by clustering
    hcRow <- hcopt(dist(exprs(gspDat1)),method="ward.D")
    hcCol <- hcopt(dist(t(exprs(gspDat1))),method="ward.D")
    pheatmap(exprs(gspDat1), main="sorted by clustering",
         color=heat.col,
         annotation_col = annot_col,
         annotation_colors = annotColors,
         annotation_row = annot_row,
         cluster_rows=hcRow,
         cluster_cols=hcCol,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize=fontsize_row,
         fontsize_row=fontsize_row,
         labels_row=sapply(featureNames(gspDat1),substr,1,27),
         border_color=border_color,
         scale = "none")

    return(gspDiff)
}
###################################################
## functions used in analysis08.apoe.projection.Rmd
###################################################

## extract "01" (primary), "02" (metastatic), or "11" (adjacent normal) codes from sample names
tcga.sample.code <- function( snames, tail=TRUE )
{
    if ( tail )
        sapply(snames,function(Z) rev(unlist(strsplit(Z,"-")))[1])
    else
        sapply(snames,function(Z) paste(rev(rev(unlist(strsplit(Z,"-")))[-1]),collapse="-"))
}
## wrapper function
run.projection <- function( gep, sig, out.dir, title, seed=123, do.heatmap=TRUE, verbose=TRUE, save.plot=TRUE )
{
    if ( length(cmn <- intersect(featureNames(gep),sig[[1]]))<1 ) stop( "empty signature" )
    gep <- gep[cmn,]

    VERBOSE(verbose,"Running ASSIGN ...\n")

    ## RUN ASSIGN
    dir.create(out.dir,recursive=TRUE)
    set.seed(seed)
    assign.wrapper(geneList=sig,
                   testData=exprs(gep),
                   trainingLabel=NULL,
                   outputDir=out.dir,
                   n_sigGene=NA,
                   adaptive_B=TRUE,
                   adaptive_S=TRUE,
                   mixture_beta=TRUE,
                   iter=3000)

    OUT <- load.var(file=paste(out.dir,'output.rda',sep='/'))
    OUT <- annotateAssignOutput(assignOutput=OUT)
    saveRDS(OUT,file=file.path(out.dir,'annotatedOutput.RDS'))
    VERBOSE(verbose,"done.\n")

    if (do.heatmap)
    {
        VERBOSE(verbose,"Generating heatmap ... ")

        if (save.plot) png( file.path(out.dir,'assign.heatmap.%02d.png'), width=900, height=900)
        assign.heatmap.wrapper(EXP=exprs(gep),
                               assignOutput=OUT,
                               colSideCols=NULL,
                               do.barplot=TRUE)
        #plotAll(eSet=gep, output.data=OUT, gene_list=sig[[1]], do.log=FALSE, simple.heatmap=FALSE,
        #        title=title,ylab=paste("APOE signature",names(sig)), do.barplot=TRUE)
        if (save.plot) dev.off()
        VERBOSE(verbose,"done.\n")
    }
    OUT
}
#################################################
## functions used in analysis09.apoe.survival.Rmd
#################################################

## preprocess survival data: extract from the tcga file and create a
## .. Surv data object
process.survival <- function( filen )
{
    Tab <- read.delim(filen,sep="\t",comment.char="",skip=1,header=TRUE,row.names=1)
    TIME <- Tab$last_followup_days
    STATUS <- Tab$Status
    print(table(STATUS),useNA="ifany")
    if ( grepl("tumor",STATUS[1]) ) {
        STATUS <- as.integer(STATUS=="with tumor")
        print(table(STATUS))
    }
    else {
        STATUS <- as.integer(STATUS=="dead")
        print(table(STATUS))
    }
    Srv <- Surv(TIME,STATUS)
    rownames(Srv) <- rownames(Tab)
    Srv
}
## preprocess ASSIGN projection output: extract 'activity' scores
process.tcga.projection <- function( filen )
{
    proj <- readRDS(filen)
    proj <- proj$mcmc.pos.mean.testData$kappa_pos
    rownames(proj) <- sapply(rownames(proj),tcga.sample.code,tail=FALSE)
    proj
}
## test for the stratification of survival as a function of projection score
test.tcga.survival <- function( Prj=NULL, Prj.file=NULL, surv, main=NULL,
                                threshold=c("median","mean"), col=c("green","red") )
{
    if ( !xor(is.null(Prj),is.null(Prj.file)) ) stop( "only Prj or Prj.file can be non-NULL" )
    threshold <- match.arg(threshold)
    if ( !is.null(Prj.file) ) {
        Prj <- process.tcga.projection(Prj)
    }
    Srv <- process.survival(surv)

    cmn <- intersect(rownames(Srv),rownames(Prj))
    if ( length(cmn)<1 ) stop( "no matching between survival and projection rownames" )
    Prj1 <- Prj[match.nona(cmn,rownames(Prj)),]
    Srv1 <- Srv[match.nona(cmn,rownames(Srv))]
    if ( any(rownames(Prj1)!=rownames(Srv1)) ) stop( "rownames(Prj1)!=rownames(Srv1)" )

    thresh.fun <- match.fun(threshold)
    I01 <- as.numeric(Prj1>=thresh.fun(Prj1))
    plot(survfit(Srv1 ~ I01),col=col, mark.time=FALSE, main=main)
    pvl <- survdiff(Srv1 ~ I01)
    legend("bottomleft",col=col,lty=1,lwd=2,
           legend=paste0("score",c(">=","<"),threshold, " (n=", table(I01),")"))
    legend("topright",legend=paste0("pvalue=",format.pval(survdiff.pvalue(pvl))))

    survdiff.pvalue(pvl)
}
## generate the TCGA survival file path
tcga.surv.path <- function(type, path=file.path(PATH,"data/tcga/survival"),
                           sty=c("survival","diseaseFree"))
{
    sty <- match.arg(sty)
    stub <- file.path(path,type,gsub("_","",tcga.datasets[type,"date"]),type)
    if ( sty=="survival" )
        paste0( stub, ".clin.merged.txt_survival_table.txt" )
    else if ( sty=="diseaseFree" )
        paste0( stub, ".clin.merged.txt_DiseaseFree_survival_table.txt" )
    else
        stop( "unrecognized sty:", sty )
}
## SURVDIFF PVALUE
survdiff.pvalue <- function( sobj )
{
  if ( !is.list(sobj) | names(sobj)[5]!="chisq" ) stop( "Survdiff object expected" )
  return( 1-pchisq(sobj$chisq,length(sobj$n)-1) )
}
## PHEATMAP GROUP

pheatmap.group <- function(
    mat,
    group,
    color,
    annotation_col,
    annotation_row=NULL,
    annotation_colors,
    cluster_rows=TRUE,
    border_color=NA,
    ...
)
{
    ## create a new 'meta-group'
    if ( length(group)>1 ) {
        if ( length(group)>2 ) stop( "length(group)>2" )
        newGroup <- factor(apply(pData(dat)[,group],1,paste,collapse="."),
                           levels=as.vector(sapply( levels(pData(dat)[,group[1]]),
                               function(z) paste(z,levels(pData(dat)[,group[2]]),sep="."))))
        dat$group <- newGroup
        group <- "group"

    }
    ## BEGIN input checks
    if ( !is.null(annotation_row) && any(rownames(annotation_row)!=featureNames(mat)) )
        stop( "row mismatch")
    if ( nrow(annotation_col)!=ncol(mat) )
        stop( "nrow(annotation_col)!=ncol(mat)" )
    if ( any(rownames(annotation_col)!=sampleNames(mat)) )
        stop( "rownames(annotation_col)!=sampleNames(mat)" )
    if ( !(group %in% colnames(pData(mat))) )
        stop( "group must be one of the pData's columns" )
    if ( !is.factor(pData(mat)[,group]) )
        stop( "!is.factor(pData(mat)[,group])" )
    if ( is.null(levels(pData(mat)[,group])) )
        stop( "is.null(levels(pData(mat)[,group]))" )
    ## END input checks

    hc.row <- {
        if (cluster_rows)
            hcopt(as.dist(1-cor(t(exprs(mat)),use="pairwise.complete.obs")),method="ward.D")
        else
            FALSE
    }
    ORD <- NULL

    for ( grp in levels(pData(mat)[,group]) )
    {
        matI <- mat[,pData(mat)[,group]==grp]
        hc.col <- hcopt(dist(t(exprs(matI))),method="ward.D")
        ORD <- c(ORD,sampleNames(matI)[hc.col$order])
    }
    mat <- mat[,ORD]
    pheatmap(exprs(mat),
             color=color,
             annotation_col=annotation_col,
             annotation_row=annotation_row,
             annotation_colors=annotation_colors,
             cluster_cols=FALSE,
             cluster_rows=hc.row,
             border_color=border_color,
             ...)

}
## MODIFIED ksGenescore

ksGenescore <- function
(
 n.x,               # length of ranked list
 y,                 # positions of geneset items in ranked list (basically, ranks)
 do.pval=T,         # compute asymptotic p-value
 alternative=c("two.sided","greater","less"),
 do.plot=FALSE,     # draw the ES plot
 bare=FALSE,        # return score & p-value only (a 2-tuple)
 weight=NULL,       # weights for weighted score (see Subramanian et al.) (usually, sort(score))
 weight.p=1,        # weights' exponent
 cls.lev=c(0,1),    # class labels to display
 absolute=FALSE,    # takes max - min score rather than the maximum deviation from null
 plot.labels=FALSE, # hits' labels
 ...                # additional plot arguments
)
{
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative)
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  x.axis <- y.axis <- NULL

  ## weighted GSEA score
  ##
  if ( !is.null(weight) )
  {
    weight <- abs(weight[y])^weight.p

    Pmis <- rep(1, n.x); Pmis[y] <- 0; Pmis <- cumsum(Pmis); Pmis <- Pmis/(n.x-n.y)
    Phit <- rep(0, n.x); Phit[y] <- weight; Phit <- cumsum(Phit); Phit <- Phit/Phit[n.x]
    z <- Phit-Pmis

    score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]

    x.axis <- 1:n.x
    y.axis <- z
  }
  ## KS score
  ##
  else
  {
      y <- sort(y)
      n <- n.x * n.y/(n.x + n.y)
      hit <- 1/n.y
      mis <- 1/n.x

      ## to compute score, only the y positions and their immediate preceding
      ## ..positions are needed
      ##
      Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y)
      D <- rep( 0, length(Y) ); D[y.match] <- (1:n.y)
      zero <- which(D==0)[-1]; D[zero] <- D[zero-1]

      z <- D*hit - Y*mis

      score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]

      x.axis <- Y;
      y.axis <- z;
      if(Y[1]>0) {
          x.axis <- c(0,x.axis);
          y.axis <- c(0,y.axis);
      }
      if ( max(Y)<n.x ) {
          x.axis <- c(x.axis,n.x)
          y.axis <- c(y.axis,0)
      }
  }
  if (do.plot) {
      plot(x.axis, y.axis, type="l",
           #xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
           #    "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
           xlab="",ylab="",...)
      abline(h=0)
      abline(v=n.x/2,lty=3)
      axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
      i.max <- which.max(abs(y.axis))
      #points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
      points( x.axis[i.max], y.axis[i.max], pch=20,cex=2.5, col="red")
      text(x.axis[i.max]+n.x/20,y.axis[i.max]-0.025,round(y.axis[i.max],2),cex=2)
  }
  if ( !do.pval ) {
      return( as.numeric(score) )
  }
  ## ELSE compute p-value as in function ks.test but return signed statistic
  ##
  tmp <- suppressWarnings(ks.test(1:n.x,y,alternative=alternative))
  tmp$statistic <- score # use the signed statistic
  return( if (bare) c(score=as.numeric(tmp$statistic), p.value=tmp$p.value) else tmp )
}
colGradient <- function( cols, length, cmax=255 )
{
  ## e.g., to create a white-to-red gradient with 10 levels
  ##
  ##   colGradient(cols=c('white','red'),length=10)
  ##
  ## or, to create a blue-to-white-to-red gradients with 9 colors (4 blue's, white, 4 red's)
  ##
  ##   colGradient(cols=c('blue','white','red'),length=9)
  ##
  ramp <- colorRamp(cols)
  rgb( ramp(seq(0,1,length=length)), maxColorValue=cmax )
}
## replicate a matrix or data.frame row- or column-wise
repmat <- function(MX,MARGIN=2,n,each=FALSE)
{
  m <- if ( MARGIN==1 ) ncol(MX) else nrow(MX)
  idx <- if ( each ) rep(1:m,each=n) else rep(1:m,times=n)

  if ( MARGIN==1 )
    return(MX[,idx])
  else if ( MARGIN==2 )
    return(MX[idx,])
  else
    stop("MARGIN:",MARGIN)
}
