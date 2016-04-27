####################################################
#
#		samplingProbs <- function(...)
#
#		tree = object of class phylo
#		cladeTable = a dataframe with 1 column of species names and a second column of group assignment
#			Must either be a table for all species in tree, or a table of greater species richness, including those species in the tree.
#		cladeRichness = either NULL or a vector of species counts, named by group names.
#		globalSampling = percent sampling of the backbone of the phylogeny
#		output = path + output file name (.txt)
#		writeToDisk = boolean, should the table be written to disk, defaults to TRUE

##' @title Creates clade-specific sampling fractions
##'
##' @description If the user would like to specify species sampling on a
##'     clade-by-clade basis, a sampling probability table can be provided to
##'     \code{BAMM}.
##'
##' @param tree An object of class \code{phylo}.
##' @param cladeTable A dataframe with one column of species names and a
##'     second column of group assignment. 
##' @param cladeRichness Either \code{NULL} or a vector of species counts,
##'     named by clade names.
##' @param globalSampling percent sampling for the backbone of the phylogeny.
##' @param output Path + output file name.
##' @param writeToDisk A logical, should the table be written to disk,
##'     defaults to \code{TRUE}.
##'
##' @details This function handles two types of input: The cladeTable can
##'     either contain the species found in the phylogeny, along with clade
##'     assignment of those species, or it can contain more species than found
##'     in the phylogeny. If the table only contains those species in the
##'     phylogeny, then a vector \code{cladeRichness} must be provided that
##'     contains known clade richness. If the cladeTable contains more than
##'     the species in the phylogeny, then cladeRichness should be set to
##'     \code{NULL}. The \code{globalSampling} value should represent the
##'     overall completeness of species sampling in terms of major clades. See
##'     \url{http://bamm-project.org} for additional details.
##'
##' @return If \code{writeToDisk = TRUE}, then no object is returned. If
##'     \code{writeToDisk = FALSE}, then a dataframe is returned. The
##'     resultant table must contain one row for each species in the
##'     phylogeny, along with clade assignment, and sampling fraction. The
##'     first line must contain the overall sampling fraction for the
##'     phylogeny and must be written as tab-delimited, with no headers.
##'
##' @author Pascal Title
##'
##' @examples
##' # Generate dummy data
##' tree <- read.tree(text="(((t1:2,(t2:1,t3:1):1):1,((t4:1,t5:1):1,t6:2):1)
##'                   :1,(t7:3,(t8:2,t9:2):1):1);")
##' tree$tip.label <- paste(rep('Species',9),1:9,sep='')
##'
##' spTable <- as.data.frame(matrix(nrow=9,ncol=2))
##' spTable[,1] <- tree$tip.label
##' spTable[1:3,2] <- 'cladeA'
##' spTable[4:6,2] <- 'cladeB'
##' spTable[7:9,2] <- 'cladeC'
##' richnessVec <- c(cladeA=5, cladeB=4, cladeC=12)
##'
##' # Option 1: We have a table of clade assignment for the species in the
##' #           tree, along with a vector of known clade richness
##' spTable
##' richnessVec
##' samplingProbs(tree, cladeTable = spTable, cladeRichness = richnessVec,
##'               globalSampling = 1, writeToDisk = FALSE)
##' 
##' # Option 2: We have a table of known species, beyond the sampling in the
##' #           phylogeny
##' spTable <- rbind(spTable, c('Species10','cladeA'),c('Species11','cladeA'),
##'                  c('Species12','cladeC'), c('Species13','cladeC'),
##'                  c('Species14','cladeC'))
##'
##' spTable
##'
##' samplingProbs(tree, cladeTable = spTable, cladeRichness = NULL, 
##'               globalSampling = 0.9, writeToDisk = FALSE)
##' @export
samplingProbs <- function(tree, cladeTable, cladeRichness = NULL, globalSampling, output, writeToDisk = TRUE) {
	
	if (length(intersect(tree$tip.label,cladeTable[,1])) != length(tree$tip.label)) {
		stop("Not all species from tree are in cladeTable.");
	}
	
	if (nrow(cladeTable) == length(tree$tip.label)) {
		if (is.null(cladeRichness)) {
			stop("If cladeTable only contains species from tree, then cladeRichness must be provided.");
		}
	}
	
	if (nrow(cladeTable) > length(tree$tip.label)) {
		if (!is.null(cladeRichness)) {
			warning("cladeTable contains more species than in tree, so cladeRichness vector will be ignored.");
		}
	}
	
	if (!is.null(cladeRichness)) {
		if (length(cladeRichness) != length(unique(cladeTable[,2]))) {
			stop("The cladeRichness vector must contain the same number of clades as are described in the cladeTable.");
		}
	}
	
	if (!is.vector(cladeRichness) & !is.null(cladeRichness)) {
		stop("Error: cladeRichness must either be NULL or an integer vector named with clade names.");
	}
	
	if (ncol(cladeTable) > 2) {
		stop("cladeTable must contain 2 columns: one of species, and one of clade assignment.");
	}
	
	if (is.matrix(cladeTable)) {
		cladeTable <- as.data.frame(cladeTable, stringsAsFactors=FALSE);
	}
	
	if (nrow(cladeTable) > length(tree$tip.label)) {
		probs <- as.data.frame(matrix(nrow=length(tree$tip.label),ncol=3));
		colnames(probs) <- c('sp','clade','prob');
		for (i in 1:length(tree$tip.label)) {
			probs[i,1] <- tree$tip.label[i];
			clade <- cladeTable[cladeTable[,1] == tree$tip.label[i],2];
			inTree <- intersect(cladeTable[cladeTable[,2] == clade,1],tree$tip.label);
			probs[i,2] <- clade;
			probs[i,3] <- length(inTree) / length(cladeTable[cladeTable[,2] == clade,1]);
		}
		probs <- rbind(c(globalSampling,'','',''),probs);
	}
	
	if (nrow(cladeTable) == length(tree$tip.label) & !is.null(cladeRichness)) {
		probs <- as.data.frame(matrix(nrow=length(tree$tip.label),ncol=3));
		colnames(probs) <- c('sp','clade','prob');
		for (i in 1:length(tree$tip.label)) {
			probs[i,1] <- tree$tip.label[i];
			clade <- cladeTable[cladeTable[,1] == tree$tip.label[i],2];
			probs[i,2] <- clade;
			probs[i,3] <- nrow(cladeTable[cladeTable[,2] == clade,]) / cladeRichness[clade];
		}
		probs <- rbind(c(globalSampling,'','',''),probs);
	}
	if (writeToDisk) {
		write.table(probs, file=output, quote=F, col.names=F, row.names=F, sep='\t');
	} else {
		return(probs);
	}
}








