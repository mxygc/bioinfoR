#' Scripts of the factor analysis based pipeline
#' for the integrative omics
if (!exists("FactAnal") || !is.environment(FactAnal)) FactAnal <- new.env(parent = emptyenv())
local ({
    normalize.exprs <- function(exprs, normalization) {
        dat <- exprs
        method <- normalization
        fa.r.norm.minmax <- function(dat.row) {
        	min.row <- min(dat.row)
        	max.row <- max(dat.row)
        	return ((dat.row - min.row) / (max.row - min.row))
        }
        
        fa.r.norm.mean <- function(dat.row) {
        	mean.row <- mean(dat.row)
        	return ((dat.row - mean.row) / mean.row)
        }
        
        fa.r.norm.z <- function(dat.row) {
        	mean.row <- mean(dat.row)
        	sd.row <- sd(dat.row)
        	return ((dat.row - mean.row) / sd.row)
        }
    	if (method == "minmax") {
    		fa.dat <- apply(dat, 1, FUN = fa.r.norm.minmax)
    	}
    	else if (method == "z") {
    		fa.dat <- apply(dat, 1, FUN = fa.r.norm.z)
    	}
    	else if (method == "mean") {
    		fa.dat <- apply(dat, 1, FUN = fa.r.norm.mean)
    	}
    	return (t(fa.dat))
    }

        scree.plot <- function(dat, rep, cent) {
    	ev <- eigen(cor(dat)) # get eigenvalues
    	ap <- parallel(subject = nrow(dat), var = ncol(dat), rep = 100, cent = .05)
    	nS <- nScree(ev$values, ap$eigen$qevpea)
    	plotnScree(nS)
    }

    get.fa.comb <- function(fa.number) {
    	fa.combs <- c()
    	if (fa.number == 1) {
    		fa.combs <- c(fa.combs, "F1")
    		return (fa.combs)
    	}
    	else {
    		for (i in 1:(fa.number - 1)) {
    			fa.comb <- c(paste("F", i, sep = ""))
    			for (j in (i + 1):fa.number ) {
    				fa.combs <- c(fa.combs, fa.comb)
    				fa.comb <- paste(fa.comb, "+", "F", j, sep = "")
    			}
    			fa.combs <- c(fa.combs, fa.comb)
    		}
    		fa.combs <- c(fa.combs, paste("F", fa.number, sep = ""))
    		return (fa.combs)
    	}
    }
    
    get.lda.predict.group <- function(fa.combination, fa.loading, group, test.method, cv = TRUE) {
    	lda.data <- data.frame(fa.loading, group)
    	lda.formula <- as.formula(paste("group", "~", fa.combination))
        if (cv) {
    	    lda.predict <- lda(lda.formula, data = lda.data, CV = T)$class
        } else {
            lda.predict <- predict(lda(lda.formula, data = lda.data, CV = F))$class
        }
    	lda.cmp <- table(group, lda.predict)
    	lda.predict.accuracy <- sum(diag(lda.cmp)) / length(group)
    	if (sum(grep("fisher*", test.method))) {
    		lda.predict.p <- fisher.test(lda.cmp)$p.value
    	}
    	if (sum(grep("chi*", test.method))) {
    		lda.predict.p <- chisq.test(lda.cmp)$p.value
    	}
    	if (sum(grep("mcnemar*", test.method))) {
    		lda.predict.p <- mcnemar.test(lda.cmp, correct = FALSE)$p.value
    	}
    	return (list(accuracy = lda.predict.accuracy, p = lda.predict.p))
    }
    
    get.lda.predict <- function(fa.loading, fa.combinations, fa.GRP, test.method = "fisher", cv = TRUE) {
    	len.col <- length(fa.GRP)
    	len.row <- length(fa.combinations)
    	lda.predict.accuracy.table <- matrix(NA, len.row, len.col)
    	lda.predict.p.table <- matrix(NA, len.row, len.col)
    	for (i in 1:len.col) {
    		lda.predict <- lapply(fa.combinations, FUN = get.lda.predict.group, fa.loading, fa.GRP[[i]], test.method, cv)
    		lda.predict.accuracy <- lapply(lda.predict, function(x) {return (x$accuracy)})
    		lda.predict.p <- lapply(lda.predict, function(x) {return (x$p)})
    		lda.predict.accuracy.table[, i] <- unlist(lda.predict.accuracy)
    		lda.predict.p.table[, i] <- unlist(lda.predict.p)
    	}
    	group.names <- names(fa.GRP)
    	rownames(lda.predict.accuracy.table) <- fa.combinations
    	colnames(lda.predict.accuracy.table) <- group.names
    	rownames(lda.predict.p.table) <- fa.combinations
    	colnames(lda.predict.p.table) <- group.names
    
    	return (list(accuracy = lda.predict.accuracy.table, p = lda.predict.p.table))
    }
    
    get.lda.best <- function(lda.predict, fa.number, fa.loading, p.cutoff) {
    	lda.predict.accuracy <- lda.predict$accuracy
    	lda.predict.p <- lda.predict$p
    	fa.groups <- colnames(lda.predict.accuracy)
    	names(fa.groups) <- fa.groups
    	fa.combs <- rownames(lda.predict.accuracy)
    	fa.comb.lengths <- unlist(lapply(fa.combs, function(fa.comb) {length(unlist(strsplit(fa.comb, "+", fixed = TRUE)))}))
    	get.lda.best.index <- function(fa.group, p.cutoff) {
    		lda.predict.accuracy.group <- lda.predict.accuracy[, fa.group]
    		lda.predict.p.group <- lda.predict.p[, fa.group]
    		small.p.index <- (1:length(fa.combs))[lda.predict.p.group <= p.cutoff]
    		if(length(small.p.index) > 0) {
    			accuracy.with.small.p <- lda.predict.accuracy.group[small.p.index]
    			max.accuracy.index <- small.p.index[which(max(accuracy.with.small.p) == accuracy.with.small.p)]
    			if(length(max.accuracy.index) > 1) {
    				comb.length <- fa.comb.lengths[max.accuracy.index]
    				min.comb.length.index <- max.accuracy.index[which.min(comb.length)]
    				return (min.comb.length.index)
    			}
    			else {
    				return (max.accuracy.index)
    			}
    		}
    		else {
    			return (NA)
    		}
    	}
    	get.lda.best.accuracy <- function(fa.group, lda.best.index) {
    		if(is.na(lda.best.index[[fa.group]])) {
    			return (NA)
    		}
    		else {
    			lda.best.accuracy <- lda.predict.accuracy[lda.best.index[[fa.group]], fa.group]
    			return (lda.best.accuracy)
    		}
    	}
    
    	get.lda.best.p <- function(fa.group, lda.best.index) {
    		if(is.na(lda.best.index[[fa.group]])) {
    			return (NA)
    		}
    		else {
    			lda.best.p <- lda.predict.p[lda.best.index[[fa.group]], fa.group]
    			return (lda.best.p)
    		}
    	}
    	get.lda.best.model <- function(fa.group, lda.best.index) {
    		if(is.na(lda.best.index[[fa.group]])) {
    			return ("")
    		}
    		else {
    			lda.best.model <- fa.combs[lda.best.index[[fa.group]]]
    		}
    	}
    
    	lda.best.index <- lapply(fa.groups, FUN = get.lda.best.index, p.cutoff)
    	lda.best.accuracy <- lapply(fa.groups, FUN = get.lda.best.accuracy, lda.best.index)
    	lda.best.p <- lapply(fa.groups, FUN = get.lda.best.p, lda.best.index)
    	lda.best.model <- lapply(fa.groups, FUN = get.lda.best.model, lda.best.index)
    
    	result <- list(accuracy = unlist(lda.best.accuracy), p = unlist(lda.best.p), model = unlist(lda.best.model))
    	return (result)
    
    }
    
    get.lda.summary <- function(lda.best.accuracy, lda.best.p, lda.best.model) {
    	n.col <- ncol(lda.best.accuracy)
    	n.row <- nrow(lda.best.accuracy)
    	lda.summary <- matrix("", nrow = n.row, ncol = n.col)
    	rownames(lda.summary) <- rownames(lda.best.accuracy)
    	colnames(lda.summary) <- colnames(lda.best.accuracy)
    	for (i in 1:n.row) {
    		for (j in 1:n.col) {
    			lda.summary[i,j] <- paste(lda.best.model[i,j], "(", round(lda.best.accuracy[i,j], digits = 3), ",", round(lda.best.p[i,j], digits = 3), ")", sep = "")
    		}
    	}
    	return (lda.summary)
    }
    
    do.factor.analysis <- function(exprs, normalization = "mean", p.cutoff = 0.05, fa.numbers = NULL, fa.groups, fa.factors, test.method = "fisher.test", cv = TRUE, dir.out, file.loading = "loading.csv", file.score = "score.csv", file.p = "p.csv", file.accuracy = "accuracy.csv", file.summary = "summary.csv", file.scree = "screeplot.pdf", plot = TRUE) {
        nsample <- ncol(exprs)
        ngroup <- nlevels(as.factor(fa.groups))
        fa.dat <- normalize.exprs(exprs, normalization)
        fscree <- file.path(dir.out, file.scree)
        dir.create(dir.out, F, T)
        if (plot) { 
        	pdf(file = file.scree)
        	tryCatch(scree.plot(fa.dat), error = function(e) {
                     dev.off()
                     stop("error in scree plot")
                }, finally = { 
                     cat("scree plot done\n")
                     dev.off()
                }     
            )
        }
        if (is.null(fa.numbers)) fa.numbers <- 1:ngroup
        lda.best.accuracy <- matrix(NA, nrow = length(fa.numbers), ncol = length(fa.factors))
        lda.best.p <- matrix(NA, nrow = length(fa.numbers), ncol = length(fa.factors))
        lda.best.model <- matrix("", nrow = length(fa.numbers), ncol = length(fa.factors))
    	lda.summary <- matrix(NA, nrow = length(fa.numbers), ncol = length(fa.factors))
        for(fa.number in fa.numbers) {
        	fa.result <- factor.pa.ginv(fa.dat, nfactor = fa.number, residuals = T, prerotate = T, rotate = "promax", scores = T)
        	fa.loading <- fa.result$loadings[1:nsample, , drop = F]
        	fa.score <- fa.result$scores
            dout <- file.path(dir.out, fa.number)
            dir.create(dout, F, T)
            floading <- file.path(dout, file.loading)
        	write.csv(fa.loading, file = floading)
            fscore <- file.path(dout, file.score)
        	write.csv(fa.score, file = fscore)
        	fa.combinations <- get.fa.comb(fa.number)
     		lda.predict <- get.lda.predict(fa.loading, fa.combinations, fa.factors, test.method = test.method, cv = cv) # test.method = "chisq.test", "mcnema.test"
            faccuracy <- file.path(dout, file.accuracy)
            write.csv(lda.predict$accuracy, file = faccuracy)
        	fp <- file.path(dout, file.p)
            write.csv(lda.predict$p, file = fp)
        	lda.best <- get.lda.best(lda.predict, fa.number, fa.loading, p.cutoff)
        	lda.best.accuracy[fa.number, ] <- lda.best$accuracy
        	lda.best.p[fa.number, ] <- lda.best$p
        	lda.best.model[fa.number, ] <- lda.best$model
        	}
    	rownames(lda.best.accuracy) <- paste("model", fa.numbers, sep = "")
    	colnames(lda.best.accuracy) <- names(fa.factors)
    	rownames(lda.best.p) <- paste("model", fa.numbers, sep = "")
    	colnames(lda.best.p) <- names(fa.factors)
    	rownames(lda.best.model) <- paste("model", fa.numbers, sep = "")
    	colnames(lda.best.model) <- names(fa.factors)
    	lda.summary <- get.lda.summary(lda.best.accuracy, lda.best.p, lda.best.model)
        fsummary <- file.path(dir.out, file.summary)
    	write.csv(lda.summary, file = fsummary)
    }
    
    get.fa.mols <- function(dir.in, file.score = "score.csv", file.scoreplot = "score.pdf", fa.number, sd.cutoff = 2, plot = TRUE, save = TRUE, idtypes = c("asis", "symbol"), id.maps = NULL, target.types = c("predicted", "validated", "both")) {
    	fscore <- file.path(dir.in, fa.number, file.score)
    	fa.score <- read.csv(fscore, header = T, check.names = F, row.names = 1)
    	fa.score <- as.matrix(fa.score)
    	fa.labels <- paste("F", c(1:fa.number), sep = "")
    	plot.fa.score.density <- function(fa.labels, fa.scores, fa.number) {
    		fa.score.density.ymaxs <- lapply(fa.labels, function(fa.label) {max(density(fa.scores[, fa.label])$y)})
    		fa.score.density.xmaxs <- lapply(fa.labels, function(fa.label) {max(density(fa.scores[, fa.label])$x)})
    		fa.score.density.xmins <- lapply(fa.labels, function(fa.label) {min(density(fa.scores[, fa.label])$x)})
    		fa.score.density.xmin <- min(unlist(fa.score.density.xmins))
    		fa.score.density.xmax <- max(unlist(fa.score.density.xmaxs))
    		fa.score.density.ymax.index <- which.max(fa.score.density.ymaxs)
    		fa.score.density.ymax <- max(unlist(fa.score.density.ymaxs))
    		plot(density(fa.scores[, fa.score.density.ymax.index]), xlim = c(fa.score.density.xmin, fa.score.density.xmax), main = paste("Factor score distributions of model", fa.number, sep = " "), type = "n")
    		cols <- rainbow(fa.number)
    		for (fa.num in 1:fa.number) {
    			col.index <- fa.num 
                fa.label <- fa.labels[col.index]
    			lines(density(fa.scores[, fa.label]), col = cols[col.index], lwd = 3)
    		}
    		legend(fa.score.density.xmin, fa.score.density.ymax, legend = fa.labels, col = cols, lty = rep(1, fa.num), lwd = rep(2, fa.num))
    	}
    	get.fa.mol <- function(fa.label, fa.scores, sd.cutoff) {
    		fa.score <- fa.scores[, fa.label]
    		fa.score.sd <- sd(fa.score)
    		mols.plus <- names(fa.score[fa.score >= fa.score.sd * sd.cutoff])
    		mols.minus <- names(fa.score[fa.score <= -fa.score.sd * sd.cutoff])
    		fa.label.plus <- paste(fa.label, "+", sep = '')
    		fa.label.minus <- paste(fa.label, "-", sep = '')
    		mols <- list(mols.plus, mols.minus)
    		names(mols) <- c(fa.label.plus, fa.label.minus)
    		return (mols)
    	}
        fscoreplot <- file.path(dir.in, fa.number, file.scoreplot)
        if (plot) {
            pdf(file = fscoreplot)
    	    tryCatch(plot.fa.score.density(fa.labels, fa.score, fa.number),
                     error = function(e) {
                         dev.off()
                         stop("error in plotting")
                     }, finally = {
                         dev.off()
                         cat("plot ok\n")
                     })
        }
    	mols <- lapply(fa.labels, FUN = get.fa.mol, fa.scores = fa.score, sd.cutoff = sd.cutoff)
        if (save) {
            lab.sd <- paste("sigma", sd.cutoff, sep = "_")
            if (is.null(id.maps) && (!all(idtypes %in% "asis"))) {
                stop("no id map")
            }
            mols <- do.call(c, mols)
            dout <- file.path(dir.in, fa.number, lab.sd)
            dir.create(dout, F, T)
            fmols.asis <- file.path(dout, "mols_asis.csv")
            write.csv(list2dataframe(mols, ""), file = fmols.asis, row.names = F)
            for (target.type in target.types) {
                file.mols.symbol <- paste("mols_symbol_", target.type, ".csv", sep = "")
                fmols.symbol <- file.path(dout, file.mols.symbol)
                write.csv(list2dataframe(lapply(mols, function(mol, id.maps, target.type) unique(unlist(id.maps[[target.type]][mol], use.names = F)), id.maps = id.maps, target.type = target.type), ""), file = fmols.symbol, row.names = F)
            }
        }
    }
    
    ## v2 output all molecules in factors
    get.fa.scores <- function(dir.in, file.score = "score.csv", fa.number, sort = TRUE) {
    	fscore <- file.path(dir.in, fa.number, file.score)
    	fa.score <- read.csv(fscore, header = T, check.names = F, row.names = 1)
    	fa.score <- as.matrix(fa.score)
    	fa.labels <- paste("F", c(1:fa.number), sep = "")
    	get.fa.score <- function(fa.label, fa.scores, sort = TRUE) {
    		fa.score <- fa.scores[, fa.label]
            if (isTRUE(sort)) fa.score <- sort(fa.score, decreasing = T)
            fa.score
    	}
    	scores <- sapply(fa.labels, FUN = get.fa.score, fa.scores = fa.score, sort = sort, USE.NAMES = T, simplify = F)
    }
    
    
    get.fa.sigmols <- function(scores, fa.labels, sd.cutoff = 2, sort = TRUE) {
        get.fa.sigmol <- function(fa.label, scores, sd.cutoff, sort = TRUE) {
            score <- scores[[fa.label]]
            if (isTRUE(sort)) score <- sort(score, decreasing = T)
            score.sd <- sd(score, na.rm = T)
            mol <- names(score)
            sigmol.plus <- mol[score >= sd.cutoff * score.sd]
            sigmol.minus <- mol[score <= -sd.cutoff * score.sd]
            sigmol <- list(sigmol.plus, sigmol.minus)
            names(sigmol) <- paste(fa.label, c("+", "-"), sep = "")
            sigmol
        }
        sigmols <- sapply(fa.labels, get.fa.sigmol, scores = scores, sd.cutoff = sd.cutoff, sort = sort, USE.NAMES = T, simplify = F)
    }

       for (obj in ls()) 
       assign(obj, get(obj), envir = Enrichments)
}
