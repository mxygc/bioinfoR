if (!exists("Enrichments") || !is.environment(Enrichments)) Enrichments <- new.env(parent = emptyenv())

local({
	
    get.genesets <- function(file, sets = NULL, id = "offsym", type = "basic", map = NULL, uppercase = F) {
    # the sheet is first letter capitalized only
        genesets <- read.csv(file = file, head = T, as.is = T, check.names = F)
        genesets <- lapply(genesets, function(x, uppercase) { x <- x[x != ""]; if (uppercase) x <- toupper(x); x <- unique(x) }, uppercase = uppercase)
        if (is.null(sets)) sets <- sub(".MAJOR$", "", grep("MAJOR$", names(genesets), value = T))
        sets.basic <- sets
        sets.paralog <- paste(sets, "MAJOR", sep = ".")
        genesets.basic <- sapply(sets.basic, function(s, genesets) list(basic = genesets[[s]], paralog = character(0)), genesets = genesets, simplify = F, USE.NAMES = T)
        genesets.paralog <- sapply(sets.basic, function(s, genesets) list(basic = genesets[[s]], paralog = setdiff(genesets[[paste(s, "MAJOR", sep = ".")]], genesets[[s]])), genesets = genesets, simplify = F, USE.NAMES = T)
        names(genesets.paralog) <- sets.paralog
        genesets <- c(genesets.basic, genesets.paralog)
    } 
    
    get.hpg <- function(s, bg, l, n = 3000, seed = 1, upper.case = T) {
        get.perm.lists <- function(bg, g, n){
        	lg <- length(g)
        	lbg <- length(bg)
        	perms.lists <- replicate(n, bg[sample(lbg, replace = FALSE, size = lg)])
            if (lg == 1) perms.lists <- matrix(perms.lists, nrow = lg)
            perms.lists
        }
        get.hypergeom.p <- function(bg, l, s) {
        	q <- sum(l %in% s)
        	m <- sum(bg %in% s)
        	n <- sum(!bg %in% s)
        	k <- length(l)
        	phyper(q = q - 1, m = m, n = n, k = k, lower.tail = F)
        }
        if (is.list(s)) {
            b <- s$basic
            e <- s$paralog
            i <- intersect(b, e)
            if (any(i)) {
                if (!setequal(i, b)) 
                    warning("the extened list has some but not whole overlapping genes with the basic one")
                s <- e
                e <- setdiff(e, b)
            }
            else s <- union(b, e)
        }
        else {
            b <- s
            e <- character()
        }
      	bg <- unique(bg)
       	l <- unique(l)
       	s <- unique(s)
        if (upper.case) {
           	bg <- toupper(bg)
        	l <- toupper(l)
        	s <- toupper(s)
            b <- toupper(b)
            e <- toupper(e)
        }
    	ib <- paste(intersect(b, l), collapse = ";")
        ie <- paste(intersect(e, l), collapse = ";")
        if (length(l) > 0) {
        	p.real <- get.hypergeom.p(bg = bg, l = l, s = s)
        	if (!missing(seed)) seed = seed
            set.seed(seed)
        	l.perms <- get.perm.lists(bg = bg, g = l, n = n)
            if (n >= 1) {
            	colnames(l.perms) <- paste("perm_", 1:n, sep = "")
                p.perms <- apply(l.perms, 2, function(l.perm, bg, s) get.hypergeom.p(bg = bg, l = l.perm, s = s), bg = bg, s = s)
            } else p.perms <- NULL
        } else {
            p.real <- NA
            p.perms <- rep(NA, n)
        }
    	list(p = p.real, p.perms = p.perms, basic = ib, paralog = ie)
    }
    
    get.hypergeom <- function(genes, mods = NULL, genesets, bg, n = 3000, seed = 1, correct.methods = c("BH", "BY", "Holm", "Hochberg", "Hommel", "Bonferroni", "bootstrap"), correct = "row", sep.prlg = FALSE, upper.case = TRUE, fun.get.hpg = get.hpg) {
    # for each comparison: genes = list(mod_1 = character(), mod_2 = character(), ...)
    # if paralog = F, intersection = c(intersect(l, s))
    # else, intersection = c(intersect(l, s), intersect(e, s))
    # correct defines way to do multiple correction: NULL (no correction), by row, col, or both
        if (is.null(mods)) mods <- names(genes)
        hpg <- NULL
        hpgs <- NULL
        nrow <- length(mods)
        ncol <- length(genesets)
        sets <- names(genesets)
        p.perms <- list()
        p <- matrix(NA, nrow, ncol, dimnames = list(mods, sets))
        basic <- matrix("", nrow, ncol, dimnames = list(mods, sets))
        paralog <- basic
        overlap <- basic
        hpg.mc <- NULL
        get.mc <- function(method, p, mods, sets, p.perms = NULL, correct = "row") {
            .methods <- c(BH = "BH", BY = "BY", Holm = "holm", Hochberg = "hochberg", Hommel = "hommel", Bonferroni = "bonferroni", bootstrap = "bootstrap")
            m <- .methods[method]
            if (m != "bootstrap") {
                if (correct == "both") {
                    mc <- as.vector(p)
                    mc <- p.adjust(mc, method = m)
                    mc <- matrix(mc, length(mods), length(sets))
                    colnames(mc) <- sets
                } 
                else if (correct == "row") {
                    mc <- t(apply(p, 1, p.adjust, method = m))
                }
                else if (correct == "col") {
                    mc <- apply(p, 2, p.adjust, method = m)
                }
            } 
            else {
                if (is.null(p.perms)) 
                    stop("permutation derived p values don't exist")
                mc <- p
                p.perm <- list()
                if (correct == "both")
                    p.perm <- unlist(p.perms)
                else if (correct == "row")
                    for (mod in mods) p.perm[[mod]] <- unlist(p.perms[[mod]])
                else if (correct == "col") {
                    for (set in sets) for (mod in mods) p.perm[[set]] <- c(p.perm[[set]], p.perms[[mod]][[set]])
                }
                for (mod in mods) 
                    for (set in sets) {
                        if (correct == "both")
                            mc[mod, set] <- mean(p.perm <= p[mod, set], na.rm = T)
                        else if (correct == "row")
                            mc[mod, set] <- mean(p.perm[[mod]] <= p[mod, set], na.rm = T)
                        else if (correct == "col")
                            mc[mod, set] <- mean(p.perm[[set]] <= p[mod, set], na.rm = T)
                    }
            }
            mc
        }
        if (is.null(correct)) n <- 0
        for (mod in mods) {
            gene <- genes[[mod]]
            cat(mod, "\n")
            hpg <- c(hpg, list(sapply(genesets, fun.get.hpg, bg = bg, l = gene, n = n, seed = seed, upper.case = upper.case, USE.NAMES = T, simplify = F)))
            names(hpg)[length(hpg)] <- mod
        }
        for (mod in mods) {
            for (set in sets) {
                p[mod, set] <- hpg[[mod]][[set]][["p"]]
                basic[mod, set] <- hpg[[mod]][[set]][["basic"]]
                paralog[mod, set] <- hpg[[mod]][[set]][["paralog"]]
                overlap[mod, set] <- ifelse(paralog[mod, set] == "", basic[mod, set], paste(basic[mod, set], paralog[mod, set], sep = ";"))
                if (!is.null(correct)) p.perms[[mod]][[set]] <- hpg[[mod]][[set]][["p.perms"]]
            }
        } 
        if (!is.null((correct))) hpg.mc <- sapply(correct.methods, get.mc, p = p, mods = mods, sets = sets, p.perms = p.perms, correct = correct, USE.NAMES = T, simplify = F)
        if (sep.prlg) hpgs <- c(list(p = p), hpg.mc, list(basic = basic), list(paralog = paralog))
        else hpgs <- c(list(p = p), hpg.mc, list(overlap = overlap))
        rownames <- paste(rep(mods, each = length(hpgs)), names(hpgs), sep = "_")
        hypergeom <- matrix("", length(rownames), length(sets), dimnames = list(rownames, sets))
        for (mod in mods)
            for (set in sets)
                for (r in names(hpgs)) {
                    rowname <- paste(mod, r, sep = "_")
                    hypergeom[rowname, set] <- hpgs[[r]][mod, set]
                }
        hypergeom
    }
    
    # use parallel computing for percs
    get.hypergeoms <- function(genes, nregs = NULL, regs = NULL, percs = NULL, mods = NULL, genesets, bg, n = 3000, seed = 1, correct.methods = c("BH", "BY", "Holm", "Hochberg", "Hommel", "Bonferroni", "bootstrap"), correct = "row", sep.prlg = FALSE, parallel = FALSE, cl = NULL, upper.case = TRUE, fun.get.hypergeom = get.hypergeom, fun.get.hpg = get.hpg) {
        if (is.null(nregs)) nregs <- names(genes)
        if (is.null(regs)) regs <- names(genes[[1]])
        if (is.null(percs)) percs <- names(genes[[1]][[1]])
        if (is.null(mods)) mods <- names(genes[[1]][[1]][[1]])
        hypergeoms <- NULL
        if (parallel) {
            require("parallel")
            cl <- makeCluster(length(percs))
        }
        for (nreg in nregs)
            for (reg in regs) {
                cat(nreg, reg, "\n")
                gene <- genes[[nreg]][[reg]]
                if (parallel && !is.null(cl)) {
                    hypergeoms[[nreg]][[reg]] <- parSapply(cl, percs, function(perc, gene, mods, genesets, bg, n, seed, correct.methods, correct, sep.prlg, upper.case, fun.get.hypergeom, fun.get.hpg) { g <- gene[[perc]]; fun.get.hypergeom(g, mods = mods, genesets = genesets, bg = bg, n = n, seed = seed, correct.methods = correct.methods, correct = correct, sep.prlg = sep.prlg, upper.case = upper.case, fun.get.hpg = fun.get.hpg) }, gene = gene, mods = mods, genesets = genesets, bg = bg, n = n, seed = seed, correct.methods = correct.methods, correct = correct, sep.prlg = sep.prlg, upper.case = upper.case, fun.get.hypergeom = fun.get.hypergeom, fun.get.hpg = fun.get.hpg, USE.NAMES = T, simplify = F)
                } else {
                    hypergeoms[[nreg]][[reg]] <- sapply(percs, FUN = function(perc, gene, mods, genesets, bg, n, seed, correct.methods, correct, sep.prlg, upper.case, fun.get.hpg) { cat(perc, "\n"); g <- gene[[perc]]; get.hypergeom(g, mods = mods, genesets = genesets, bg = bg, n = n, seed = seed, correct.methods = correct.methods, correct = correct, sep.prlg = sep.prlg, upper.case = upper.case, fun.get.hpg = fun.get.hpg) }, gene = gene, mods = mods, genesets = genesets, bg = bg, n = n, seed = seed, correct.methods = correct.methods, correct = "row", sep.prlg = sep.prlg, upper.case = upper.case, fun.get.hpg = fun.get.hpg, USE.NAMES = T, simplify = F)
                    }
                }
        if (parallel && !is.null(cl)) stopCluster(cl)
        hypergeoms
    }
    get.hypergeoms.rep <- function(genes, reps = NULL, nregs = NULL, regs = NULL, percs = NULL, mods = NULL, genesets, bg, n = 3000, seed = 1, correct.methods = c("BH", "BY", "Holm", "Hochberg", "Hommel", "Bonferroni", "bootstrap"), correct = "row", sep.prlg = FALSE, parallel = TRUE, cl = NULL, upper.case = TRUE, get.hypergeom = get.hypergeom, get.hpg = get.hpg) {
        if (is.null(reps)) reps <- names(genes)
        if (is.null(nregs)) nregs <- names(genes[[1]])
        if (is.null(regs)) regs <- names(genes[[1]][[1]])
        if (is.null(percs)) percs <- names(genes[[1]][[1]][[1]])
        if (is.null(mods)) mods <- names(genes[[1]][[1]][[1]][[1]])
        hypergeoms <- NULL
        for (rp in reps) {
            cat(rp, "\n")
            gene <- genes[[rp]]
            hypergeoms[[rp]] <- get.hypergeoms(gene, nregs = nregs, regs = regs, percs = percs, mods = mods, genesets = genesets, bg = bg, n = n, seed = seed, correct.methods = correct.methods, correct = correct, sep.prlg = sep.prlg, parallel = parallel, cl = cl, upper.case = upper.case, get.hypergeom = get.hypergeom, get.hpg = get.hpg)
            }
        hypergeoms
    }
    
    get.hpg.sigset <- function(hpg, mods = NULL, pval = 0.05, overlap = F) {
        if (is.null(mods)) {
            rnames <- rownames(hpg)
            mods <- paste("mod", sort(as.numeric(unique(unlist(lapply(strsplit(rnames, "_"), function(x) x[2]))))), sep = "_")
        }
        get.hpg.ss <- function(mod, hpg, pval, overlap) {
            sn <- colnames(hpg)
            rp <- paste(mod, "p", sep = "_")
            ri <- paste(mod, "overlap", sep = "_")
            p <- hpg[rp, ]
            i <- hpg[ri, ]
            idx <- which(p < pval)
    #       if (any(idx)) sigset <- paste(sn[idx], " (", p[idx], ":", i[idx], ")", collapse = "; ", sep = "")
            if (overlap) {
                if (any(idx)) sigset <- paste(sn[idx], " (", i[idx], ")", collapse = "; ", sep = "")
                else sigset <- ""
            } else {
                if (any(idx)) sigset <- paste(sn[idx], collapse = "; ")
                else sigset <- ""
            }
            sigset
        }
        sigset <- sapply(mods, get.hpg.ss, hpg = hpg, pval = pval, overlap = overlap, simplify = T, USE.NAMES = T)
    }
    get.hpg.sigsets <- function(hpgs, nregs = NULL, regs = NULL, percs = NULL, mods = NULL, pval = 0.05, overlap = F) {
        if (is.null(nregs)) nregs <- names(hpgs)
        if (is.null(regs)) regs <- names(hpgs[[1]])
        if (is.null(percs)) percs <- names(hpgs[[1]][[1]])
        if (is.null(mods)) {
            rnames <- rownames(hpgs[[1]][[1]][[1]])
            mods <- paste("mod", sort(as.numeric(unique(unlist(lapply(strsplit(rnames, "_"), function(x) x[2]))))), sep = "_")
        }
        sigsets <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    sigsets[[nreg]][[reg]][[perc]] <- get.hpg.sigset(hpgs[[nreg]][[reg]][[perc]], mods = mods, pval = pval, overlap = overlap)
                }
        sigsets
    }
    
    get.hpg.significance <- function(hpg, mods = NULL, pval = 0.05, binarize = T) {
        if (is.null(mods)) {
            rnames <- rownames(hpg)
            mods <- paste("mod", sort(as.numeric(unique(unlist(lapply(strsplit(rnames, "_"), function(x) x[2]))))), sep = "_")
        }
        get.hpg.sig <- function(mod, hpg, pval = 0.05, binarize = T) {
            sn <- colnames(hpg)
            rp <- paste(mod, "p", sep = "_")
            p <- as.numeric(hpg[rp, ])
            p[is.na(p)] <- 1
            if (binarize) p <- 1 - as.numeric(p >= pval)
            'names<-'(p, sn)
        }
        significance <- t(sapply(mods, get.hpg.sig, hpg = hpg, pval = pval, binarize = binarize, simplify = T, USE.NAMES = T))
    }
    get.hpg.significances <- function(hpgs, nregs = NULL, regs = NULL, percs = NULL, mods = NULL, pval = 0.05, binarize = T) {
        if (is.null(nregs)) nregs <- names(hpgs)
        if (is.null(regs)) regs <- names(hpgs[[1]])
        if (is.null(percs)) percs <- names(hpgs[[1]][[1]])
        if (is.null(mods)) {
            rnames <- rownames(hpgs[[1]][[1]][[1]])
            mods <- paste("mod", sort(as.numeric(unique(unlist(lapply(strsplit(rnames, "_"), function(x) x[2]))))), sep = "_")
        }
        sigsets <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    sigsets[[nreg]][[reg]][[perc]] <- get.hpg.significance(hpgs[[nreg]][[reg]][[perc]], mods = mods, pval = pval, binarize = binarize)
                }
        sigsets
    }
    get.hpg.significances.rep <- function(hpgs, reps, nregs = NULL, regs = NULL, percs = NULL, mods = NULL, pval = 0.05, binarize = T) {
        if (is.null(reps)) reps <- names(hpgs)
        if (is.null(nregs)) nregs <- names(hpgs[[1]])
        if (is.null(regs)) regs <- names(hpgs[[1]][[1]])
        if (is.null(percs)) percs <- names(hpgs[[1]][[1]][[1]])
        if (is.null(mods)) {
            rnames <- rownames(hpgs[[1]][[1]][[1]][[1]])
            mods <- paste("mod", sort(as.numeric(unique(unlist(lapply(strsplit(rnames, "_"), function(x) x[2]))))), sep = "_")
        }
        sigsets <- NULL
        for (rp in reps) {
            cat(rp, "\n")
            sigsets[[rp]] <- get.hpg.significances(hpgs = hpgs[[rp]], nregs = nregs, regs = regs, percs = percs, mods = mods, pval = pval, binarize = binarize)
        }
        sigsets
    }
    
    # ============================================================
    # GSEA preranked analysis
    # ============================================================
    get.gsea.gmx <- function(genesets, dir = ".", filename = "gsea", upper.case = T) {
        if (upper.case) genesets <- rapply(genesets, function(x) unique(toupper(x)), how = "replace")
        get.gmx.m <- function(genesets, type = "basic") {
            if (type == "basic") {
                m <- sapply(1:length(genesets), function(i) c("", genesets[[i]]$basic), simplify = F)
                names(m) <- names(genesets)
            }
            else if (type == "paralog") {
                m <- sapply(1:length(genesets), function(i) c("", union(genesets[[i]]$basic, genesets[[i]]$paralog)), simplify = T)
                names(m) <- paste(names(genesets), "MAJOR", sep = ".")
            }   
            else if (is.null(type)) {
                m <- sapply(1:length(genesets), function(i) c(names(genesets)[i], "", genesets[[i]]), simplify = F)
                names(m) <- names(genesets)
            }
            m <- list2dataframe(m)
        }
        gmx <- list()
        if (is.list(genesets[[1]])) {
            basic <- get.gmx.m(genesets, type = "basic")
            paralog <- get.gmx.m(genesets, type = "paralog")
        } else {
            gene <- get.gmx.m(genesets, type = NULL)
        }
        dir.create(dir, recursive = T, showWarnings = F)
        if (is.list(genesets[[1]])) {
            for (type in c("basic", "paralog")) {
                m <- get.gmx.m(genesets, type)
                gmx[[type]] <- m
                file <- paste(filename, type, sep = "_")
                file <- paste(file, "gmx", sep = ".")
                fout <- paste(dir, file, sep = "/")
                write.table(m, file = fout, col.names = T, row.names = F, quote = F, sep = "\t")
            }
        } else {
            m <- get.gmx.m(genesets, NULL)
            gmx <- m
            file <- paste(file, "gmx", sep = ".")
            fout <- paste(dir, file, sep = "/")
            write.table(m, file = fout, col.names = T, row.names = F, quote = F, sep = "\t")
        }
        invisible(gmx)
    }
    # the input is regweight.fpr, in which the duplicate genes are already filtered
    get.gsea.rnk <- function(regweight, mods = NULL, id.map, upper.case = TRUE) {
        if (is.null(mods)) mods <- names(regweight)
        collapse.rw <- function(rw, map) {
            id0 <- names(rw)
            id1 <- map[id0]
            l0 <- length(id0)
            l1 <- length(unlist(id1))
            rw1 <- data.frame(reg = rep("", l1), weight = rep(NA, l1), stringsAsFactors = F)
            j <- 1
            for (i in 1:l0) {
                if (i %% 1000 == 0) cat(i, "\n")
                id <- id1[[i]] 
                l <- length(id)
                if (l == 0) next
                k <- j + l - 1
                rw1[j:k, 1] <- id
                rw1[j:k, 2] <- rep(rw[i], l)
                j <- k + 1 
            } # for identical gene correspoding to multiple IDs and weights, select the max weight
            rw2 <- tapply(rw1[, 2], rw1[, 1], max, simplify = T)
        }
        rnk <- sapply(mods, function(mod, rw, map, uc) {
            cat(mod, "\n")
            rw <- regweight[[mod]]
            rw <- collapse.rw(rw, map)
            if (uc) rw <- "names<-"(rw, toupper(names(rw)))
            rw <- sort(rw, decreasing = T)
            }, rw = regweight, map = id.map, uc = upper.case, USE.NAMES = T, simplify = T)
    }
    
    get.gsea.rnks <- function(regweights, nregs = NULL, regs = NULL, percs = "perc_100", mods = NULL, id.map, upper.case = T) {
        if (is.null(nregs)) nregs <- names(regweights)
        if (is.null(regs)) regs <- names(regweights[[1]])
        if (is.null(percs)) percs <- names(regweights[[1]][[1]])
        if (is.null(mods)) mods <- names(regweights[[1]][[1]][[1]])
        rnks <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    rnks[[nreg]][[reg]][[perc]] <- list()
                    rnks[[nreg]][[reg]][[perc]] <- get.gsea.rnk(regweights[[nreg]][[reg]][[perc]], mods = mods, id.map = id.map, upper.case = upper.case)
                }
        rnks
    }
    get.gsea.rnk.fpr <- function(dir, nreg, reg, perc, convert = FALSE, id.map = NULL, upper.case = T) {
        file <- paste("reg", "p", perc, sep = "_")
        file <- paste(file, "txt", sep = ".")
        file <- paste(dir, paste("nreg", nreg, sep = "_"), reg, file, sep = "/")
        regweight <- read.csv(file, sep = "\t", head = F, as.is = T)
        regweight <- by(data = regweight, INDICES = regweight$V2, FUN = function(x) { y <- tapply(x$V3, x$V1, FUN = max, simplify = T); z <- as.vector(y); names(z)  <- names(y); z }, simplify = F)
        names(regweight) <- paste("mod", names(regweight), sep = "_") 
        regweight <- lapply(regweight, function(x) sort(unlist(x), decreasing = T))
        if (convert && !is.null(id.map)) {
            count <- 0
            regweight <- lapply(regweight, function(rw, map) {
                count <<- count + 1
                cat("mod_", count, "\n", sep = "")
                id0 <- names(rw)
                id1 <- map[id0]
                l0 <- length(id0)
                l1 <- length(unlist(id1))
                rw1 <- data.frame(reg = rep("", l1), weight = rep(NA, l1), stringsAsFactors = F)
                j <- 1
                sapply(1:l0, function(i, rw, id1) { 
                       if (i %% 1000 == 0) cat("    ", i, "\n")
                       id <- id1[[i]]; 
                       l <- length(id); 
                       if (l == 0) next
                       k <- j + l - 1; 
                       rw1[j:k, 1] <<- id; 
                       rw1[j:k, 2] <<- rep(rw[i], l); 
                       j <<- k + 1 }, rw = rw, id1 = id1, simplify = T)
                rw2 <- tapply(rw1[, 2], rw1[, 1], max, simplify = T)
            }, map = id.map)
        }
        if (upper.case) regweight <- lapply(regweight, function(x) { names(x) <- toupper(names(x)); x})
        regweight
    }
    
    gen.gsea.prerank <- function(jar = "/picb/clingenet/Applications/GSEA/gsea2-2.0.10.jar", rnk, file.gmx, xmx = 1024, collapse = FALSE, nperm = 3000, rpt.label, seed = NULL, set.max = 500, set.min = 5, dir.out, dir.log, script, run = FALSE) {
        dir.create(path = dir.out, recursive = T, showWarnings = F)
        file.out <- paste(rpt.label, "rnk", sep = ".")
        fout <- paste(dir.out, file.out, sep = "/")
        file.log <- paste(rpt.label, "log", sep = ".")
        flog <- paste(dir.log, file.log, sep = "/")
        write.table(rnk, file = fout, row.names = T, col.names = F, sep = "\t", quote = F)
    	cmd = paste(
            "java", "-cp", jar,
    		paste("-Xmx", xmx, sep = ""), 
    		"xtools.gsea.GseaPreranked",
            "-rnk", fout,
    		"-gmx", file.gmx,
    		"-collapse", ifelse(collapse, "true", "false"),
    		"-nperm", nperm,
    		"-rpt_label", rpt.label, 
            ifelse(is.null(seed), "", paste("-rnd_seed", seed)), 
    		"-set_max", set.max,
    		"-set_min", set.min,
    		"-out", dir.out, 
            ">", flog, 
            "2>&1 &")
        cat(cmd, file = script)
    	if (run) {
            exit <- system(cmd)
            if (exit != 0) warning(rpt.label, ": has problem")
            else T
        }
    }
    
    gen.gsea.preranks <- function(jar = "/picb/clingenet/Applications/GSEA/gsea2-2.0.10.jar", rnks, rpt.labels = NULL, file.gmx, xmx = "1024m", collapse = FALSE, nperm = 3000, seed = NULL, set.max = 500, set.min = 5, dir.out, dir.log, dir.script, run = FALSE, procmax = 10) {
        dir.create(dir.log, recursive = T, showWarnings = F)
        dir.create(dir.script, recursive = T, showWarnings = F)
        if (is.null(rpt.labels)) rpt.labels <- names(rnks)
        i <- 0
        for (label in rpt.labels) {
            i <- i + 1
            rnk <- rnks[[label]]
            file.script <- paste("gsea", label, sep = "_")
            file.script <- paste(file.script, "bash", sep = ".")
            fscript <- paste(dir.script, file.script, sep = "/")
            running <- gen.gsea.prerank(jar = jar, rnk = rnk, rpt.label = label, file.gmx = file.gmx, xmx = xmx, collapse = collapse, nperm = nperm, seed = seed, set.max = set.max, set.min = set.min, dir.out = dir.out, dir.log = dir.log, script = fscript, run = run)
            if (run && running && i %% 10 == 0) {
                while (running) {
                    ps <- system2("ps", args = "-u $USER", stdout = T, stderr = F)
                    if (any(grepl("java", ps))) Sys.sleep(5)
                    else running <- F
                }
            }
        }
    }
    gen.gseas.prerank <- function(jar = "/picb/clingenet/Applications/GSEA/gsea2-2.0.10.jar", ranks, nregs = "nreg_all", regs = NULL, percs = "perc_100", mods = NULL, rpt.labels = NULL, file.gmx, xmx = "1024m", collapse = FALSE, nperm = 3000, seed = NULL, set.max = 500, set.min = 5, dir.out, dir.log, dir.script, run = FALSE, procmax = 10) {
        if (is.null(nregs)) nregs <- names(ranks)
        if (is.null(regs)) regs <- names(ranks[[1]])
        if (is.null(percs)) percs <- names(ranks[[1]][[1]])
        if (is.null(mods)) mods <- names(ranks[[1]][[1]][[1]])
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    rnks <- ranks[[nreg]][[reg]][[perc]]
                    dout <- paste(dir.out, nreg, reg, perc, sep = "/")
                    dlog <- paste(dir.log, nreg, reg, perc, sep = "/")
                    dscript <- paste(dir.script, nreg, reg, perc, sep = "/")
                    gen.gsea.preranks(jar = jar, rnks = rnks, rpt.labels = NULL, file.gmx = file.gmx, xmx = xmx, collapse = collapse, nperm = nperm, seed = seed, set.max = set.max, set.min = set.min, dir.out = dout, dir.log = dlog, dir.script = dscript, run = run, procmax = procmax)
                }
    }
    read.gsea <- function(dir.in, mods, regs = c("up", "down"), sets, setsizes) {
        # an ugly fix for the historical variable name
        compares <- mods
        get.suffix <- function(dir) {
            sf <- strsplit(basename(list.dirs(dir, recursive = F, full.names = F)),"\\.")
            names(sf) <- sapply(sf, function(x) x[1])
            sf <- sapply(sf, function(x) x[3], USE.NAMES = T, simplify = T)
        }
        get.cmp.failed <- function(compares, suffix) {
            cmp.f <- setdiff(compares, names(suffix))
        }
        get.ce <- function(set, dir) {
            f <- paste(set, "xls", sep = ".")
            f <- paste(dir, f, sep = "/")
            m <- read.csv(file = f, head = T, sep = "\t", as.is = T, check.names = F)
            core <- paste(m[m[, "CORE ENRICHMENT"] == "Yes", "PROBE"], collapse = ";")
            overlap <- paste(m[, "PROBE"], collapse = ";")
            c(overlap = overlap, core = core)
        }
        get.set.size <- function(dir) {
            f <- paste(dir, "gene_set_sizes.xls", sep = "/")
            m <- read.csv(file = f, head = T, sep = "\t", as.is = T)
            rownames(m) <- m[, "NAME"]
            m <- m[, -ncol(m)] # remove the last column
        }
        get.stat <- function(cmp, reg, dir, suf) {
            REG <- c(up = "pos", down = "neg")
            f <- paste("gsea", "report", "for", "na", REG[reg], suf, sep = "_")
            f <- paste(f, "xls", sep = ".")
            f <- paste(dir, f, sep = "/")
            m <- read.csv(file = f, head = T, sep = "\t", as.is = T, check.names = F)
            set.size <- get.set.size(dir)
            sets <- rownames(set.size)
            n <- length(sets)
            if (nrow(m) == 0) {
                stat <- data.frame("NAME" = sets,
                                   "SIZE" = set.size[, "AFTER.RESTRICTING.TO.DATASET"],
                                   "ES" = rep(NA, n),
                                   "NES" = rep(NA, n),
                                   "NOM p-val" = rep(NA, n),
                                   "FDR q-val" = rep(NA, n),
                                   "FWER p-val" = rep(NA, n),
                                   "overlap" = rep("", n), 
                                   "core" = rep("", n), 
                                   stringsAsFactors = F)
                rownames(stat) <- sets
            } else {
                stat <- m[, c("NAME", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val")]
                ce <- t(sapply(stat[, "NAME"], get.ce, dir = dir, USE.NAMES = T, simplify = T))
                rownames(stat) <- stat[, "NAME"]
                stat <- data.frame(stat, ce, stringsAsFactors = F)
                stat <- t(sapply(sets, function(s) {
                            if (s %in% stat[, "NAME"]) unlist(stat[s, ]) #use unlist to coerce the slice of a dataframe - which is a 1xn sub dataframe - to become a true vector (say A)
                            else c("NAME" = s, 
                                   "SIZE" = set.size[s, "AFTER.RESTRICTING.TO.DATASET"],
                                   "ES" = NA,
                                   "NES" = NA,
                                   "NOM p-val" = NA,
                                   "FDR q-val" = NA,
                                   "FWER p-val" = NA,
                                   "overlap" = "",
                                   "core" = "") # this is a vector anyway (say B)
                                           }, USE.NAMES = T, simplify = T))
    # the final assembly of As and Bs won't be screwed up if each is the same type (vector)
            }
            stat <- stat[, -1]
            data.frame("ORIGINAL.SIZE" = set.size[sets, "ORIGINAL.SIZE"], stat[sets, ], stringsAsFactors = F)
        }
        get.stats <- function(cmps, regs, dir, sets, setsizes) {
            suffix <- get.suffix(dir)
            cmps.f <- get.cmp.failed(compares, suffix)
            stats <- list()
            n <- length(sets)
            for (cmp in cmps) {
                cat(cmp, "\n")
                if (!(cmp %in% cmps.f)) {
                    suf <- suffix[cmp]
                    d <- paste(cmp, "GseaPreranked", suf, sep = ".")
                    d <- paste(dir, d, sep = "/")
                    for (reg in regs)
                        stats[[cmp]][[reg]] <- get.stat(cmp = cmp, reg = reg, dir = d, suf = suf)
                } else {
                    for (reg in regs) {
                        stat <- data.frame(
                               "NAME" = sets,
                               "ORIGINAL.SIZE" = setsizes,
                               "SIZE" = rep(0, n), 
                               "ES" = rep(NA, n),
                               "NES" = rep(NA, n),
                               "NOM p-val" = rep(NA, n),
                               "FDR q-val" = rep(NA, n),
                               "FWER p-val" = rep(NA, n),
                               "overlap" = rep("", n), 
                               "core" = rep("", n), 
                               stringsAsFactors = F)
                        rownames(stat) <- sets
                        stats[[cmp]][[reg]] <- stat
                    }
                }
            }
            stats
        }
        stats <- get.stats(cmps = compares, regs = regs, dir = dir.in, sets = sets, setsizes = setsizes)
    }
    read.gseas <- function(dir.in, nregs, regs, percs, mods, regulations = c("up", "down"), sets, setsizes) {
        gseas <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    din <- paste(dir.in, nreg, reg, perc, sep = "/")
                    gseas[[nreg]][[reg]][[perc]] <- read.gsea(din, mods, regulations, sets, setsizes)
                }
        gseas
    }
    
    get.gsea.sigset <- function(gsea, mods, pval = 0.05, qval = NULL, fwer = NULL, n = 2) {
        compares <- mods
        sigset <- matrix("", ncol = 2, nrow = length(compares))
        rownames(sigset) <- compares
        colnames(sigset) <- c("up", "down")
        get.gs <- function(cmp, gsea, pval = 0.05, qval = NULL, fwer = NULL, n = 2) {
            cat(cmp, "\n")
            stat.up <- gsea[[cmp]][["up"]]
            stat.down <- gsea[[cmp]][["down"]]
            flag.up <- rep(T, nrow(stat.up))
            if (!is.null(pval)) flag.up <- flag.up & stat.up[, "NOM.p.val"] < pval
            if (!is.null(qval)) flag.up <- flag.up & stat.up[, "FDR.q.val"] <= qval
            if (!is.null(fwer)) flag.up <- flag.up & stat.up[, "FWER.p.val"] <= fwer
            if (!is.null(n)) stat.up[, "SIZE"] >= n
            flag.up[is.na(flag.up)] <- F
            if(sum(flag.up) == 0) gs.up <- ""
            else {
                ss.up <- rownames(stat.up)[flag.up]
                gs.up <- sapply(ss.up, function(s) { paste(s, " (", stat.up[s, "core"], ")", sep = "") }, USE.NAMES = T)
                gs.up <- unlist(gs.up)
                gs.up <- paste(gs.up, collapse = "; ")
            }
            flag.down <- rep(T, nrow(stat.down))
            if (!is.null(pval)) flag.down <- flag.down & stat.down[, "NOM.p.val"] < pval
            if (!is.null(qval)) flag.down <- flag.down & stat.down[, "FDR.q.val"] <= qval
            if (!is.null(fwer)) flag.down <- flag.down & stat.down[, "FWER.p.val"] <= fwer
            if (!is.null(n)) stat.down[, "SIZE"] >= n
            flag.down[is.na(flag.down)] <- F
            if(sum(flag.down) == 0) gs.down <- ""
            else {
                ss.down <- rownames(stat.down)[flag.down]
                gs.down <- sapply(ss.down, function(s) { paste(s, " (", stat.down[s, "core"], ")", sep = "") }, USE.NAMES = T)
                gs.down <- unlist(gs.down)
                gs.down <- paste(gs.down, collapse = "; ")
            }
            c(up = gs.up, down = gs.down)
        }
        for (cmp in compares) {
            sigset[cmp, ] <- get.gs(cmp = cmp, gsea = gsea, pval = pval, qval = qval, fwer = fwer, n = n)
        }
        sigset
        #as.data.frame(sigsets, stringsAsFactors = F)
    }
    
    get.gsea.sigsets <- function(gseas, nregs = "nreg_all", regs, percs = "perc_100", mods, pval = 0.05, qval = NULL, fwer = NULL, n = 2) {
        sigsets <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    gsea <- gseas[[nreg]][[reg]][[perc]]
                    sigsets[[nreg]][[reg]][[perc]] <- list()
                    sigsets[[nreg]][[reg]][[perc]] <- get.gsea.sigset(gsea, mods, pval, qval, fwer, n)
                }
        sigsets
    }
    get.gsea.table <- function(gsea, mods, sets, regulations = c("up", "down")) {
        # name conflict: regs = up/down conflicts with the global candidate regulator list names: regs = c(APE2RAB ...)
        cmps <- mods
        regs <- regulations
        gst <- NULL
        n <- length(sets)
        for (cmp in cmps) {
            cat(cmp, "\n")
            for (reg in regs) {
                rnb <- paste(cmp, reg, sep = "_")
                gs <- gsea[[cmp]][[reg]]
                if (is.null(gs)) {
                    m <- data.frame(pval = rep(NA, n),
                                    qval = rep(NA, n),
                                    core = rep(NA, n),
                                    stringsAsFactors = F)
                    rownames(m) <- sets
                } else {
                    m <- gs[, c("NOM.p.val", "FDR.q.val", "core")]
                    colnames(m) <- c("pval", "qval", "core")
                }
                m <- t(m)
                rn <- paste(rnb, rownames(m), sep = "_")
                gst <- rbind(gst, m)
                rownames(gst)[(nrow(gst) - 2):nrow(gst)] <- rn
            }
        }
        gst
     #   as.data.frame(gst, stringsAsFactors = F)
    }
    get.gsea.tables <- function(gseas, nregs = "nreg_all", regs, percs = "perc_100", mods, sets, regulations = c("up", "down")) {
        gsea.tables <- NULL
        for (nreg in nregs)
            for (reg in regs)
                for (perc in percs) {
                    cat(nreg, reg, perc, "\n")
                    gsea <- gseas[[nreg]][[reg]][[perc]]
                    gsea.tables[[nreg]][[reg]][[perc]] <- list()
                    gsea.tables[[nreg]][[reg]][[perc]] <- get.gsea.table(gsea, mods = mods, sets = sets, regulations = regulations)
                }
        gsea.tables
    }
   for (obj in ls()) 
       assign(obj, get(obj), envir = Enrichments)
})
