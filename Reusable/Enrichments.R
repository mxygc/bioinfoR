if (!exists("Enrichments") || !is.environment(Enrichments)) Enrichments <- new.env(parent = emptyenv())

local({

    get.hpg <- function(s, bg, l, n = 3000, seed = 1) {
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
        	phyper(q = q - 1, m = m, n = n, k = k, lower.tail = FALSE)
        }
    	bg <- unique(bg)
    	l <- unique(l)
    	s <- unique(s)
    	i <- paste(intersect(l, s), collapse = ";")
        if (length(l) > 0) {
        	p.real <- get.hypergeom.p(bg = bg, l = l, s = s)
            if (n > 0) {
        	    if (!missing(seed)) { seed = seed; set.seed(seed) }
        	    l.perms <- get.perm.lists(bg = bg, g = l, n = n)
                p.perms <- apply(l.perms, 2, function(l.perm, bg, s) get.hypergeom.p(bg = bg, l = l.perm, s = s), bg = bg, s = s)
            } else 
                p.perms <- rep(NA, n)
        } else {
            p.real <- NA
            p.perms <- rep(NA, n)
            i <- ""
        }
    	list(p = p.real, p.perms = p.perms, i = i)
    }
    
    get.hypergeom <- function(degs, genesets, bg, cl = NULL, n = 3000, seed = 1, cmps = NULL, mc.methods = c("bh", "by", "hlm", "hcbg", "hmml", "bfrn", "btsp"), mc.way = "row") {
    ##degs = list(cmp1 = vector(genes), cmp2 = vector(genes), ...)
    ##genesets = list(set1 = vectors(genes), set2 = vector(geens))
    ##multicorrect.method defines way to do correction: by row, col, or both
        hypergeoms <- NULL
        hpgs <- NULL
        if (is.null(cmps)) cmps <- names(degs)
        ncmp <- length(cmps)
        nset <- length(genesets)
        sets <- names(genesets)
        p.perms <- list()
        hypergeoms$p <- matrix(NA, nrow = ncmp, ncol = nset, dimnames = list(cmps, sets))
        hypergeoms$i <- matrix("", nrow = ncmp, ncol = nset, dimnames = list(cmps, sets))
        get.mc <- function(method, p, cmps, sets, ncmp, nset, p.perms = NULL, way = "both") {
            .methods <- c(bh = "BH", by = "BY", hlm = "holm", hcbg = "hochberg", hmml = "hommel", bfrn = "bonferroni", btsp = "bootstrap")
            m <- .methods[method]
            if (m != "bootstrap") {
                if (way == "both") {
                    mc <- as.vector(p)
                    mc <- p.adjust(mc, method = m)
                    mc <- matrix(mc, nrow = ncmp, ncol = nset, dimnames = list(cmps, sets))
                } else if (way == "row") {
                    mc <- t(apply(p, 1, p.adjust, method = m))
                } else if (way == "col") {
                    mc <- apply(p, 2, p.adjust, method = m)
                }
            } else {
                if (is.null(p.perms)) stop("permutation derived p values don't exist")
                mc <- p
                p.perm <- NULL
                if (way == "both") 
                    p.perm <- unlist(p.perms)
                else if (way == "row")
                    p.perm <- lapply(cmps, function(cmp) unlist(p.perms[[cmp]]))
                else if (way == "col") {
                    p.perm <- lapply(sets, function(set) unlist(lapply(cmps, function(cmp) p.perms[[cmp]][[set]])))
                }
                for (cmp in cmps) 
                    for (set in sets) {
                        if (way == "both")
                            mc[cmp, set] <- mean(p.perm <= p[cmp, set], na.rm = T)
                        else if (way == "row")
                            mc[cmp, set] <- mean(p.perm[[cmp]] <= p[cmp, set], na.rm = T)
                        else if (way == "col")
                            mc[cmp, set] <- mean(p.perm[[set]] <= p[cmp, set], na.rm = T)
                    }
            }
            mc
        }
        rearrange <- function(hpgm, methods, cmps, sets, nset) {
            rns <- c("p", methods, "i")
            nrn <- length(rns)
            hpg <- matrix("", nrow = nrn, ncol = nset, dimnames = list(rns, sets))
            hpgs <- NULL
            for (cmp in cmps) {
                for (rn in rns) {
                    for (set in sets) {
                        hpg[rn, set] <- hpgm[[rn]][cmp, set]
                    }
                }
                hpgs <- c(hpgs, list(hpg))
                names(hpgs)[length(hpgs)] <- cmp
            }
            hpgs
        }
        for (cmp in cmps) {
            deg <- degs[[cmp]];
            cat(cmp, "\n")
            if (!is.null(cl))
                hpgs <- c(hpgs, list(parSapply(cl, genesets, get.hpg, bg = bg, l = deg, n = n, seed = seed, USE.NAMES = T, simplify = F)))
            else 
                hpgs <- c(hpgs, list(sapply(genesets, get.hpg, bg = bg, l = deg, n = n, seed = seed, USE.NAMES = T, simplify = F)))
            names(hpgs)[length(hpgs)] <- cmp
        }
        for (cmp in cmps) {
            for (set in sets) {
                hypergeoms$p[cmp, set] <- hpgs[[cmp]][[set]][["p"]]
                hypergeoms$i[cmp, set] <- hpgs[[cmp]][[set]][["i"]]
                p.perms[[cmp]][[set]] <- hpgs[[cmp]][[set]][["p.perms"]]
            }
        }
        hypergeoms.mc <- sapply(mc.methods, get.mc, p = hypergeoms$p, cmps = cmps, sets = sets, ncmp = ncmp, nset = nset, p.perms = p.perms, way = mc.way, USE.NAMES = T, simplify = F)
        hypergeoms <- c(hypergeoms, hypergeoms.mc)
        hypergeoms <- rearrange(hpgm = hypergeoms, methods = mc.methods, cmps = cmps, sets = sets, nset = nset)
        rownames <- paste(rep(cmps, each = 2 + length(mc.methods)), c("p", mc.methods, "i"), sep = ".")
        hypergeoms.df <- do.call(rbind, hypergeoms)
        rownames(hypergeoms.df) <- rownames
        hypergeoms.df
    }
    
    prepare.expr <- function(expr, annot.map, dir.out) {
    # There are GCT, RES, PCL, TXT 4 kinds of expression formats
    # here just use TXT
        affy <- rownames(expr)
        offsym <- annot.map[affy]
        expr.df <- data.frame(NAME = affy, DESCRIPTION = offsym, expr, stringsAsFactors = F)
        dir.create(dir.out, recursive = T, showWarnings = F)
        file.out <- file.path(dir.out, "expr.txt")
        write.table(file = file.out, expr.df, sep = "\t", row.names = F, quote = F)
    }
    prepare.exprs <- function(exprs, batches, annot.map, dir.out) {
        for (batch in batches) {
            dout <- file.path(dir.out, batch)
            expr <- exprs[[batch]]
            prepare.expr(expr, annot.map, dout)
        }
    }
    prepare.phenotype <- function(expr, dir.out, new.label = TRUE) {
    # CLS format
        samples <- colnames(expr)
        groups <- Labels$get.groups(samples, new.label)
        dir.create(dir.out, recursive = T, showWarnings =  F)
        file <- file.path(dir.out, "pheno.cls")
        nsample <- length(samples)
        ngroup <- length(unique(groups))
        l1 <- paste(nsample, ngroup, 1)
        l2 <- paste("#", paste(unique(groups), collapse = " "))
        l3 <- paste(groups, collapse = " ")
        lines <- paste(l1, l2, l3, sep = "\n")
        cat(file = file, lines)
    }
    prepare.phenotypes <- function(exprs, batches, dir.out, new.label = TRUE) {
        for (batch in batches) {
            expr <- exprs[[batch]]
            dout <- file.path(dir.out, batch)
            prepare.phenotype(expr, dout, new.label)
        }
    }
    prepare.genesets <- function(genesets, sets = NULL, file.out = "genesets.gmx", dir.out) {
    # gmx format
        dir.create(dir.out, recursive = T, showWarnings = F)
        if (is.null(sets)) sets <- colnames(genesets)
        gsets <- genesets[, sets, drop = F]
        gsets <- lapply(gsets, function(x) { x <- x[x != ""]; x <- toupper(x); x <- unique(x); x <- c("", x) })
        gsets <- list2dataframe(gsets, "")
        file <- file.path(dir.out, file.out)
        write.table(gsets, file = file, sep = "\t", quote = F, row.names = F)
    }
    prepare.annot <- function(annot.map, file.out = "annot.chip", dir.out) {
        dir.create(dir.out, recursive = T, showWarnings = F)
        file <- file.path(dir.out, file.out)
        affy <- names(annot.map)
        offsym <- toupper(unname(annot.map))
        annot <- data.frame("Probe Set ID" = affy, 
                            "Gene Symbol" = offsym,
                            "Gene Title" = "na",
                            check.names = F)
        write.table(file = file, annot, row.names = F, sep = "\t", quote = F)
    }
    # for permutation by label (phenotype), the number of replicates in any group should be >= 3
    run.gsea <- function(memsize = 2000, res, cls, gmx, comparison, chip = "ftp.broadinstitute.org://pub/gsea/annotations/Rat230_2.chip", collapse = TRUE, nperm = 10000, permute = "phenotype", rnd_type = "no_balance", seed = "timestamp", set_max = 500, set_min = 5, plot_top_x = 20, norm = "meandiv", save_rnd_lists = FALSE, median = FALSE, num = 100, scoring_scheme = "weighted", make_sets = TRUE, mode = "Max_probe", metric = "Signal2Noise", order = "descending", include_only_symbols = TRUE, sort = "real", zip_report = FALSE, run = T, verbose = F, dir.out) { 
        cmp <- sub("2", "_versus_", comparison)
        cls = paste(cls, "#", cmp, sep = "")
        dir.create(dir.out, recursive = T, showWarnings = F)
        command = paste("java", paste("-Xmx", memsize, "m", sep = ""), 
    					"-cp",
    					"~/Applications/GSEA/gsea2-2.0.10.jar",
    					"xtools.gsea.Gsea", 
    					"-res", res, 
    					"-cls", cls,
    					"-gmx", gmx,
    					"-chip", chip,
    					"-collapse", ifelse(collapse, "true", "false"),
    					"-nperm", nperm,
    					"-permute", permute,
    					"-rpt_label", comparison,
    					"-rnd_seed", seed, 
    					"-save_rnd_lists", ifelse(save_rnd_lists, "true", "false"),
    					"-set_max", set_max, 
    					"-set_min", set_min,
                        "-plot_top_x", plot_top_x, 
                        "-norm", norm, 
                        "-median", ifelse(isTRUE(median), "true", "false"), 
                        "-num", num, 
                        "-scoring_scheme", scoring_scheme,
                        "-make_sets", ifelse(make_sets, "true", "false"),  
                        "-mode", mode, 
                        "-metric", metric, 
                        "-order", order, 
                        "-include_only_symbols", ifelse(include_only_symbols, "true", "false"), 
                        "-sort", sort, 
    					"-zip_report", ifelse(zip_report, "true", "false"), 
    					"-gui", "false", 
    					"-out", dir.out
    					)
    	if (verbose) print(command)
    	if (run) system(command)
    }
    run.gsea.batch <- function(memsize = 2000, res, cls, gmx, comparisons, chip = "ftp.broadinstitute.org://pub/gsea/annotations/Rat230_2.chip", collapse = TRUE, nperm = 10000, permute = "phenotype", rnd_type = "no_balance", seed = "timestamp", set_max = 500, set_min = 5, plot_top_x = 20, norm = "meandiv", save_rnd_lists = FALSE, median = FALSE, num = 100, scoring_scheme = "weighted", make_sets = TRUE, mode = "Max_probe", metric = "Signal2Noise", order = "descending", include_only_symbols = TRUE, sort = "real", zip_report = FALSE, run = T, verbose = F, dir.out) {
        for (comparison in comparisons)
            run.gsea(memsize = memsize, res = res, cls = cls, gmx = gmx, comparison = comparison, chip = chip, collapse = collapse, nperm = nperm, permute = permute, rnd_type = rnd_type, seed = seed, set_max = set_max, set_min = set_min, plot_top_x = plot_top_x, norm = norm, save_rnd_lists = save_rnd_lists, median = median, num = num, scoring_scheme = scoring_scheme, make_sets = make_sets, mode = mode, metric = metric, order = order, include_only_symbols = include_only_symbols, sort = sort, zip_report = zip_report, run = run, verbose = verbose, dir.out = dir.out)
    }
    run.gsea.batches <- function(memsize = 2000, file.res = "expr.txt", file.cls = "pheno.cls", file.chip, file.gmx, comparisons.batch, batches, collapse = TRUE, nperm = 10000, permute = "phenotype", rnd_type = "no_balance", seed = "timestamp", set_max = 500, set_min = 5, plot_top_x = 20, norm = "meandiv", save_rnd_lists = FALSE, median = FALSE, num = 100, scoring_scheme = "weighted", make_sets = TRUE, mode = "Max_probe", metric = "Signal2Noise", order = "descending", include_only_symbols = TRUE, sort = "real", zip_report = FALSE, run = T, verbose = F, dir.in, dir.out) {
        for (batch in batches) {
            bat <- strsplit(batch, "\\.")[[1]][1]
            comparisons <- comparisons.batch[[bat]]
            din <- file.path(dir.in, batch)
            res <- file.path(din, file.res)
            cls <- file.path(din, file.cls)
            dout <- file.path(dir.out, permute, batch)
            dir.create(dout, recursive = T, showWarnings = F)
            run.gsea.batch(memsize = memsize, res = res, cls = cls, gmx = file.gmx, comparisons = comparisons, chip = file.chip, collapse = collapse, nperm = nperm, permute = permute, rnd_type = rnd_type, seed = seed, set_max = set_max, set_min = set_min, plot_top_x = plot_top_x, norm = norm, save_rnd_lists = save_rnd_lists, median = median, num = num, scoring_scheme = scoring_scheme, make_sets = make_sets, mode = mode, metric = metric, order = order, include_only_symbols = include_only_symbols, sort = sort, zip_report = zip_report, run = run, verbose = verbose, dir.out = dout)
        }
    }
    get.gsea <- function(regs = c("up", "down")) { ... }
    
    run.gsea.preranked <- function(jar = "~/Applications/GSEA/gsea2-2.0.10.jar", rnk, gmx, memsize = 1024, collapse = FALSE, nperm = 10000, rpt.label = "my_analysis", rnd.seed = "timestamp", set.max = 500, set.min = 5, save_rnd_lists = FALSE, dir.out, run = TRUE, verbose = FALSE) {
            dir.create(path = dir.out, recursive = T, showWarnings = F)
    	    command = paste(
                        "java", "-cp", jar,
    					paste("-Xmx", memsize, "m", sep = ""), 
    					"xtools.gsea.GseaPreranked",
                        "-rnk", rnk,
    					"-gmx", gmx,
    					"-collapse", ifelse(collapse, "true", "false"), 
    					"-nperm", nperm,
                        "-save_rnd_lists", ifelse(save_rnd_lists, "true", "false"),
    					"-rpt_label", rpt.label, 
    					"-rnd_seed", rnd.seed,
    					"-set_max", set.max,
    					"-set_min", set.min,
    					"-out", dir.out)
    	    if (verbose) print(command)
    	    if (run) system(command)
    }
    run.gsea.preranked.batch <- function(jar = "~/Applications/GSEA/gsea2-2.0.10.jar", degs.preranked.batch, comparisons = NULL, gmx, memsize = 1024, collapse = FALSE, nperm = 10000, set.max = 500, set.min = 5, seed = "timestamp", save_rnd_lists = FALSE, dir.out = "Report/blood/gsea_preranked", run = TRUE, verbose = FALSE) {
        if (is.null(comparisons)) comparisons <- names(degs.rnk)
        for (comparison in comparisons) {
            rnk <- degs.preranked.batch[[comparison]]
            rnk[, 1] <- toupper(rnk[, 1])
            file.out <- paste(comparison, "rnk", sep = ".")
            fout <- file.path(dir.out, fout)
            write.table(rnk, file = fout, row.names = F, col.names = F, sep = "\t", quote = F)
            run.gsea.preranked(jar = jar, rnk = fout, gmx = gmx, memsize = memsize, collapse = collapse, nperm = nperm, rpt.label = comparison, rnd.seed = seed, set.max = set.max, set.min = set.min, dir.out = dir.out, save_rnd_lists = save_rnd_lists, run = run, verbose = verbose)
        }
    }
    run.gsea.preranked.batches <- function(jar = "~/Applications/GSEA/gsea2-2.0.10.jar", degs.preranked.batches, batches, comparisons.batches = NULL, gmx, memsize = 1024, collapse = FALSE, nperm = 10000, set.max = 500, set.min = 5, seed = "timestamp", save_rnd_lists = FALSE, dir.out = "Report/blood/gsea_preranked", run = TRUE, verbose = FALSE) {
        for (batch in batches) {
            degs.preranked.batch <- degs.preranked.batches[[batch]]
            bat <- strsplit(batch, "\\.")[[1]][1]
            comparisons.batch <- comparisons.batches[[bat]]
            run.gsea.preranked.batch(jar = jar, degs.preranked.batch = degs.preranked.batch, comparisons = comparisons.batch, gmx = gmx, memsize = memsize, collapse = collapse, nperm = nperm, set.max = set.max, set.min = set.min, seed = seed, save_rnd_lists = save_rnd_lists, dir.out = dir.out, run = run, verbose = verbose)
        }
    }
    
    get.gsea.batch <- function(dir.in, comparisons, regs = c("up", "down"), genesets, tool = "GseaPreranked") {
        # tool: "Gsea" or "GseaPreranked"
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
        get.stat <- function(cmp, reg, dir, suf, tool = "GseaPreranked") {
            .regs <- c(up = "pos", down = "neg")
            if (tool == "GseaPreranked") {
                f <- paste("gsea", "report", "for", "na", .regs[reg], suf, sep = "_") 
            } else if (tool == "Gsea") {
                .cmps <- strsplit(cmp, "2")[[1]]
                names(.cmps) <- c("up", "down")
                f <- paste("gsea", "report", "for", .cmps[reg], suf, sep = "_")
            }
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
        get.stats <- function(cmps, regs, dir, genesets, tool = "GseaPreranked") {
            suffix <- get.suffix(dir)
            cmps.f <- get.cmp.failed(cmps, suffix)
            stats <- list()
            n <- ncol(genesets)
            for (cmp in cmps) {
                if (!(cmp %in% cmps.f)) {
                    suf <- suffix[cmp]
                    d <- paste(cmp, tool, suf, sep = ".")
                    d <- paste(dir, d, sep = "/")
                    for (reg in regs)
                        stats[[cmp]][[reg]] <- get.stat(cmp = cmp, reg = reg, dir = d, suf = suf, tool = tool)
                } else {
                    for (reg in regs) {
                        stat <- data.frame(
                               "NAME" = colnames(genesets),
                               "ORIGINAL.SIZE" = sapply(genesets, function(x) { x <- x[x != ""]; length(x) }),
                               "SIZE" = rep(0, n), 
                               "ES" = rep(NA, n),
                               "NES" = rep(NA, n),
                               "NOM p-val" = rep(NA, n),
                               "FDR q-val" = rep(NA, n),
                               "FWER p-val" = rep(NA, n),
                               "overlap" = rep("", n), 
                               "core" = rep("", n), 
                               stringsAsFactors = F)
                        rownames(stat) <- colnames(genesets)
                        stats[[cmp]][[reg]] <- stat
                    }
                }
            }
            stats
        }
        stats <- get.stats(cmps = comparisons, regs = regs, dir = dir.in, genesets = genesets, tool = tool)
    }
    get.gsea.batches <- function(dir.in, batches, comparisons.batches, regs = c("up", "down"), genesets, tool = "GseaPreranked") {
        gsea.batches <- list()
        for (batch in batches) {
            bat <- strsplit(batch, "\\.")[[1]][1]
            comparisons <- comparisons.batches[[bat]]
            din <- file.path(dir.in, batch)
            gsea.batches[[batch]] <- get.gsea.batch(dir.in = din, comparisons = comparisons, regs = regs, genesets = genesets, tool = tool)
        }
        gsea.batches
    }
    
    get.gsea.sigsets <- function(gsea, compares, pval = 0.05, qval = NULL, fwer = NULL, n = 2) {
        sigsets <- matrix("", ncol = 2, nrow = length(compares))
        rownames(sigsets) <- compares
        colnames(sigsets) <- c("up", "down")
        get.gs <- function(cmp, gsea, pval = 0.05, qval = NULL, fwer = NULL, n = 2) {
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
            sigsets[cmp, ] <- get.gs(cmp = cmp, gsea = gsea, pval = pval, qval = qval, fwer = fwer, n = n)
    
        }
        sigsets
    }
    
    get.gsea.table.batch <- function(gsea, comparisons, sets, regs = c("up", "down")) {
        gst <- NULL
        n <- length(sets)
        if (is.null(comparisons)) comparisons <- names(gsea)
        for (cmp in comparisons) {
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
    }
    get.gsea.table.batches <- function(gsea.batches, batches = NULL, comparisons.batches = NULL, sets, regs = c("up", "down")) {
        gsea.table.batches <- list()
        if (is.null(batches)) batches <- names(gsea.batches)
        for (batch in batches) {
            bat <- strsplit(batch, "\\.")[[1]][1]
            comparisons <- comparisons.batches[[bat]]
            gsea.batch <- gsea.batches[[batch]]
            gsea.table.batches[[batch]] <- get.gsea.table(gsea.batch, comparisons, sets = sets, regs = regs)
        }
        gsea.table.batches
    }
    assemble.gsea.tables <- function(gsea.tables, headers, labels.rev = c("rev", "norev"), labels.rm = c("norm", "rm"), sets, rows = c("pval", "qval", "core"), regs = c("up", "down")) {
    # headers is a 4-column matrix
    # col1, assay: blood, tissue, joint
    # col2, batch, batch1, batch2, merged
    # col3, comparisons: AP_A2RA_B, ...
    # col4, disp.names: "tissue early batch1"
        gsea.tables.flatted <- list()
        .maps.label <- c(norev = "NOREV", rev = "REV", norm = "HAS_RAB1", rm = "NO_RAB1")
        assemble.gsea.table <- function(i, headers, label.rev, label.rm, sets, rows, regs, gsea.tables) {
            header <- headers[i, ]
            assay <- header[1]
            bat <- header[2]
            comparison <- header[3]
            if (assay == "blood") {
                if (bat == "batch1") {
                    batch <- bat
                } else {
                    batch <- paste(bat, label.rev, label.rm, sep = ".")
                }
            } else if (assay == "tissue") {
                if (bat == "batch1") {
                    batch <- bat
                } else {
                    batch <- paste(bat, label.rm, sep = ".")
                }
            } else if (assay == "joint") {
                batch <- paste(bat, label.rm, sep = ".")
            }
            nreg <- length(regs)
            nrow <- length(rows)
            rownames <- paste(rep(comparison, nreg * nrow), 
                              rep(regs, each = nrow),
                              rows, sep = "_")
            if (is.null(gsea.tables[[assay]][[batch]]) || 
                ! all(rownames %in% rownames(gsea.tables[[assay]][[batch]]))) {
                gsea.table <- matrix("", nrow = length(rownames), ncol = length(sets), dimnames = list(rownames, sets))
            } else {
                gsea.table <- gsea.tables[[assay]][[batch]][rownames, sets]
            }
            data.frame(assay = header[4], comparison = rownames, gsea.table, stringsAsFactors = F, row.names = NULL)
        }
        for (label.rev in labels.rev) 
            for (label.rm in labels.rm) {
                label <- paste(.maps.label[label.rev], .maps.label[label.rm], sep = "_")
                gsea.table.flatted <- lapply(1:nrow(headers), assemble.gsea.table, headers = headers, label.rev = label.rev, label.rm = label.rm, sets = sets, rows = rows, regs = regs, gsea.tables)
                gsea.table.flatted <- do.call(rbind, gsea.table.flatted)
                gsea.tables.flatted[[label]] <- gsea.table.flatted
            }
        gsea.tables.flatted
    }

   for (obj in ls()) 
       assign(obj, get(obj), envir = Enrichments)
})
