#' RNAseq differential analyses using LIMMA, edgeR, DESeq, baySeq

if (!exists("DiffAnal") || !is.environment(DiffAnal)) DiffAnal <- new.env(parent = emptyenv())

local( {
    get.limma.stat <- function(exprs, factors, comparison, adjust = "BH", number = Inf, file = NULL, save = FALSE) {
        require(limma)
        comps <- rev(strsplit(comparison, "2")[[1]])
        comp1 <- comps[1]; comp2 <- comps[2]
        idx1 <- which(factors == comp1); idx2 <- which(factors == comp2)
        expr <- exprs[, c(idx1, idx2)]
        design <- model.matrix(~ 0 + factor(rep(c(1, 2), c(length(idx1), length(idx2)))))
        colnames(design) <- comps
        fit <- lmFit(expr, design)
        contrasts <- paste(comp2, comp1, sep = "-")
        contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit3 <- eBayes(fit2)
        coef <- which(colnames(fit2$contrasts) == contrasts)
        limma <- topTable(fit3, coef = coef, adjust = adjust, number = number)
        id <- rownames(limma)
        reg <- c("down", "nochange", "up")[sign(limma[, "logFC"]) + 2]
        stat <- data.frame(expr[id, ], limma, reg, stringsAsFactors = F)
        if (save && !is.null(file)) {
            write.csv(stat, file = file)
        }
        stat
    }
    get.limma.deg <- function(stat, avgexpr = NULL, pval = 0.05, qval = NULL, fc = 2, file = NULL, save = FALSE) {
        idx <- rep(T, nrow(stat))
        ids <- rownames(stat)
        if (is.numeric(avgexpr) && avgexpr > 0) idx <- idx & stat[, "AveExpr"] >= avgexpr
        if (is.numeric(pval) && pval >= 0 && pval <= 1) idx <- idx & stat[, "P.Value"] < pval
        if (is.numeric(qval) && qval >= 0 && qval <= 1) idx <- idx & stat[, "adj.P.Val"] < qval
        if (is.numeric(fc) && fc != 0) idx <- idx & abs(stat[, "logFC"]) >= log2(abs(fc))
        up <- ids[idx & stat[, "reg"] == "up"]
        down <- ids[idx & stat[, "reg"] == "down"]
        both <- ids[idx]
        deg <- list(up = up, down = down, both = both)
        if (save && !is.null(file)) {
            df.deg <- list2dataframe(deg, "")
            tryCatch(write.csv(df.deg, file = file, row.names = F),
                     error = function(e) {
                         browser() }, finally = cat("ok\n")
                    )
        }
        deg
    }
    # ============================================================
    # CPM after scale normalization + voom transformation
    # + limma
    # input should be not log
    # ============================================================
    get.diffs.limma <- function(counts, cpms, global = FALSE, nrms, modes, levels, annots, aligns, groups, compares, batches = NULL, delta = 1, rs = TRUE, fc.geom = FALSE, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, lm.adjust = "BH", lm.fc = NULL, lm.avgexpr = NULL, lm.pval = 0.05, lm.qval = NULL, lm.bval = NULL, rm.be = FALSE, ftds = NULL, stats = NULL, update = FALSE, multicore = FALSE, verbose = FALSE) {
        require("limma")
        require("edgeR")
        get.foldchange <- function(e1, e2, logged = FALSE, delta = 1, geom = FALSE) {
            if (geom) {
                if (logged) {
                    fc <- rowMeans(e2, na.rm = T) - rowMeans(e1, na.rm = T)
                    avgexpr <- rowMeans(cbind(2 ^ e1, 2 ^ e2), na.rm = T) - delta
                }
                else {
                    fc <- rowMeans(log2(e2 + delta), na.rm = T) - rowMeans(log2(e1 + delta), na.rm = T)
                    avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
                }
            } else {
               if (logged) {
                    e1 <- 2 ^ e1 - delta
                    e2 <- 2 ^ e2 - delta
                }
                fc <- (rowMeans(e2, na.rm = T) + delta) / (rowMeans(e1, na.rm = T) + delta)
                avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
            }
            data.frame(fc = fc, avgexpr = avgexpr)
        }
        get.ranksum <- function(e1, e2) {
            get.wilcox.ranksum <- function(x, y, n1 = NULL, shift = NULL) { # 0.009
                if (is.null(n1)) n1 <- length(x)
                if (is.null(shift)) shift <- n1 * (n1 + 1) / 2
                z <- c(x, y)
                r <- rank(z)
                sum(r[1:n1]) - shift # make it behave as R
            }
            e <- cbind(e1, e2)
            n1 <- ncol(e1)
            n2 <- ncol(e2)
            n <- n1 + n2
            rs.null <- (n1 * n2) / 2
            shift <- (n1 * (n1 + 1)) / 2
            get.rs <- function(e, n1, n) get.wilcox.ranksum(e[1:n1], e[(n1 + 1):n], n1, shift) 
            ranksum <- apply(e, 1, get.rs, n1 = n1, n = n)
            dranksum <- rs.null - ranksum
            data.frame(ranksum = ranksum, dranksum = dranksum)
        }
        get.limma <- function(count, cpm, nrm, group, compare, batch = NULL, verbose = FALSE, delta = 1, adjust = "BH", rm.be = FALSE) {
            grps <- rev(strsplit(compare, "2")[[1]])
            grpns <- sapply(grps, function(g) sum(group == g))
            if ((!is.null(batch)) & rm.be) {
                design <- model.matrix(~0+factor(rep(c(1, 2), grpns))+batch)
                colnames(design) <- c(grps, "BATCH")
            } else {
                design <- model.matrix(~0+factor(rep(c(1, 2), grpns))) 
                colnames(design) <- grps
            }
            cds <- edgeR::DGEList(counts = count, group = group)
            cds <- edgeR::calcNormFactors(cds, method = c(tmm = "TMM", rle = "RLE", uq = "upperquartile")[nrm])
            v <- voom(counts = cds, design = design, lib.size = NULL, plot = T)
            fit <- lmFit(v, design)
            compare <- sub(pattern = "2", x = compare, replacement = "-")
            contrast.matrix <- makeContrasts(contrasts = compare, levels = design)
            fit2 <- contrasts.fit(fit, contrast.matrix)
            fit3 <- eBayes(fit2)
            coef <- which(colnames(fit2$contrasts) == sub(pattern = "2", replacement = "-", x = compare))
            limma <- topTable(fit3, coef = coef, adjust = adjust, number = Inf)
            if (verbose) {
                elog <- v$E
                colnames(elog) <- paste("log", colnames(elog), sep = "_")
                id <- rownames(limma)
                elog <- elog[id, ]
                limma <- data.frame(elog, limma)
            }
            limma
        }
        get.stat <- function(compare, count, cpm, batches = NULL, elog = NULL, nrm, fit = NULL, global = FALSE, groups, delta = 1, rs = TRUE, fc.geom = FALSE, lm.adjust = "BH", rm.be = FALSE, ftd = NULL, verbose = FALSE) {
            cat(compare, "\n")
            grps <- rev(strsplit(compare, "2")[[1]])
            cols.1 <- which(groups == grps[1])
            cols.2 <- which(groups == grps[2])
            group <- groups[c(cols.1, cols.2)]
            batch <- batches[c(cols.1, cols.2)]
            id <- rownames(count)
            if (global & (!is.null(fit))) {
                f <- NULL
                left <- id[!id %in% f]
                coef <- which(colnames(fit$contrasts) == sub(pattern = "2", replacement = "-", x = compare))
                limma <- topTable(fit, coef = coef, adjust = lm.adjust, number = Inf)
                if (verbose) {
                    id <- rownames(limma)
                    elog <- elog[id, c(cols.1, cols.2)]
                    colnames(elog) <- paste("log", colnames(elog), sep = "_")
                    limma <- data.frame(elog, limma)
                }
            } else {
                if (!is.null(ftd)) f <- ftd[[compare]]
                else f <- NULL
                left <- id[!id %in% f]
                cnt <- count[left, c(cols.1, cols.2)]
                cm <- cpm[left, c(cols.1, cols.2)]
                limma <- get.limma(count = cnt, cpm = cm, nrm = nrm, group = group, compare = compare, batch = batch, delta = delta, adjust = lm.adjust, rm.be = rm.be, verbose = verbose)
            }
            cm1 <- cpm[left, cols.1]
            cm2 <- cpm[left, cols.2]
            foldchange <- get.foldchange(cm1, cm2, delta = delta, logged = F, geom = fc.geom)
            if (isTRUE(rs)) ranksum <- get.ranksum(cm1, cm2)
            reg <- c("down", "nc", "up")[sign(log2(foldchange[, "fc"])) + 2]
            names(reg) <- rownames(foldchange)
            id <- rownames(limma)
            stat <- data.frame(foldchange[id, ])
            if (isTRUE(rs)) stat <- data.frame(stat, ranksum[id, ])
            stat <- data.frame(stat, limma[id, ])
            stat <- data.frame(stat, reg = reg[id])
        }
        get.stats <- function(counts, cpms, fit = NULL, global = FALSE, nrms, modes, levels, annots, aligns, groups, compares, batches = NULL, rs = TRUE, delta = 1, fc.geom = FALSE, lm.adjust = "BH", rm.be = FALSE, ftds = NULL, cl = NULL, verbose = F) {
            stats <- NULL
            for (nrm in nrms)
                for (mode in modes) 
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns) {
                                cat(nrm, mode, level, annot, align, "\n")
                                count <- counts[[mode]][[level]][[annot]][[align]]
                                cpm <- cpms[[nrm]][[mode]][[level]][[annot]][[align]]
                                if (global)  {
                                    if (!is.null(ftds)) ftd <- ftds[["global"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    else ftd <- NULL
                                    ids <- rownames(count)
                                    left <- ids[!ids %in% ftd]
                                    cnt <- count[left, ]
                                    cm <- cpm[left, ]
                                    cds <- edgeR::DGEList(counts = cnt, group = groups)
                                    cds <- edgeR::calcNormFactors(cds, method = c(tmm = "TMM", rle = "RLE", uq = "upperquartile")[nrm]) 
                                    grps <- unique(groups)
                                    grpns <- sapply(grps, function(g) sum(groups == g))
                                    if ((!is.null(batches)) & rm.be) {
                                        design <- model.matrix(~0+factor(rep(1:length(grps), grpns))+batches)
                                        colnames(design) <- c(as.character(grps), "BATCH")
                                    } else {
                                        design <- model.matrix(~0+factor(rep(1:length(grps), grpns)))
                                        colnames(design) <- as.character(grps)
                                    }
                                    v <- voom(counts = cds, design = design, lib.size = NULL, plot = T)
                                    elog <- v$E
                                    fit <- lmFit(v, design)
                                    contrast.matrix <- makeContrasts(contrasts = sub(pattern = "2", x = compares, replacement = "-"), levels = design)
                                    fit <- contrasts.fit(fit, contrast.matrix)
                                    fit <- eBayes(fit)
                                } else {
                                    if (!is.null(ftds)) ftd <- ftds[["local"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    else ftd <- NULL
                                    fit <- NULL
                                    elog <- NULL
                                }
                                stats[[nrm]][[mode]][[level]][[annot]][[align]] <- list()
                                if (!is.null(cl)) {
                                    stats[[nrm]][[mode]][[level]][[annot]][[align]] <- parSapply(cl = cl, compares = compares, batches = batches, FUN = get.stat, count = count, cpm = cpm, elog = elog, nrm = nrm, fit = fit, global = global, groups = groups, rs = rs, delta = delta, fc.geom = fc.geom, lm.adjust = lm.adjust, rm.be = rm.be, ftd = ftd, verbose = verbose, USE.NAMES = T, simplify = F)
                                } else {
                                    stats[[nrm]][[mode]][[level]][[annot]][[align]] <- sapply(compares, FUN = get.stat, count = count, cpm = cpm, batches = batches, elog = elog, nrm = nrm, fit = fit, global = global, groups = groups, rs = rs, delta = delta, fc.geom = fc.geom, lm.adjust = lm.adjust, rm.be = rm.be, ftd = ftd, verbose = verbose, USE.NAMES = T, simplify = F)
                                }
                            }
            stats
        }
        get.deuid <- function(stat, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, lm.fc = NULL, lm.avgexpr = NULL, lm.pval = 0.05, lm.qval = NULL, lm.bval = NULL) {
            flag <- rep(T, nrow(stat))
            if (!is.null(fc.fc)) flag <- flag & (abs(log2(stat[, "fc"])) >= log2(fc.fc))
            if (!is.null(fc.avgexpr)) flag <- flag & (stat[, "avgexpr"] >= fc.avgexpr)
            if (!is.null(rs.dranksum)) flag <- flag & abs(stat[, "dranksum"]) >= rs.dranksum
            if (!is.null(lm.fc)) flag <- flag & (abs(stat[, "logFC"]) >= log2(lm.fc))
            if (!is.null(lm.avgexpr)) flag <- flag & (stat[, "AveExpr"] >= log2(lm.avgexpr))
            if (!is.null(lm.pval)) flag <- flag & (stat[, "P.Value"] < lm.pval)
            if (!is.null(lm.qval)) flag <- flag & (stat[, "adj.P.Val"] < lm.qval)
            if (!is.null(lm.bval)) flag <- flag & (stat[, "B"] >= lm.bval)
            flag[is.na(flag)] <- F
            reg <- stat[flag, "reg"]
            deuid <- rownames(stat)[flag]
            up <- deuid[reg == "up"]
            down <- deuid[reg == "down"]
            list(up = up, down = down, both = deuid)
        }
        get.deuids <- function(stats, nrms, modes, levels, annots, aligns, compares, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, lm.fc = NULL, lm.avgexpr = NULL, lm.pval = 0.05, lm.qval = NULL, lm.bval = NULL) {
            deuids <- NULL
            for (nrm in nrms)
                for (mode in modes)
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns)
                                for (compare in compares) {
                                   stat <- stats[[nrm]][[mode]][[level]][[annot]][[align]][[compare]]
                                   deuid <- get.deuid(stat, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, lm.fc = lm.fc, lm.avgexpr = lm.avgexpr, lm.pval = lm.pval, lm.qval = lm.qval, lm.bval = lm.bval)
                                   deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- list()
                                   deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- deuid
                                }
            deuids
        }
        if (!update | is.null(stats)) {
            if (multicore & (length(compares) > 1)) {
                require("parallel")
                cl <- makeCluster(length(compares))
            } else cl <- NULL
            stats <- get.stats(counts = counts, cpms = cpms, global = global, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, groups = groups, compares = compares, batches = batches, rs = rs, delta = delta, fc.geom = fc.geom, lm.adjust = lm.adjust, rm.be = rm.be, ftds = ftds, cl = cl, verbose = verbose)
            if (multicore & (length(compares) > 1)) stopCluster(cl)
        }
        deuids <- get.deuids(stats = stats, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, compares = compares, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, lm.fc = lm.fc, lm.avgexpr = lm.avgexpr, lm.pval = lm.pval, lm.qval = lm.qval, lm.bval = lm.bval)
        diffs <- list(stats = stats, deuids = deuids)
        new.diffs <- NULL
        for (nrm in nrms)
            for (mode in modes) {
                new.diffs[[nrm]][[mode]] <- list(stats = diffs[["stats"]][[nrm]][[mode]], deuids = diffs[["deuids"]][[nrm]][[mode]])
            }
        return(new.diffs)
    }
    
    get.diffs.edger <- function(counts, cpms, global = FALSE, nrms, modes, levels, annots, aligns, groups, compares, logged = FALSE, delta = 1, rs = FALSE, fc.geom = FALSE, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, edgr.fc = NULL, edgr.avgcpm = NULL, edgr.pval = 0.05, edgr.qval = NULL, ftds = NULL, stats = NULL, update = FALSE, multicore = FALSE) {
        get.foldchange <- function(e1, e2, logged = FALSE, delta = 1, geom = FALSE) {
            if (geom) {
                if (logged) {
                    fc <- rowMeans(e2, na.rm = T) - rowMeans(e1, na.rm = T)
                    avgexpr <- rowMeans(cbind(2 ^ e1, 2 ^ e2), na.rm = T) - delta
                }
                else {
                    fc <- rowMeans(log2(e2 + delta), na.rm = T) - rowMeans(log2(e1 + delta), na.rm = T)
                    avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
                }
            } else {
               if (logged) {
                    e1 <- 2 ^ e1 - delta
                    e2 <- 2 ^ e2 - delta
                }
                fc <- (rowMeans(e2, na.rm = T) + delta) / (rowMeans(e1, na.rm = T) + delta)
                avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
            }
            data.frame(fc = fc, avgexpr = avgexpr)
        }
        get.ranksum <- function(e1, e2) {
            get.wilcox.ranksum <- function(x, y, n1 = NULL, shift = NULL) { # 0.009
                if (is.null(n1)) n1 <- length(x)
                if (is.null(shift)) shift <- n1 * (n1 + 1) / 2
                z <- c(x, y)
                r <- rank(z)
                sum(r[1:n1]) - shift # make it behave as R
            }        
            e <- cbind(e1, e2)
            n1 <- ncol(e1)
            n2 <- ncol(e2)
            n <- n1 + n2
            rs.null <- (n1 * n2) / 2
            shift <- (n1 * (n1 + 1)) / 2
            get.rs <- function(e, n1, n) get.wilcox.ranksum(e[1:n1], e[(n1 + 1):n], n1, shift) 
            ranksum <- apply(e, 1, get.rs, n1 = n1, n = n)
            dranksum <- rs.null - ranksum
            data.frame(ranksum = ranksum, dranksum = dranksum)
        }
        get.edger <- function(e, group, compare, nrm) {
            grps <- rev(strsplit(compare, "-")[[1]])
            cds <- edgeR::DGEList(counts = e, group = group)
            cds <- edgeR::calcNormFactors(cds, method = c(tmm = "TMM", rle = "RLE", uq = "upperquartile")[nrm])
            cds <- edgeR::estimateCommonDisp(cds)
            cds <- edgeR::estimateTagwiseDisp(cds)
            disp <- cds$tagwise.dispersion
            names(disp) <- rownames(cds$counts)
            et <- edgeR::exactTest(cds, pair = grps)
            tt <- edgeR::topTags(et, n = nrow(et))
            tt.table <- tt$table
            rn <- rownames(tt.table)
            edger <- data.frame(tagwiseDisp = disp[rn], tt.table[rn, ])
            edger
        }
        get.stat <- function(compare, cds = NULL, count, cpm, global = FALSE, nrm = NULL, groups, logged = FALSE, delta = 1, rs = FALSE, fc.geom = FALSE, ftd = NULL) {
            require("edgeR")
            cat(compare, "\n")
            grps <- rev(strsplit(compare, "-")[[1]])
            cols.1 <- which(groups == grps[1])
            cols.2 <- which(groups == grps[2])
            group <- groups[c(cols.1, cols.2)]
            if (global & (!is.null(cds))) {
                left <- rownames(cds$counts)
                disp <- cds$tagwise.dispersion
                names(disp) <- left
                et <- edgeR::exactTest(cds, pair = grps)
                tt <- edgeR::topTags(et, n = nrow(et))
                tt.table <- tt$table
                rn <- rownames(tt.table)
                edger <- data.frame(tagwiseDisp = disp[rn], tt.table[rn, ])
            } else {
                if (!is.null(ftd)) {
                    f <- ftd[[compare]]
                } else f <- NULL
                id <- rownames(count)
                left <- id[!id %in% f]
                c <- count[left, c(cols.1, cols.2)]
                edger <- get.edger(c, group, compare, nrm)
            }
            cm <- cpm[left, c(cols.1, cols.2)]
            cm1 <- cpm[left, cols.1]
            cm2 <- cpm[left, cols.2]
            foldchange <- get.foldchange(cm1, cm2, delta = delta, logged = logged, geom = fc.geom)
            if (isTRUE(rs)) ranksum <- get.ranksum(cm1, cm2)
            reg <- c("down", "nc", "up")[sign(log2(foldchange[, "fc"])) + 2]
            names(reg) <- rownames(foldchange)
            id <- rownames(edger)
            stat <- data.frame(foldchange[id, ])
            if (isTRUE(rs)) stat <- data.frame(stat, ranksum[id, ])
            stat <- data.frame(stat, edger[id, ])
            stat <- data.frame(stat, reg = reg[id])
        }
        get.stats <- function(counts, cpms, global = FALSE, nrms, modes, levels, annots, aligns, groups, compares, rs = FALSE, logged = FALSE, delta = 1, fc.geom = FALSE, ftds = NULL, cl = NULL) {
            stats <- NULL
            for (nrm in nrms) 
                for (mode in modes)
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns) {
                                cat(nrm, mode, level, annot, align, "\n")
                                cpm <- cpms[[nrm]][[mode]][[level]][[annot]][[align]]
                                count <- counts[[mode]][[level]][[annot]][[align]]
                                ftd <- NULL
                                if (global) {
                                    if (!is.null(ftds)) ftd <- ftds[["global"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    id <- rownames(count)
                                    left <- id[!id %in% ftd]
                                    count <- count[left, ]
                                    cmp <- cpm[left, ]
                                    cds <- edgeR::DGEList(counts = count, group = groups)
                                    cds <- edgeR::calcNormFactors(cds, method = c(tmm = "TMM", rle = "RLE", uq = "upperquartile")[nrm])
                                    cds <- edgeR::estimateCommonDisp(cds)
                                    cds <- edgeR::estimateTagwiseDisp(cds)
                                } else {
                                    if (!is.null(ftds)) ftd <- ftds[["local"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    cds <- NULL
                                }
                                stats[[nrm]][[mode]][[level]][[annot]][[align]] <- list()
                                if (is.null(cl)) {
                                    stats[[nrm]][[mode]][[level]][[annot]][[align]] <- sapply(compares, get.stat, cds = cds, count = count, cpm = cpm, global = global, nrm = nrm, groups = groups, rs = rs, logged = logged, delta = delta, fc.geom = fc.geom, ftd = ftd, USE.NAMES = T, simplify = F)
                                } else {
                                    stats[[nrm]][[mode]][[level]][[annot]][[align]] <- parSapply(cl, compares, get.stat, cds = cds, count = count, cpm = cpm, global = global, nrm = nrm, groups = groups, rs = rs, logged = logged, delta = delta, fc.geom = fc.geom, ftd = ftd, USE.NAMES = T, simplify = F)
                                }
                            }
            stats
        }
        get.deuid <- function(stat, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, edgr.fc = NULL, edgr.avgcpm = NULL, edgr.pval = 0.05, edgr.qval = NULL) {
            flag <- rep(T, nrow(stat))
            if (!is.null(fc.fc)) flag <- flag & (abs(log2(stat[, "fc"])) >= log2(fc.fc))
            if (!is.null(fc.avgexpr)) flag <- flag & (stat[, "avgexpr"] >= fc.avgexpr)
            if (!is.null(rs.dranksum)) flag <- flag & abs(stat[, "dranksum"]) >= rs.dranksum
            if (!is.null(edgr.fc)) flag <- flag & (stat[, "FoldChange"] >= edgr.fc)
            if (!is.null(edgr.avgcpm)) flag <- flag & (stat[, "logCPM"] > log2(edgr.avgcpm))
            if (!is.null(edgr.pval)) flag <- flag & (stat[, "PValue"] < edgr.pval)
            if (!is.null(edgr.qval)) flag <- flag & (stat[, "FDR"] < edgr.qval)
            flag[is.na(flag)] <- F
            reg <- stat[flag, "reg"]
            deuid <- rownames(stat)[flag]
            up <- deuid[reg == "up"]
            down <- deuid[reg == "down"]
            list(up = up, down = down, both = deuid)
        }
        get.deuids <- function(stats, nrms, modes, levels, annots, aligns, compares, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, edgr.fc = NULL, edgr.avgcpm = NULL, edgr.pval = 0.05, edgr.qval = NULL) {
            deuids <- NULL
            for (nrm in nrms)
                for (mode in modes)
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns) 
                                for (compare in compares) {
                                    stat <- stats[[nrm]][[mode]][[level]][[annot]][[align]][[compare]]
                                    deuid <- get.deuid(stat, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, edgr.fc = edgr.fc, edgr.avgcpm = edgr.avgcpm, edgr.pval = edgr.pval, edgr.qval = edgr.qval)
                                    deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- list()
                                    deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- deuid
                                }
            deuids
        }
        if (!update | is.null(stats)) {
            if (multicore & (length(compares) > 1)) {
                require("parallel")
                cl <- makeCluster(length(compares))
            } else cl <- NULL
            stats <- get.stats(counts = counts, cpms = cpms, global = global, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, groups = groups, compares = compares, rs = rs, logged = logged, delta = delta, fc.geom = fc.geom, ftds = ftds, cl = cl)
            if (multicore & (length(compares) > 1)) stopCluster(cl)
        } else {
            stats.1 <- NULL
            for (nrm in nrms)
                for (mode in modes) {
                    stats.1[[nrm]][[mode]] <- stats[[nrm]][[mode]][["stats"]]
                }
            stats <- stats.1
        }
        deuids <- get.deuids(stats = stats, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, compares = compares, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, edgr.fc = edgr.fc, edgr.avgcpm = edgr.avgcpm, edgr.pval = edgr.pval, edgr.qval = edgr.qval)
        diffs <- list(stats = stats, deuids = deuids)
        new.diffs <- NULL
        for (nrm in nrms)
            for (mode in modes) {
                new.diffs[[nrm]][[mode]] <- list(stats = diffs[["stats"]][[nrm]][[mode]], deuids = diffs[["deuids"]][[nrm]][[mode]])
            }
        new.diffs
    }
    # ============================================================
    # DESeq
    # note: only RLE is available for DESeq
    # ============================================================
    get.diffs.deseq <- function(counts, cpms, global = FALSE, nrms = "rle", modes = c("union", "internonempty", "interstrict"), levels, annots, aligns, groups, compares, logged = FALSE, delta = 1, rs = FALSE, fc.geom = FALSE, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, dsq.fc = NULL, dsq.avgcpm = NULL, dsq.pval = 0.05, dsq.qval = NULL, ftds = NULL, stats = NULL, update = FALSE) {
        require("DESeq")
        get.foldchange <- function(e1, e2, logged = FALSE, delta = 1, geom = FALSE) {
            if (geom) {
                if (logged) {
                    fc <- rowMeans(e2, na.rm = T) - rowMeans(e1, na.rm = T)
                    avgexpr <- rowMeans(cbind(2 ^ e1, 2 ^ e2), na.rm = T) - delta
                }
                else {
                    fc <- rowMeans(log2(e2 + delta), na.rm = T) - rowMeans(log2(e1 + delta), na.rm = T)
                    avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
                }
            } else {
               if (logged) {
                    e1 <- 2 ^ e1 - delta
                    e2 <- 2 ^ e2 - delta
                }
                fc <- (rowMeans(e2, na.rm = T) + delta) / (rowMeans(e1, na.rm = T) + delta)
                avgexpr <- rowMeans(cbind(e1, e2), na.rm = T)
            }
            data.frame(fc = fc, avgexpr = avgexpr)
        }
        get.ranksum <- function(e1, e2) {
            e <- cbind(e1, e2)
            n1 <- ncol(e1)
            n2 <- ncol(e2)
            n <- n1 + n2
            rs.null <- (n1 * n2) / 2
            shift <- (n1 * (n1 + 1)) / 2
            get.rs <- function(e, n1, n) get.wilcox.ranksum(e[1:n1], e[(n1 + 1):n], n1, shift) 
            ranksum <- apply(e, 1, get.rs, n1 = n1, n = n)
            dranksum <- rs.null - ranksum
            data.frame(ranksum = ranksum, dranksum = dranksum)
        }
        get.deseq <- function(count, group, compare) {
            grps <- rev(strsplit(compare, "-")[[1]])
            ids <- rownames(count)
            cds <- DESeq::newCountDataSet(countData = count, conditions = group)
            cds <- DESeq::estimateSizeFactors(cds)
            cds <- DESeq::estimateDispersions(cds)
            fit <- fitInfo(cds)
            disp <- data.frame(perGeneDispEsts = fit[["perGeneDispEsts"]], fittedDispEsts = fit[["fittedDispEsts"]], pooledDisp = fData(cds)[["disp_pooled"]])
            rownames(disp) <- ids
            nbt <- nbinomTest(cds = cds, condA = grps[1], condB = grps[2])
            nbt <- nbt[order(nbt[, "pval"]), ]
            ids <- nbt[, "id"]
            rownames(nbt) <- ids
            nbt <- nbt[, -1]
            deseq <- data.frame(disp[ids, ], nbt[ids, ])
        }
        get.stat <- function(compare, cds = NULL, count, cpm, global = FALSE, groups, logged = FALSE, delta = 1, rs = FALSE, fc.geom = FALSE, ftd = NULL) {
            cat(compare, "\n")
            grps <- rev(strsplit(compare, "-")[[1]])
            cols.1 <- which(groups == grps[1])
            cols.2 <- which(groups == grps[2])
            group <- groups[c(cols.1, cols.2)]
            if (global & (!is.null(cds))) {
                left <- rownames(cds@assayData$counts)
                fit <- fitInfo(cds)
                disp <- data.frame(perGeneDispEsts = fit[["perGeneDispEsts"]], fittedDispEsts = fit[["fittedDispEsts"]], pooledDisp = fData(cds)[["disp_pooled"]])
                rownames(disp) <- left
                nbt <- nbinomTest(cds = cds, condA = grps[1], condB = grps[2])
                nbt <- nbt[order(nbt[, "pval"]), ]
                ids <- nbt[, "id"]
                rownames(nbt) <- ids
                nbt <- nbt[, -1]
                deseq <- data.frame(disp[ids, ], nbt[ids, ])
            } else {
                if (!is.null(ftd)) {
                    f <- ftd[[compare]]
                } else f <- NULL
                id <- rownames(count)
                left <- id[!id %in% f]
                cnt <- count[left, c(cols.1, cols.2)]
                deseq <- get.deseq(count = cnt, group = group, compare = compare)
            }
            cm <- cpm[left, c(cols.1, cols.2)]
            cm1 <- cpm[left, cols.1]
            cm2 <- cpm[left, cols.2]
            foldchange <- get.foldchange(cm1, cm2, delta = delta, logged = logged, geom = fc.geom)
            if (isTRUE(rs)) ranksum <- get.ranksum(cm1, cm2)
            # here the regulation direction is by FC
            reg <- c("down", "nc", "up")[sign(log2(foldchange[, "fc"])) + 2]
            names(reg) <- rownames(foldchange)
            id <- rownames(deseq)
            stat <- data.frame(foldchange[id, ])
            if (isTRUE(rs)) stat <- data.frame(stat, ranksum[id, ])
            stat <- data.frame(stat, deseq[id, ])
            stat <- data.frame(stat, reg = reg[id])
        }
        get.stats <- function(counts, cpms, global = FALSE, nrms, modes, levels, annots, aligns, groups, compares, rs = FALSE, logged = FALSE, delta = 1, fc.geom = FALSE, ftds = NULL) {
            stats <- NULL
            for (nrm in nrms) 
                for (mode in modes)
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns) {
                                cat(nrm, mode, level, annot, align, "\n")
                                cpm <- cpms[[nrm]][[mode]][[level]][[annot]][[align]]
                                count <- counts[[mode]][[level]][[annot]][[align]]
                                ftd <- NULL
                                if (global) {
                                    if (!is.null(ftds)) ftd <- ftds[["global"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    id <- rownames(count)
                                    left <- id[!id %in% ftd]
                                    count <- count[left, ]
                                    cpm<- cpm[left, ]
                                    cds <- DESeq::newCountDataSet(countData = count, conditions = groups)
                                    cds <- DESeq::estimateSizeFactors(cds)
                                    cds <- DESeq::estimateDispersions(cds) 
                                    count <- NULL
                                } else {
                                    if (!is.null(ftds)) ftd <- ftds[["local"]][[nrm]][[mode]][[level]][[annot]][[align]]
                                    cds <- NULL
                                }
                                stats[[nrm]][[mode]][[level]][[annot]][[align]] <- list()
                                stats[[nrm]][[mode]][[level]][[annot]][[align]] <- sapply(compares, get.stat, cds = cds, count = count, cpm = cpm, global = global, groups = groups, rs = rs, logged = logged, delta = delta, fc.geom = fc.geom, ftd = ftd, USE.NAMES = T, simplify = F)
                            }
            stats
        }
        get.deuid <- function(stat, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, dsq.fc = NULL, dsq.avgcpm = NULL, dsq.pval = 0.05, dsq.qval = NULL) {
            flag <- rep(T, nrow(stat))
            if (!is.null(fc.fc)) flag <- flag & (abs(log2(stat[, "fc"])) >= log2(fc.fc))
            if (!is.null(fc.avgexpr)) flag <- flag & (stat[, "avgexpr"] >= fc.avgexpr)
            if (!is.null(rs.dranksum)) flag <- flag & abs(stat[, "dranksum"]) >= rs.dranksum
            if (!is.null(dsq.fc)) flag <- flag & (abs(stat[, "log2FoldChange"]) >= log2(dsq.fc))
            if (!is.null(dsq.avgcpm)) flag <- flag & (stat[, "baseMean"] > dsq.avgcpm)
            if (!is.null(dsq.pval)) flag <- flag & (stat[, "pval"] < dsq.pval)
            if (!is.null(dsq.qval)) flag <- flag & (stat[, "padj"] < dsq.qval)
            flag[is.na(flag)] <- F
            reg <- stat[flag, "reg"]
            deuid <- rownames(stat)[flag]
            up <- deuid[reg == "up"]
            down <- deuid[reg == "down"]
            list(up = up, down = down, both = deuid)
        }
        get.deuids <- function(stats, nrms, modes, levels, annots, aligns, compares, fc.fc = 3, fc.avgexpr = 1, rs.dranksum = 2, dsq.fc = NULL, dsq.avgcpm = NULL, dsq.pval = 0.05, dsq.qval = NULL) {
            deuids <- NULL
            for (nrm in nrms)
                for (mode in modes)
                    for (level in levels)
                        for (annot in annots)
                            for (align in aligns) 
                                for (compare in compares) {
                                    stat <- stats[[nrm]][[mode]][[level]][[annot]][[align]][[compare]]
                                    deuid <- get.deuid(stat, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, dsq.fc = dsq.fc, dsq.avgcpm = dsq.avgcpm, dsq.pval = dsq.pval, dsq.qval = dsq.qval)
                                    deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- list()
                                    deuids[[nrm]][[mode]][[level]][[annot]][[align]][[compare]] <- deuid
                                }
            deuids
        }
        if (!update | is.null(stats)) {
            stats <- get.stats(counts = counts, cpms = cpms, global = global, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, groups = groups, compares = compares, rs = rs, logged = logged, delta = delta, fc.geom = fc.geom, ftds = ftds)
        } else {
            stats.1 <- NULL
            for (nrm in nrms)
                for (mode in modes) {
                    stats.1[[nrm]][[mode]] <- stats[[nrm]][[mode]][["stats"]]
                }
            stats <- stats.1
        }
        deuids <- get.deuids(stats = stats, nrms = nrms, modes = modes, levels = levels, annots = annots, aligns = aligns, compares = compares, fc.fc = fc.fc, fc.avgexpr = fc.avgexpr, rs.dranksum = rs.dranksum, dsq.fc = dsq.fc, dsq.avgcpm = dsq.avgcpm, dsq.pval = dsq.pval, dsq.qval = dsq.qval)
        diffs <- list(stats = stats, deuids = deuids)
        new.diffs <- NULL
        for (nrm in nrms)
            for (mode in modes) {
                new.diffs[[nrm]][[mode]] <- list(stats = diffs[["stats"]][[nrm]][[mode]], deuids = diffs[["deuids"]][[nrm]][[mode]])
            }
        new.diffs
    }
    
    get.diffs.cuffdiff <- function(dir.in, levels, annots, aligns, compares, status = c("OK"), qval = NULL, pval = 0.05, check = TRUE) {
        # check the logs of cuffdiff results
        .levels <- paste(c("gene_exp", "isoform_exp", "promoters", "splicing"), "diff", sep = ".")
        names(.levels) <- c("gene", "transcript", "promoter", "splicing")
        check.cuffdiff <- function(dir.in, levels, annots, aligns, compares) {
            for (level in levels)
                for (annot in annots)
                    for (align in aligns)
                        for (compare in compares) {
                            cat(level, annot, align, compare, "\n")
                            din <- paste(dir.in, annot, align, compare, sep = "/")
                            fin <- paste(din, .levels[level], sep = "/")
                            cmd <- paste("wc -l", fin)
                            system(cmd)
                        }
        }
        if (check) check.cuffdiff(dir.in, levels, annots, aligns, compares)
        # test status: OK, NOTEST, LOWDATA, HIDATA, FAIL
        # fdr TRUE: only significant == yes
        # fdr FALSE: p-value < cutoff
        diffs <- NULL
        # ensembl doesn't have promoter, splicing results
        get.deuid <- function(stat, status, qval = 0.05, pval = 0.10) {
            rows <- rep(FALSE, nrow(stat))
            for (s in status) 
                rows <- rows | (stat[, "status"] == s)
            if (!is.null(pval))
                rows <- rows & (stat[, "p_value"] < pval)
            if (!is.null(qval))
                rows <- rows & (stat[, "q_value"] < qval)
            deuid <- rownames(stat)[rows]
            up <- deuid[stat[deuid, "reg"] == "up"]
            down <- deuid[stat[deuid, "reg"] == "down"]
            list(up = up, down = down, both = deuid)
        }
        for (level in levels)
            for (annot in annots)
                for (align in aligns) 
                    for (compare in compares) {
                        cat(level, annot, align, compare, "\n")
                        din <- paste(dir.in, annot, align, compare, sep = "/")
                        fin <- paste(din, .levels[level], sep = "/")
                        stat <- read.csv(fin, sep = "\t", head = TRUE, as.is = TRUE, check.names = FALSE) 
                        rownames(stat) <- stat[, "test_id"]
                        # add regulation column to gene/isoform exp stats
                        if (level == "gene" | level == "transcript") {
                            reg <- c("down", "nc", "up")[sign(stat[, "log2(fold_change)"]) + 2]
                            stat <- stat[, c("status", "value_1", "value_2", "log2(fold_change)", "test_stat", "p_value", "q_value")]
                            stat <- stat[order(stat[, "p_value"]), ]
                            stat <- data.frame(stat, reg = reg)
                        } 
                        deuid <- get.deuid(stat = stat, status = status, pval = pval, qval = qval)
                        diffs[["stats"]][[level]][[annot]][[align]][[compare]] <- list()
                        diffs[["deuids"]][[level]][[annot]][[align]][[compare]] <- list()
                        diffs[["stats"]][[level]][[annot]][[align]][[compare]] <- stat
                        diffs[["deuids"]][[level]][[annot]][[align]][[compare]] <- deuid
                    }
        diffs
    }
    # ============================================================
    # ANOVA
    # X = treatment + timepoint
    # ============================================================
    get.diffs.anova <- function(exprs, nrms = NULL, modes = NULL, levels, annots, designs, pval = 0.05, qval = NULL, adjust = "BH",  ftds = NULL, core = 5) {
        require("parallel")
        get.stat <- function(dat, designs) {
            dat.aov <- data.frame(dat = dat, designs)
            lm.aov <- lm(dat ~ tp + treat, data = dat.aov)
            stat <- anova(lm.aov)
            stat.pval <- stat[names(designs),"Pr(>F)"]
        }
        get.stats <- function(data, level, annot, designs, adjust = NULL, flags = NULL, absent = 1) {
            if (!is.null(flags)) {
                if (nrow(flags) == nrow(data)) {
                    rows <- rowSums(flags) >= absent
                    data <- data[-rows, ]
                }
            }
            rows <- apply(data, 1, sd) == 0
            data <- data[!rows, ]
            cl <- makeCluster(core)
            stats <- t(parApply(cl = cl, X = data, MARGIN = 1, FUN = function(d, designs) get.stat(d, designs), designs = designs))
            stopCluster(cl)
            #		stats <- t(apply(data, 1, get.stat, designs)) 
            colnames(stats) <- paste(names(designs), "pval", sep = ".")
            if (!is.null(adjust)) {
                stats.q <- apply(stats, 2, p.adjust, method = adjust)
                colnames(stats.q) <- paste(names(designs), "qval", sep = ".")
                stats <- cbind(stats, stats.q)
            }
            if (level == "isoforms") {
                stats <- data.frame(transript.id = rownames(data), stats)
            } else if (level == "genes") {
                if (annot == "igenome") {
                    ids <- t(sapply(rownames(data), function(x) strsplit(x, ";")[[1]], simplify = TRUE, USE.NAMES = FALSE))
                    stats <- data.frame(gene.id = ids[, 1], locus = ids[, 2], stats)
                }
                else 
                    stats <- data.frame(gene.id = rownames(data), stats)
            }
            stats
        }
        get.degs <- function(stats, designs, pval = 0.05, qval = NULL) {
            get.deg <- function(d, stats, pval, qval) {
                rows <- rep(TRUE, nrow(stats))
                if (!is.null(pval)) {
                    pcol <- paste(d, "pval", sep = ".")
                    rows <- rows & (stats[, pcol] < pval)
                }
                if (!is.null(qval)) {
                    qcol <- paste(d, "qval", sep = ".")
                    rows <- rows & (stats[, qval] < qval)
                }
                rs <- rownames(stats)[rows]
                rs <- rs[!is.na(rs)]
            }
            sapply(names(designs), get.deg, stats, pval, qval)
        }
        stats <- NULL
        degs <- NULL
        for (level in levels) 
            for (annot in annots)
                for (align in aligns) {
                    stats[[level]][[annot]][[align]] <- get.stats(exprs[[level]][[annot]][[align]], level = level, annot = annot, designs = designs, adjust = adjust, flags = flags[[level]][[annot]][[align]], absent = absent)
                    degs[[level]][[annot]][[align]] <- get.degs(stats[[level]][[annot]][[align]], pval = pval, qval = qval, designs = designs)
                }
        list(stats = stats, degs = degs)
    }

    get.sam.samobj <- function(exprs, factors, comparison, nperms = 100, nvals = 50, logged2 = TRUE) {
        require(samr)
        comps <- rev(strsplit(comparison, "2")[[1]])
        comp1 <- comps[1]; comp2 <- comps[2]
        idx1 <- which(factors == comp1); idx2 <- which(factors == comp2)
        expr <- exprs[, c(idx1, idx2)]
        y <- rep(c(1, 2), c(length(idx1), length(idx2)))
        ids <- rownames(exprs)
        data <- list(x = expr, y = y, geneid = ids, genenames = ids, logged2 = logged2)
        obj <- samr(data = data, resp.type = "Two class unpaired", nperms = nperms)
        dtable <- samr.compute.delta.table(obj)
        list(obj = obj, dtable = dtable, data = data)
    }
    get.sam.stat <- function(obj, dtable, data, delta, fc = 0, plot = FALSE, file.plot = NULL, save = FALSE, file.save = NULL) {
        sgtable <- samr.compute.siggenes.table(obj, delta, data, dtable, fc, all.genes = T)
        stat <- rbind(sgtable$genes.up, sgtable$genes.lo)
        ids <- stat[, "Gene ID"]
        df.stat <- data.frame(data$x[ids, ], 
                              "Score(d)" = as.numeric(stat[, "Score(d)"]),
                              "Numerator(r)" = as.numeric(stat[, "Numerator(r)"]),
                              "Denominator(s+s0)" = as.numeric(stat[, "Denominator(s+s0)"]),
                              "logFC" = log2(as.numeric(stat[, "Fold Change"])),
                              "qval" = as.numeric(stat[, "q-value(%)"]) / 100,
                              check.names = F)
        rownames(df.stat) <- ids
        df.stat <- df.stat[order(df.stat[, "Score(d)"], decreasing = T), ]
        if (plot && !is.null(file.plot)) {
            pdf(file = file.plot)
            tryCatch(samr.plot(obj, delta), error = function(e) {
                dev.off(); stop("error in plotting") },
                finally = { cat("ok plot\n"); dev.off() })
        }
        if (save && !is.null(file.save)) {
            write.csv(df.stat, file = file.save)
        }
        df.stat
    }
    
    get.sam.deg <- function(stat, score = NULL, qval = 0.4, fc = 2, file = NULL, save = FALSE) {
        idx <- rep(T, nrow(stat))
        ids <- rownames(stat)
        if (is.numeric(score) && score > 0) idx <- idx & abs(stat[, "Score(d)"]) >= score
        if (is.numeric(qval) && qval >= 0 && qval <= 1) idx <- idx & stat[, "qval"] < qval
        if (is.numeric(fc) && fc != 0) idx <- idx & abs(stat[, "logFC"]) >= log2(abs(fc))
        up <- ids[idx & stat[, "logFC"] > 0]
        down <- ids[idx & stat[, "logFC"] < 0]
        both <- ids[idx]
        deg <- list(up = up, down = down, both = both)
        if (save && !is.null(file)) {
            df.deg <- list2dataframe(deg, "")
            tryCatch(write.csv(df.deg, file = file, row.names = F),
                     error = function(e) {
                         browser() }, finally = cat("ok\n")
                    )
        }
        deg
    }
   for (obj in ls()) 
       assign(obj, get(obj), envir = Enrichments)
}
