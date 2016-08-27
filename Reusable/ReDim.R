#' Scripts for demension reduction

reduce.dims <- function(exprs, dir, file, grpcols, levels, annots, aligns, type = "pdf") {
    # ============================================================
    # MDS
    # ============================================================
    plot.mds <- function(data, grpcols, level, annot, align) {
        require("MASS")
        require("stats")
        dist <- as.dist(1 - cor(data))
        mds.classic <- cmdscale(d = dist, k = 2)
        mds.nonmetric <- isoMDS(d = dist, k = 2)$points
        mains <- paste(c("classic MDS", "nonmetric MDS"), paste(level, annot, align, sep = " "), sep = "\n")
        plot(mds.classic[, 1], -mds.classic[, 2], type = "n", xlab = "", ylab = "", xlim = range(mds.classic[, 1]) + c(-0.02, 0.02), ylim = range(-mds.classic[, 2]) + c(-0.02, 0.02), asp = 1, axes = TRUE)
        text(x = mds.classic[, 1], y = -mds.classic[, 2], colnames(data), col = grpcols)
        title(main = mains[1], sub = "")
        plot(mds.nonmetric, xlim = range(mds.nonmetric[, 1]) + c(-0.02, 0.02), ylim = range(mds.nonmetric[, 2]) + c(-0.02, 0.02), type = "n", xlab = "", ylab = "", asp = 1, axes = TRUE)
        text(mds.nonmetric, colnames(data), col = grpcols)
        title(main = mains[2], sub = "")
        NULL
    }
    # ============================================================
    # PCA
    # ============================================================
    plot.pca <- function(data, nf = 5, grpcols, level, annot, align) {
        require("psych")
        mains <- paste(c("PCA classic", "PCA varimax"), paste(level, annot, align, sep = " "), sep = "\n")
        pca.cl <- princomp(data, cor = TRUE)
        plot(pca.cl$loadings, xlim = range(pca.cl$loadings[, 1]) + c(-0.02, 0.02), ylim = range(pca.cl$loadings[, 2]) + c(-0.02, 0.02), type = "n")
        text(pca.cl$loadings, colnames(data), col = grpcols)
        title(main = mains[1], sub = "")
        pca.rt <- principal(data, nfactors = nf, rotate = "varimax")
        plot(pca.rt$loadings, xlim = range(pca.rt$loadings[, 1]) + c(-0.02, 0.02), ylim = range(pca.rt$loadings[, 2]) + c(-0.02, 0.02), type = "n")
        text(pca.rt$loadings, colnames(data), col = grpcols)
        title(main = mains[2], sub = "")
        NULL
    }
    # ============================================================
    # FA
    # ============================================================
    plot.fa <- function(data, nf = 5, grpcols, level, annot, align) {
        require("HDMD")
        main <- paste(paste("FactAnal varimax", nf, "factors", sep = " "), paste(level, annot, align, sep = " "), sep = "\n")
#        fa <- factanal(data, nf, rotation = "varimax") 
# this has problem with merged htseq data, error: 
# unable to optimize from these starting value(s)
        fa <- factor.pa.ginv(data, nfactors = nf, rotate = "Promax")
        loading <- fa$loadings[, 1:2]
        plot(loading, xlim = range(loading[, 1]) + c(-0.02, 0.02), ylim = range(loading[, 2]) + c(-0.02, 0.02), type = "n")
        text(loading, colnames(data), col = grpcols)
        title(main = main, sub = "")
        NULL
    }
    reduce.dim <- function(expr, grpcols, level, annot, aligns) {
        for (align in aligns) {
            e <- expr[[align]]
            plot.mds(e, grpcols, level, annot, align)
            plot.pca(e, nf = 5, grpcols, level, annot, align)
            plot.fa(e, nf = 5, grpcols, level, annot, align)
        }
    }
    for (level in levels) 
        for (annot in annots) {
            expr <- exprs[[level]][[annot]]
            if(!is.null(type)) {
                dir.create(path = dir, recursive = TRUE, showWarnings = FALSE)
                dn <- paste(dir, paste(file, level, annot, sep = "_"), sep = "/")
                fn <- paste(dn, type, sep = ".")
                m <- length(aligns)
                if(type == "ps") {
                    postscript(file = fn, width = 5 * m, height = 5 * 5)
                } else if (type == "pdf") {
                    pdf(file = fn, width = 5 * m, height = 5 * 5)
                }
                par(mfcol = c(5, m))
                reduce.dim(expr, grpcols, level, annot, aligns)
                dev.off()
            }
        }
    NULL
}
# ============================================================
# metrics: pcc, euclidiean, maximun, manhattan, canberra, 
# binary, minkowski
# ============================================================
cluster.exprs <- function(exprs, dir, file, levels, annots, aligns, metrics = NULL, cols = NULL, type = "pdf") {
    get.hcluster <- function(exprs, metrics, cols = NULL, type = NULL, dir = NULL, level, annot, align) {
        if(!is.null(cols)) exprs <- exprs[, cols]
        hclust.1 <- function(x, m) { 
            if(m == "pcc") hclust(as.dist(1 - cor(x)))
            else hclust(dist(t(x), method = m))
        }
        hc <- sapply(metrics, function(m) hclust.1(exprs, m), USE.NAMES = TRUE, simplify = FALSE)
        if (!is.null(type)) { 
            for (m in names(hc)) {
                plot(hc[[m]], main = paste(paste(level, annot, align, sep = " "), m, sep = "\n"), sub = "", cex.sub = 1.2, xlab = "")
            }
        }
        hc
    }
    hcluster <- NULL
    if (is.null(metrics)) {
        metrics <- c("pcc", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
    }
    for (level in levels) 
        for (annot in annots) {
            hcluster[[level]][[annot]] <- list()
            if(!is.null(type)) {
                dir.create(path = dir, recursive = TRUE, showWarnings = FALSE)
                dn <- paste(dir, paste(file, level, annot, sep = "_"), sep = "/")
                fn <- paste(dn, type, sep = ".")
                n <- length(metrics)
                m <- length(aligns)
                if(type == "ps") {
                    postscript(file = fn, width = 5 * m, height = 5 * n)
                } else if (type == "pdf") {
                    pdf(file = fn, width = 5 * m, height = 5 * n)
                }
                par(mfcol = c(n, m))
            }
            for (align in aligns){
                hcluster[[level]][[annot]][[align]] <- get.hcluster(exprs[[level]][[annot]][[align]], type = type, level = level, annot = annot, align = align, cols = cols, metrics = metrics)
            }
            if (!is.null(type)) {
                dev.off()
            }
        }
    hcluster
}

