local({
    read.blood.mrna.celfile <- function(dir.in, samples, celfiles, save = FALSE, log2 = TRUE, dir.out = NULL) {
        require(affy)
        affy <- ReadAffy(celfile.path = dir.in, filenames = celfiles, sampleNames = samples)
        eset.mas5 <- mas5(affy)
        call.mas5 <- mas5calls(affy)
        calls <- exprs(call.mas5)
        exprs <- exprs(eset.mas5)
        if (isTRUE(log2)) exprs <- log2(exprs)
        if (isTRUE(save) && !is.null(dir.out)) {
            dir.create(dir.out, F, T)
            fout.call <- file.path(dir.out, "calls.csv")
            fout.expr <- file.path(dir.out, "exprs.csv")
            write.csv(file = fout.call, calls, row.names = T)
            write.csv(file = fout.expr, exprs, row.names = T)
        }
        list(exprs = exprs, calls = calls)
    }
    read.exprs.blood.mrna <- function(dir.in, file.expr, file.call, samples, allow.abs = 0.1) {
        fexpr <- file.path(dir.in, file.expr)
        fcall <- file.path(dir.in, file.call)
        exprs <- read.csv(fexpr, row.names = 1, as.is = T, check.names = F)
        calls <- read.csv(fcall, row.names = 1, as.is = T, check.names = F)
        absent <- apply(calls, 1, function(x) sum(x == "A")) 
        exprs <- exprs[absent <= as.integer(length(samples) * allow.abs), samples]
        exprs
    }
    read.exprs.seq.mrna <- function(expr.file, samples, avgexpr = NULL, delta = 1, log2 = TRUE) {
        expr <- read.csv(expr.file, row.names = 1, head = T, check.names = F, as.is = T)
        if (!is.null(avgexpr)) {
            flag <- rowMeans(expr) >= avgexpr
            expr <- expr[flag, ]
        }
        if (isTRUE(log2)) expr <- log2(expr + delta)
        expr[, samples]
    }
#    source("Source/functions/RA/DEGs.R", local = T) # equally
    sys.source("Source/functions/RA/DEGs.R", envir = environment())
    for (obj in c("read.blood.mrna.celfile", "read.exprs.blood.mrna", "read.exprs.seq.mrna", "DEGs"))
        assign(obj, get(obj), envir = RA)
})
