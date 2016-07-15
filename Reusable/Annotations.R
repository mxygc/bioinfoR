if (!exists("Annotations") || !is.environment(Annotations)) Annotations <- new.env(parent = emptyenv())
local({
    get.annot.affy <- function(annot.file = "Report/annotation/affy/rat230.2.csv") {
        annot <- read.csv(file = annot.file, as.is = T, check.names = F, head = T, row.names = NULL)
        annot <- annot[, c(1, 3, 5, 4)]
        colnames(annot) <- c("affy", "offsym", "refseq", "entrez")
        annot
    }
    get.annot.seq <- function(dir.in = "Report/annotation", assembly = "rn4", annotation = "ucsc") {
        annot.file <- paste(annotation, "csv", sep = ".")
        annot.file <- paste(dir.in, assembly, annot.file, sep = "/")
        annot <- read.csv(file = annot.file, as.is = T, check.names = F, head = T, row.names = NULL)
        annot <- annot[, c(1, 2)]
        colnames(annot) <- c("refseq", "offsym")
        annot
    }
    ##############################################################
    # updated get.annot.maps()
    # input: annot (a matrix: no rownames yet, columns are IDs; 
    ## key column, value columns) 
    # output: annot.maps(a list mapping key-value1, key-value2...)
    ##############################################################
    expand.ann <- function(ann, k, v, sep) {
        kids <- strsplit(ann[, k], sep)
        vids <- ann[, v]
        ann <- do.call(rbind, lapply(1:length(kids), function(i, kids, vids) 
                cbind(kids[[i]], vids[[i]]), kids = kids, vids = vids))
        ann <- "colnames<-"(ann, c(k, v))
    }
    get.amap <- function(annot, key, value, ksep = NULL, vsep = " /// ", osep = " /// ", kvoid = NULL, vvoid = "---", ovoid = "") {
        colnames <- colnames(annot)
        if (!key %in% colnames) stop("key column not found!")
        if (!value %in% colnames) stop("value column not found!")
        if (key == value) stop("key column is identical to value column!")
        # filter the empty keys
        if (!is.null(kvoid)) annot <- annot[annot[, key] != kvoid, ]
        if (!is.null(ksep)) annot <- expand.ann(ann = annot, k = key, v = value, sep = ksep)
        kids <- annot[, key]
        if (anyDuplicated(kids)) {
            ukids <- unique(kids)
            amap <- tapply(annot[, value], kids,
                function(x, vsep, osep, vvoid, ovoid) {
                    if (!is.null(vsep)) x <- unlist(strsplit(x, vsep))
                    if (!is.null(vvoid)) x <- x[x != vvoid]
                    ifelse(length(x) > 0, paste(unique(x), collapse = osep), ovoid)
                }, vsep = vsep, osep = osep, vvoid = vvoid, ovoid = ovoid)
            amap <- amap[ukids]
        } else {
            amap <- annot[, value]
            names(amap) <- kids
            if (vvoid != ovoid) amap <- gsub(vvoid, ovoid, amap)
        }
        amap
    }
    get.annot.map <- function(key2value, annot, sep.chars = list(affy = NULL, offsym = " /// ", entrez = " /// ", refseq = " /// "), sep.char.out = " /// ", void.chars = list(affy = "", offsym = "---", refseq = "---", entrez = "---"), void.char.out = "") {
        cat(key2value, "\n")
        kvs <- strsplit(key2value, "2")[[1]]
        key <- kvs[1]
        value <- kvs[2]
        if (all(c(key, value) %in% names(sep.chars))) {
            ksep <- sep.chars[[key]]
            vsep <- sep.chars[[value]]
        } else stop("at least key or value not in sep.chars!")
        osep <- sep.char.out
        if (all(c(key, value) %in% names(void.chars))) {
            kvoid <- void.chars[[key]]
            vvoid <- void.chars[[value]]
        } else stop("at least key or value not in void.chars!")
        ovoid <- void.char.out
        get.amap(annot = annot, key = key, value = value, 
                 ksep = ksep, vsep = vsep, osep = osep, 
                 kvoid = kvoid, vvoid = vvoid, ovoid = ovoid)
    }
    get.annot.maps <- function(annot, key2values, sep.chars = list(affy = NULL, offsym = " /// ", entrez = " /// ", refseq = " /// "), sep.char.out = " /// ", void.chars = list(affy = "", offsym = "---", refseq = "---", entrez = "---"), void.char.out = "") {
        sapply(key2values, get.annot.map, annot = annot, sep.chars = sep.chars, sep.char.out = sep.char.out, void.chars = void.chars, void.char.out = void.char.out, simplify = F, USE.NAMES = T)
    }
    get.id.maps <- function(annot.maps, sep = " /// ") {
    #in id.maps, the empty key is denoted by character(0), not by a void character (eg. "" or "---") 
        lapply(annot.maps, strsplit, split = sep)
    }
    assign("Annotations", environment(), envir = .GlobalEnv)
})

