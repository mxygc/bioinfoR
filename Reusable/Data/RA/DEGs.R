if (!exists("DEGs") || !is.environment("DEGs")) DEGs <- new.env(parent = emptyenv())#assign("DEGs", new.env(parent = emptyenv()), envir = RA)
local({
    get.yuanhua.degs.b1 <- function(dir.in, rnk = FALSE, id = "affy") { ## FC=2
        compares.1 <- c("ANE_A2ANE_M","ANE_A2NRA_A","ANE_A2PLA_A","ANE_A2RA_B","ANE_M2NRA_M","ANE_M2PLA_M","ANE_M2RA_B","APC_A2ANE_A","APC_A2APC_M","APC_A2MTX_A","APC_A2NRA_A","APC_A2PLA_A","APC_A2RA_B","APC_M2ANE_M","APC_M2MTX_M","APC_M2NRA_M","APC_M2PLA_M","APC_M2RA_B","MTX_A2ANE_A","MTX_A2MTX_M","MTX_A2NRA_A","MTX_A2PLA_A","MTX_A2RA_B","MTX_M2ANE_M","MTX_M2NRA_M","MTX_M2PLA_M","MTX_M2RA_B","NRA_A2NRA_B","NRA_A2NRA_M","NRA_M2NRA_B","PLA_A2NRA_A","PLA_A2PLA_M","PLA_A2RA_B","PLA_M2NRA_M","PLA_M2RA_B","RA_B2NRA_B")
        types.1 <- ifelse(unlist(lapply(strsplit(compares.1, "2"), function(x) strsplit(x[1], "_")[[1]][2] == strsplit(x[2], "_")[[1]][2])), yes = "between", no = "within")
        types.1[c(35, 36)] <- rep("within", 2)
        regs <- c("up", "down")
        seg1 <- "mRNA.limma."
        seg2 <- ".sum-mas5.analysistype-"
        seg3 <- ".FC-2.selected"
        seg4 <- ".pvalue-0.05.csv"
        des.1 <- NULL
        degs.1 <- NULL
        din <- paste(dir.in, "FC2", sep = "/")
        for (i in 1:length(compares.1)) {
            compare <- compares.1[i]
            type <- types.1[i]
               fin <- paste(seg1, compare, seg2, type, seg3, seg4, sep = "")
               fin <- paste(din, fin, sep = "/")
               cat(fin, "\n")
               des <- read.csv(file = fin, head = T, as.is = T, check.names = F)
               degs.1[[compare]][["both"]] <- list()
               degs.1[[compare]][["both"]] <- des[, ifelse(id == "offsym", "symbols", "probes")]
               degs.1[[compare]][["up"]] <- list()
               degs.1[[compare]][["up"]] <- des[des[, "FC"] > 0, ifelse(id == "offsym", "symbols", "probes")]
               degs.1[[compare]][["down"]] <- list()
               degs.1[[compare]][["down"]] <- des[des[, "FC"] < 0, ifelse(id == "offsym", "symbols", "probes")]
               des <- des[, c(ifelse(id == "offsym", "symbols", "probes"), "FC")]
               if (id == "offsym") {
                   des <- lapply(1:nrow(des), function (i){ g <- strsplit(des[i, 1], " /// ")[[1]]; g <- unique(g); cbind(symbols = g, FC = des[i, 2]) })
                   des1 <- NULL
                   for (i in 1:length(des)) des1 <- rbind(des1, des[[i]])
                   colnames(des1) <- c("symbols", "FC")
                   des1 <- des1[(des1[, "symbols"] != "---") & (des1[, "symbols"] != ""), ]
                   des1 <- tapply(as.numeric(des1[, "FC"]), des1[, "symbols"], FUN = "median")
                   des1 <- sort(des1, decreasing = T)
                   des2 <- data.frame(symbols = rownames(des1), FC = des1)
                   des.1[[compare]] <- des2
               } else {
                   des1 <- des[, "FC"]
                   names(des1) <- des[, "probes"]
                   des1 <- sort(des1, decreasing = T)
                   des2 <- data.frame(symbols = names(des1), FC = des1)
                   des.1[[compare]] <- des2
               }
        }
        if (rnk) return(des.1)
        else return(degs.1)
        invisible()
    }
    get.yuanhua.degs.b2.4 <- function(dir.in, FC = 2, rev.lab = TRUE, rm.apa9 = FALSE, rnk = FALSE, id = "offsym") {
    # FC 1.8 or 2
    # label reversed or not
    # APA9 removed or not
        compares.2 <- c("AP_A2RA_B", "HAP_A2NORA_B", "MAP_A2RA_B", "RA_B2NORA_B", 
                        "AP_A2HAP_A", "AP_A2MAP_A", "MAP_A2HAP_A")
        types.2 <- c(rep("within", 4), rep("between", 3))
        regs <- c("up", "down")
        seg1 <- "mRNA.limma."
        seg2 <- ".sum-mas5.analysistype-"
        seg3 <- paste(".FC-", FC, ".selected", sep = "")
        seg4 <- ".pvalue-0.05.Ath-0.3"
        seg5 <- ifelse(rev.lab, yes = "R", no = "NoR")
        seg6 <- ifelse(rm.apa9, yes = "AP_A.9-removed", no = "NULL-removed")
        des <- NULL
        degs <- NULL
        din <- paste(dir.in, paste("FC", FC, sep = ""), sep = "/")
        for (i in 1:length(compares.2)) {
            compare <- compares.2[i]
            type <- types.2[i]
            fin <- paste(seg1, compare, seg2, type, seg3, seg4, ".", seg5, ".", seg6, ".csv", sep = "")
            fin <- paste(din, fin, sep = "/")
            cat(fin, "\n")
            de <- read.csv(file = fin, head = T, as.is = T, check.names = F)
            de <- de[, c(ifelse(id == "offsym", "symbols", "probes"), "FC")]
            degs[[compare]][["up"]] <- list()
            degs[[compare]][["down"]] <- list()
            degs[[compare]][["both"]] <- list()
            degs[[compare]][["up"]]  <- de[de[, "FC"] > 0, ifelse(id == "offsym", "symbols", "probes")]
            degs[[compare]][["down"]] <- de[de[, "FC"] < 0, ifelse(id == "offsym", "symbols", "probes")]
            degs[[compare]][["both"]] <- de[, ifelse(id == "offsym", "symbols", "probes")]
            if (id == "offsym") {
                de <- lapply(1:nrow(de), function (i) { g <- strsplit(de[i, 1], " /// ")[[1]]; g <- unique(g); cbind(symbols = g, FC = de[i, 2]) })
                de1 <- NULL
                for (i in 1:length(de)) de1 <- rbind(de1, de[[i]])
                colnames(de1) <- c("symbols", "FC")
                de1 <- de1[(de1[, "symbols"] != "---") & (de1[, "symbols"] != ""), ]
                de1 <- tapply(as.numeric(de1[, "FC"]), de1[, "symbols"], FUN = "median")
                de1 <- sort(de1, decreasing = T)
                de2 <- data.frame(symbols = rownames(de1), FC = de1)
                des[[compare]] <- de2
            } else {
                de1 <- de[, "FC"]
                names(de1) <- de[, "probes"]
                de1 <- sort(de1, decreasing = T)
                de2 <- data.frame(probes = names(de1), FC = de1)
                des[[compare]] <- de2
            }
        }
        if (rnk) return(des)
        else return(degs)
        invisible()
    }
    add.early.degs.b1 <- function(cmps, degs, regs = c("up", "down", "both")) {
    # early degs blood batch 1
    # E2B = M2B-A2B = ANE_E2RA_B, APC_E2RA_B, PLA_E2RA_B, MTX_E2RA_B
        get.early <- function(y, degs, regs) {
            x <- strsplit(y, "_")[[1]][1]
            m <- paste(x, "M", sep = "_")
            a <- paste(x, "A", sep = "_")
            b <- paste("RA", "B", sep = "_")
            cmp.m <- paste(m, b, sep = "2")
            cmp.a <- paste(a, b, sep = "2")
            degs.m <- degs[[cmp.m]]
            degs.a <- degs[[cmp.a]]
            degs.e <- list()
            for (reg in regs) {
                degs.e[[reg]] <- setdiff(x = degs.m[[reg]], y = degs.a[[reg]])
            }
            degs.e
        }
        early.degs <- sapply(cmps, get.early, degs = degs, regs = regs, USE.NAMES = T, simplify = F)
        degs <- c(degs, early.degs)
    }
    add.early.cmp.degs.b1 <- function(cmps, degs, regs = c("up", "down", "both")) {
    # early time points comparisons blood batch 1
    # E2E = E2B-E2B = APC_E2MTX_E, APC_E2PLA_E, APC_E2ANE_E, MTX_E2PLA_E, MTX_E2ANE_E, ANE_E2PLA_E
        get.early.cmp <- function(z, degs, regs) {
            y <- strsplit(z, "2")[[1]][1]
            x <- strsplit(z, "2")[[1]][2]
            x <- paste(x, "RA_B", sep = "2")
            y <- paste(y, "RA_B", sep = "2")
            degs.x <- degs[[x]]
            degs.y <- degs[[y]]
            degs.e <- list()
            for (reg in regs) {
                degs.e[[reg]] <- setdiff(degs.y[[reg]], degs.x[[reg]])
            }
            degs.e
        }
        early.cmp.degs <- sapply(cmps, get.early.cmp, degs = degs, regs = regs, USE.NAMES = T, simplify = F)
        degs <- c(degs, early.cmp.degs)
    }
    get.degs.blood.mrna.batch1 <- function(id = "affy") {
        dir.in <- "Report/deg/blood/mrna/batch1"
        degs.1 <- get.yuanhua.degs.b1(dir.in = dir.in, rnk = F, id = id)
        degs.1 <- add.early.degs.b1(cmps = c("ANE_E2RA_B", "APC_E2RA_B", "PLA_E2RA_B", "MTX_E2RA_B"), degs = degs.1, regs = c("up", "down", "both"))
        degs.1 <- add.early.cmp.degs.b1(cmps = c("APC_E2MTX_E", "APC_E2PLA_E", "APC_E2ANE_E", "MTX_E2PLA_E", "MTX_E2ANE_E", "ANE_E2PLA_E"), degs = degs.1, regs = c("up", "down", "both"))
        names(degs.1) <- gsub(names(degs.1), pattern="APC", replacement = "AP")
        degs.1 <- degs.1[unlist(lapply(strsplit(names(degs.1), "2"), function(x) strsplit(x, "_")[[1]][2] != "M" & strsplit(x, "_")[[2]][2] != "M"))]
        degs.1 <- degs.1[c(18, 22, 23, 20, 19, 17, 21, 24, 25, 26, 16, 8, 12, 15, 3, 5, 4, 7, 6, 9, 10, 11, 14, 1, 2, 13)]
        degs.1
    }
    get.degs.blood.mrna.batch2 <- function(id = "affy") {
        dir.in <- "Report/deg/blood/mrna/batch2"
        degs.2 <- get.yuanhua.degs.b2.4(dir.in = dir.in, FC = 2, rev.lab = T, rm.apa9 = F, rnk = F, id = id)
    }
    get.degs.tissue.mrna.batch1 <- function(id = "refseq", cmps = c("AP_M2RA_B", "AP_A2RA_B", "AP_A2AP_M", "PLA_M2RA_B", "PLA_A2RA_B", "PLA_A2PLA_M", "AP_M2PLA_M", "AP_A2PLA_A"), regs = c("up", "down", "both")) {
        dir.in <- "Report/deg/tissue/mrna/batch1"
        file.in <- paste(id, "csv", sep = ".")
        file <- file.path(dir.in, file.in)
        degs.df <- read.csv(file = file, head = T, as.is = T)
        degs <- NULL
        for (cmp in cmps) 
            for (reg in regs) {
                l <- paste(cmp, reg, sep = "_")
                deg <- degs.df[, l]
                deg <- deg[deg != ""]
                degs[[cmp]][[reg]] <- list()
                degs[[cmp]][[reg]] <- deg
            }
        degs <- add.early.degs.b1(degs = degs, cmps = c("AP_E2RA_B", "PLA_E2RA_B"), regs = regs)
        degs <- add.early.cmp.degs.b1(degs = degs, cmps = c("AP_E2PLA_E"), regs = regs)
    }
    get.degs.tissue.mrna.batch2 <- function(id = "refseq", cmps = c("RA_B2NORA_B", "AP_A2RA_B", "AP_A2NORA_B", "AP_A2NORAP_A", "NORAP_A2NORA_B", "NORAP_A2RA_B"), regs = c("up", "down", "both")) {
        dir.in <- "Report/deg/tissue/mrna/batch2"
        file.in <- paste(id, "csv", sep = ".")
        file <- file.path(dir.in, file.in)
        degs.df <- read.csv(file = file, head = T, as.is = T)
        degs <- NULL
        for (cmp in cmps) 
            for (reg in regs) {
                l <- paste(cmp, reg, sep = "_")
                deg <- degs.df[, l]
                deg <- deg[deg != ""]
                degs[[cmp]][[reg]] <- list()
                degs[[cmp]][[reg]] <- deg
            }
        degs
    }
    for (obj in c("get.degs.blood.mrna.batch1", "get.degs.blood.mrna.batch2", 
                  "get.degs.tissue.mrna.batch1", "get.degs.tissue.mrna.batch2"))
        assign(obj, get(obj), envir = DEGs)
})
