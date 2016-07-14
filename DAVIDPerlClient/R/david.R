###############################################################################
# an R encapsulater of DAVID Web Service perl client:
# get david functional analysis results by submitting lists of gene ids
# http://david.abcc.ncifcrf.gov/webservice/sample_clients/PerlClient-1.1.zip
## run.david
## read.david
## save.david
###############################################################################

#' The Wapper
#' @param gene.ids the input gene list, either a character vector or a list of multiple character vectors
#' @param dir.out  the output directory; default is the current working directory
#' @param perl     the location of Perl excutable
#' @param auth.email   the email registed in DAVID Web Service
#' @param tasks    the query types, one or more from the options: TermEnrich, TermClust, GeneClust, GeneAnnot
#' @param species  the species of each input gene list; default is Homo sapiens
#' @param list.names the names of the input gene lists; default the same as the listnames of gene.ids
#' @param categories   
run.david <- function(
                      gene.ids, 
                      dir.out = NULL,
                      perl = "perl",
                      auth.email = "vishalrp@uci.edu", 
                      tasks = c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot"), 
                      species = "Homo sapiens", 
                      id.types = c("AFFYMETRIX_3PRIME_IVT_ID"), 
                      list.names = NULL, 
                      categories = c("OMIM_DISEASE", "COG_ONTOLOGY", "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE", "GOTERM_BP_FAT", "GOTERM_CC_FAT", "GOTERM_MF_FAT", "BBID", "BIOCARTA", "KEGG_PATHWAY", "INTERPRO", "PIR_SUPERFAMILY", "SMART"),
                      te.count = 2, te.ease = 0.1,
                      tc.similarity.gene.overlap = 3,
                      tc.initial.group.member = 3,
                      tc.final.group.member = 3,
                      tc.multi.linkage = 0.5,
                      tc.similarity = 0.5,
                      gc.similarity.term.overlap = 4,
                      gc.initial.group.member = 4,
                      gc.final.group.member = 4,
                      gc.multi.linkage = 0.5,
                      gc.similarity = 0.35,
                      verbose = FALSE,
                      envir = parent.frame()) {
    run.cmd <- function(name, ids, dir.out, perl, perl.lib, perl.exe, auth.email, tasks, species, id.type, categories,
            te.count, te.ease,
            tc.similarity.gene.overlap,
            tc.initial.group.member, 
            tc.final.group.member,
            tc.multi.linkage,
            tc.similarity,
            gc.similarity.term.overlap,
            gc.initial.group.member,
            gc.final.group.member,
            gc.multi.linkage,
            gc.similarity,
            verbose
            ) {
        if (id.type == "OFFICIAL_GENE_SYMBOL") 
            stop("sorry, 'OFFICIAL_GENE_SYMBOL' has been banned by DAVID because of its ambiguity (see http://david.abcc.ncifcrf.gov/forum/viewtopic.php?f=14&t=885&sid=d860a473f348c9bcd57e669ddf308fee)")
        if (length(ids) > 3000) stop("number of genes exceeds 3000!")
        file.ids <- paste(name, "txt", sep = ".")
        file.ids <- paste(dir.out, file.ids, sep = "/")
        dir.out <- paste(dir.out, name, sep = "/")
        dir.create(dir.out, recursive = T, showWarnings = F)
        cat(ids, file = file.ids, sep = "\n")
        if (!file.exists(file.ids)) {
            warnings("the file for input does not exist! waiting for 5 sec...")
            if (!file.exists(file.ids)) stop("the file for input still does not exist!")
        }
        .tasks <- c(TermEnrich = "term_enrich",
                    TermClust  = "term_clust",
                    GeneClust  = "gene_clust",
                    GeneAnnot  = "gene_annot")
        tasks <- .tasks[tasks]
        tasks <- unname(tasks[!is.na(tasks)])
        if (length(tasks) < 1) {
            warning("None of the 4 tasks recognized, trying default: 'TermEnrich', 'TermClust', 'GeneClust' and 'GeneAnnot'")
            tasks <- "term_enrich,term_clust,gene_clust,gene_annot"
        } else {
            tasks <- paste(tasks, collapse = ",")
        }
        categories <- paste(categories, collapse = ",")
        if (tc.similarity > 1 || tc.similarity < 0) stop("term clustering similarity parameter tc.similarity should bewteen 0 and 1!")
        if (gc.similarity > 1 || gc.similarity < 0) stop("gene clustering similarity parameter tc.similarity should bewteen 0 and 1!")
        tc.similarity <- round(tc.similarity * 100)
        gc.similarity <- round(gc.similarity * 100)
        cmd <- paste(perl,
                     "-I", shQuote(perl.lib),
                     "", shQuote(perl.exe),
                     "--input", file.ids, 
                     "--out_dir", dir.out, 
                     "--auth_email", shQuote(auth.email), 
                     "--task", shQuote(tasks),
                     "--species", shQuote(species),
                     "--id_type", shQuote(id.type),
                     "--list_name", shQuote(name),
                     "--category", shQuote(categories),
                     "--te_count", te.count,
                     "--te_ease", te.ease,
                     "--tc_similarity_gene_overlap", tc.similarity.gene.overlap,
                     "--tc_initial_group_member", tc.initial.group.member,
                     "--tc_final_group_member", tc.final.group.member,
                     "--tc_multi_linkage", tc.multi.linkage,
                     "--tc_similarity", tc.similarity,
                     "--gc_similarity_term_overlap", gc.similarity.term.overlap,
                     "--gc_initial_group_member", gc.initial.group.member,
                     "--gc_final_group_member", gc.final.group.member,
                     "--gc_multi_linkage", gc.multi.linkage,
                     "--gc_similarity", gc.similarity
                     )
        if (verbose) print(cmd)
        exit <- system(cmd)
        if (exit != 0) stop("Error happened!")
    }
    perl.lib <- file.path(path.package("DAVIDPerlClient"), "Perl")
    perl.exe <- file.path(perl.lib, "DAVIDPerlClient.pl")
    if (is.null(dir.out)) dir.out = getwd()
    else if (!any(grep("^\\/", dir.out))) dir.out <- paste(getwd(), dir.out, sep = "/")

    if (is.vector(gene.ids) & (!is.list(gene.ids))) {
        if (is.null(list.names)) list.names  <- "list_1"
        else if (length(list.names) > 1) {
            warning("Only one vector of gene ids found, but with more than one list names, trying list_1")
            list.names <- "list_1"
        }
        gene.ids <- list(gene.ids)
        names(gene.ids) <- list.names
        if (length(id.types) > 1) {
            warning("Only one vector of gene ids found, but with more than one id types, trying the first one")
            id.types <- id.types[1]
        }
        names(id.types) <- names(gene.ids)
    }
    else if (is.list(gene.ids)) {
        if (is.null(list.names)) {
            if (is.null(names(gene.ids))) list.names <- paste("list", 1:length(gene.ids), collapse = "_")
        } else {
            if (!is.null(names(gene.ids))) {
                if (length(gene.ids) == length(list.names)) { 
                    warning("names of the query list exists, also list names are provided, choosing the latter")
                    names(gene.ids) <- list.names
                } else {
                    warning("number of the list names provided do not match the number of the query lists, trying to preserve the former")
                }
            } else {
                if (length(gene.ids) == length(list.names)) names(gene.ids) <- list.names
                else {
                    warning("length of list names provided do not match the length of the query list, trying to rename them in order of list_1, list_2, ...")
                    names(gene.ids) <- paste("list", 1:length(gene.ids), collapse = "_")
                }
            }
        }
        if (length(id.types) == length(gene.ids)) names(id.types) <- names(gene.ids)
        else if (length(id.types) == 1) {
            warning("multiple gene lists but only one id type... assuming they are all ", id.types)
            id.types <- rep(id.types, length(gene.ids))
            names(id.types) <- names(gene.ids)
        } else 
            stop("the number of id types provided do not match the number of the query lists!")
    }
    for (name in names(gene.ids))
        run.cmd(name = name, ids = gene.ids[[name]], dir.out = dir.out, 
                perl = perl, perl.lib = perl.lib, perl.exe = perl.exe, 
                auth.email = auth.email, 
                tasks = tasks, species = species, id.type = id.types[name], 
                categories = categories, 
                te.count = te.count, te.ease = te.ease, 
                tc.similarity.gene.overlap = tc.similarity.gene.overlap, 
                tc.initial.group.member = tc.initial.group.member, 
                tc.final.group.member = tc.final.group.member, 
                tc.multi.linkage = tc.multi.linkage, 
                tc.similarity = tc.similarity, 
                gc.similarity.term.overlap = gc.similarity.term.overlap, 
                gc.initial.group.member = gc.initial.group.member, 
                gc.final.group.member = gc.final.group.member, 
                gc.multi.linkage = gc.multi.linkage, 
                gc.similarity = gc.similarity, 
                verbose = verbose)
}

write.xls <- function(sheets, file, sheet.names = NULL, row.names = T, col.names = T) {
    if (!require("WriteXLS")) { 
        install.packages("WriteXLS")
        if (!require("WriteXLS")) 
            stop("WriteXLS cannot be installed")
    }
    if (is.null(sheet.names)) sheet.names <- names(sheets)
    attach(sheets)
    WriteXLS(x = sheet.names, ExcelFileName = file, SheetNames = sheet.names, AdjWidth = F, row.names = row.names, col.names = col.names)
    detach(sheets)
}

convert.david <- function(david, map, task, how = "append") {
    conv <- function(old, map, how = "append", task) {
        if (task == "ga") id0 <- strsplit(old, ",")[[1]]
        else id0 <- strsplit(old, ", ")[[1]]
        id1 <- map[id0]
        if (how == "append") {
            id1 <- lapply(id1, function(x) paste(unique(x), collapse = ", "))
            new <- paste(id0, " (", id1, ")", sep = "", collapse = ", ")
        }
        else if (how == "replace") {
            id1 <- do.call(c, id1)
            id1 <- unique(id1[!id1 == ""])
            new <- paste(id1, sep = ", ", collapse = ", ")
        }
        new
    }
    if (task == "te") {
        if (nrow(david) > 0) {
            names(map) <- toupper(names(map))
#            map <- lapply(map, function(x) paste(unique(x), collapse = ", "))
            old <- david[-1, 6]
            new <- sapply(old, conv, map = map, how = how, task = task, USE.NAMES = F, simplify = T)
            david[-1, 6] <- new
        }
    } else if (task == "tc") {
        if (nrow(david) > 0) {
            names(map) <- toupper(names(map))
#            map <- lapply(map, function(x) paste(unique(x), collapse = ", "))
            col1 <- david[, 1]
            idx <- grep("Annotation Cluster", col1, fixed = T)
            flag <- rep(T, nrow(david))
            flag[idx] <- F
            flag[idx + 1] <- F
            flag[(idx - 1)[-1]] <- F
            old <- david[, 6]
            new <- sapply(1:nrow(david), function(i, flag, old, map, how, task) 
                          ifelse(flag[i], conv(old[i], map, how, task), old[i]), flag = flag, old = old, map = map, how = how, task = task, USE.NAMES = F, simplify = T)
            david[, 6] <- new
        }
    } else if (task == "gc") {
        if (nrow(david) > 0) {
#            map <- lapply(map, function(x) paste(unique(x), collapse = ", "))
            col1 <- david[, 1]
            idx <- grep("Gene Group", col1, fixed = T)
            flag <- rep(T, nrow(david))
            flag[idx] <- F
            flag[idx + 1] <- F
            flag[(idx - 1)[-1]] <- F
            old <- david[, 1]
            new <- sapply(1:nrow(david), function(i, flag, old, map, how, task) 
                          ifelse(flag[i], conv(old[i], map, how, task), old[i]), flag = flag, old = old, map = map, how = how, task = task, USE.NAMES = F, simplify = T)
            david[, 1] <- new
        }
    }
    david
}

read.david <- function(dir, tasks = c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot"), list.names, id.types = NULL, convert = F, maps = NULL, how = "append") {
    .tasks.in <- c("te", "tc", "gc", "ga")
    names(.tasks.in) <- c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot")
    .tasks.out <- c("term_enrich", "term_clust", "gene_clust", "gene_annot")
    names(.tasks.out) <- c("te", "tc", "gc", "ga")
    david <- list()
    if (convert) {
        if (is.null(id.types)) stop("No id types provided")
        if (is.null(maps)) stop("No maps provided")
        if (length(list.names) == length(id.types))
            names(id.types) <- list.names
        else if (length(id.types) == 1) {
            warning("multiple gene lists but only one id type... assuming they are all ", id.types)
            id.types <- rep(id.types, length(list.names))
            names(id.types) <- list.names
        } else 
            stop("the number of id types provided do not match the number of the query lists!")
        if (!all(id.types %in% names(maps))) 
            stop("id type: ", paste(id.types[!id.types %in% names(maps)], collapse = " "), " not in the maps")
    }
    for (task in .tasks.in[tasks])
        for (list.name in list.names) {
            dir.in <- paste(dir, list.name, sep = "/")
            file.in <- paste(.tasks.out[task], "txt", sep = ".")
            fin <- paste(dir.in, file.in, sep = "/")
            dvd <- tryCatch(read.csv(file = fin, sep = "\t", head = F, as.is = T, check.names = F, blank.lines.skip = F), error = function(e) { warning("reading an empty file: ", fin); data.frame(NULL) })
            if (convert) {
                if (is.null(id.types)) stop("No ID types provided")
                if (is.null(maps)) stop("No maps provided")
                id.type <- id.types[list.name]
                map <- maps[[id.type]]
                dvd <- convert.david(dvd, map, task, how)
            }
            david[[task]][[list.name]] <- list()
            david[[task]][[list.name]] <- dvd
        }
    david
}

save.david <- function(david, tasks = c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot"), file, indiv = TRUE, xlsx = FALSE) {
    .tasks.in <- c("te", "tc", "gc", "ga")
    names(.tasks.in) <- c("TermEnrich", "TermClust", "GeneClust", "GeneAnnot")
    .tasks.out <- c("term_enrich", "term_clust", "gene_clust", "gene_annot")
    names(.tasks.out) <- c("te", "tc", "gc", "ga")
    re <- regmatches(file, regexec("(.+)/([^/]*)$", file))[[1]]
    dir.out <- re[2]
    dir.create(dir.out, recursive = T, showWarnings = F)
    if (any(grepl("\\.xlsx$", file))) xlsx <- T
    else if (any(grepl("\\.xls$", file))) xlsx <- F
    file <- sub(pattern = "\\.xlsx?$", replacement = "", x = file, perl = T)
    if (indiv) {
        for (task in .tasks.in[tasks]) {
            fout <- paste(file, paste(.tasks.out[task], ifelse(xlsx, "xlsx", "xls"), sep = "."), sep = "_")
            dvd <- david[[task]]
            dvd <- lapply(dvd, function(m) if (!is.data.frame(m)) as.data.frame(m) else m) 
            write.xls(sheets = dvd, file = fout, row.names = F, col.names = F)
        }
    } else {
        fout <- paste(file, ifelse(xlsx, "xlsx", "xls"), sep = ".")
        listnames <- names(david[[1]])
        david1 <- list()
        for (listname in listnames)
            for (task in .tasks.in[tasks]) {
                sheetname <- paste(listname, .tasks.out[task], sep = "_")
                dvd <- as.data.frame(david[[task]][[listname]])
                david1 <- c(david1, list(dvd))
                names(david1)[length(david1)] <- sheetname
            }
        write.xls(sheets = david1, file = fout, row.names = F, col.names = F)
    }
}
