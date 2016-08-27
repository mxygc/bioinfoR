#' An example of batch downloading GEO data by GEOquery

if (!exists("PC") || !is.environment(PC)) PC <- new.env(parent = emptyenv())
local({
    get.exprs.gse <- function(gse = 'GSE14794') {
        get.metadata <- function(gsm) {
            meta <- Meta(gsm)
            with(meta, c(nchan = channel_count, 
                         cline = gsub(pattern = ".*: (.*?)", replacement = "\\1", x = characteristics_ch1[1]),
                         grade = gsub(pattern = ".*: (.*?)", replacement = "\\1", x = characteristics_ch1[2]),
                         nrow = data_row_count,
                         geoid = geo_accession,
                         platform = platform_id,
    #                     source = source_name_ch1,
                         mol = gsub(pattern = ".*\\((.*?),.*", replacement = "\\1", x = title),
                         sid = gsub(pattern = ".*\\(.*, (.*?) (.*?)\\)", replacement = "\\1_\\2", x = title)))
        }
        get.exprdata <- function(gse, meta) {
            gsms <- list(mrna = rownames(metadata)[metadata[, 'platform'] == 'GPL6102'],
                         mirna= rownames(metadata)[metadata[, 'platform'] == 'GPL8178'])
            get.expr <- function(gse, gsm, meta) {
                probe <- as.character(gse@gsms[[gsm[1]]]@dataTable@table[, 'ID_REF'])
                for (gs in gsm[-1]) probe <- intersect(probe, as.character(gse@gsms[[gs]]@dataTable@table[, 'ID_REF']))
                if (length(probe) < 1) stop('no common probes')
                expr.unord <- sapply(gsm, function(gs) { 
                       expr <- gse@gsms[[gs]]@dataTable@table; 
                       rownames(expr) <- expr[, 'ID_REF'];
                       expr }, USE.NAMES = T, simplify = F)
                expr <- sapply(gsm, function(gs) as.numeric(expr.unord[[gs]][probe, 'VALUE']), USE.NAMES = T,  simplify = T)
                sid <- meta[colnames(expr), 'sid']
                rownames(expr) <- probe
                colnames(expr) <- sid
                expr
            }
            sapply(gsms, function(gsm, gse, meta) get.expr(gse, gsm, meta), gse = gse, meta = metadata, USE.NAMES = T, simplify = F)
        }
        if (!is.object(gse)) {
            require('GEOquery')
            gse <- getGEO(GEO = 'GSE14794')
        }
        annots <- list(mrna = gse@gpls$GPL6102@dataTable@table,
                       mirna= gse@gpls$GPL8178@dataTable@table)
        metadata <- t(sapply(GSMList(gse), get.metadata, USE.NAMES = T, simplify = T))
        exprdata <- get.exprdata(gse, metadata)
        metadata <- list(mrna = metadata[metadata[, 'platform'] == 'GPL6102', ],
                         mirna = metadata[metadata[, 'platform'] == 'GPL8178', ])
        list(annots = annots, meta = metadata, exprs = exprdata)
    }
    assign("PC", environment(), envir = .GlobalEnv)
})


