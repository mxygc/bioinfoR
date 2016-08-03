if (!exists("Tools") || !is.environment(Tools)) Tools <- new.env(parent = emptyenv())

local({
    object.sizes <- function() rev(sort(sapply(ls(envir = .GlobalEnv), function(x) object.size(get(x)))))

	resample <- function(x, ...) x[sample(length(x), ...)]
	
    capfirst <- function(s, strict = FALSE) {
		cap <- function(s) paste(toupper(substring(s, 1, 1)),
					{ s <- substring(s, 2); if(strict) tolower(s) else s },
					sep = "", collapse = " " )
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
	}
	
    read.xls <- function(file, sheet) {
		require("gdata")
		df <- gdata::read.xls(xls = file, sheet = sheet)
		detach("package:gdata")
		df
	}
	
	#' to be obsolete
	#' package `xlsx` is to replace
	write.xls <- function(sheets, file, sheet.names = NULL, row.names = T, col.names = T) {
        require("WriteXLS")
        if (is.null(sheet.names)) {
            sheet.names <- names(sheets)
            idx <- which(unlist(lapply(sheet.names, nchar)) > 31)
            sheet.names[idx] <- substr(sheet.names[idx], 1, 31)
        }
        attach(sheets)
        WriteXLS(x = sheet.names, ExcelFileName = file, SheetNames = sheet.names, AdjWidth = F, row.names = row.names, col.names = col.names)
        detach(sheets)
    }
    
	list2dataframe <- function(list, fill = "") {
        if (length(as.character(fill)) == 0) stop("blank symbol")
        maxlen <- max(unlist(lapply(list, length)))
        fill.empty <- function(x, maxlen, fill) 
            unname(c(x, rep(fill, maxlen - length(x))))
        as.data.frame(sapply(list, fill.empty, maxlen = maxlen, fill = fill, simplify = T, USE.NAMES = F), stringsAsFactors = F, row.names = NULL)
    }
    join.lists <- function(a, b) {
        if (!is.list(a)) stop("'a' is not a list")
        if (!is.list(b)) stop("'b' is not a list")
        if (is.null(names(a))) stop("'a' has no names")
        if (is.null(names(b))) stop("'b' has not names")
        names.a <- names(a)
        names.b <- names(b)
        names.a.unique <- names.a[!names.a %in% names.b]
        names.b.unique <- names.b[!names.b %in% names.a]
        names.shared <- intersect(names.a, names.b)
        if (length(names.shared) > 0) shared <- sapply(names.shared, function(n, a, b) unique(c(a[[n]], b[[n]])), USE.NAMES = T, simplify = F, a = a, b = b)
        else shared <- NULL
        c(a[names.a.unique], shared, b[names.b.unique])
    }
	dfrapply <- function(L, FUN, ...) {
        if (inherits(L, "data.frame")) FUN(L, ...)
        else lapply(L, dfrapply, FUN, ...)
    }
	percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
    }
    
	# java doesn't recognize cygwin/unix path,
    # so it has to be patched by this function
    unix2dos <- function(path.unix) {
        if (R.Version()$os == "cygwin")
            path.dos <- sprintf("$(cygpath -w %s)", path.unix)
        else 
            stop("seems not running in cygwin...")
        return(path.dos)
    }
    
	assign("Tools", environment(), envir = .GlobalEnv)
})
