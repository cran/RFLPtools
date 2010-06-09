Gendist <- function (x, method = "rogers", diag = FALSE, upper = FALSE){
    if (!is.na(pmatch(method, "rogers"))) method <- "rogers"
    METHODS <- c("rogers", "sanghvi", "prevosti",
                 "nei1", "nei2", "nei3", "slatkin", "shriver")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    N <- nrow(x <- as.matrix(x))
    if(method == 1){## Rogers distance (1972)
        return(dist(x, diag = diag, upper = upper)/sqrt(2))
    }
    d <- .C("R_Gendistance", x = as.double(x), nr = N, nc = ncol(x), 
            d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
            method = as.integer(method-1), DUP = FALSE, 
            NAOK = TRUE, PACKAGE = "RFLPtools")$d      
    d <- d[lower.tri(d)]
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1L]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
