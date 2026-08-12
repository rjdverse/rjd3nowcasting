#' @include utils.R
NULL

#' @importFrom rjd3jars check_java_version
.onAttach <- function(libname, pkgname) {
    # Check Java version
    rjd3jars::check_java_version(silent = FALSE, startup = TRUE)
}

#' @importFrom rJava .jpackage
#' @importFrom rjd3jars check_java_version reload_dictionaries
.onLoad <- function(libname, pkgname) {
    # Loading dependencies
    if (!requireNamespace("rjd3jars", quietly = TRUE)) {
        stop("Loading {rjd3jars} failed", call. = FALSE)
    }
    if (!requireNamespace("rjd3toolkit", quietly = TRUE)) {
        stop("Loading {rjd3toolkit} failed", call. = FALSE)
    }

    # Loading Java class
    jar_dir <- file.path(libname, pkgname, "inst", "java")
    jars_inst <- list.files(
        jar_dir,
        pattern = "\\.jar$",
        full.names = TRUE,
        all.files = TRUE
    )
    result <- rJava::.jpackage(
        pkgname,
        lib.loc = libname,
        morePaths = jars_inst
    )
    if (!result) {
        stop("Loading Java packages failed")
    }

    # Reload dictionnaries
    has_java <- rjd3jars::check_java_version(silent = TRUE)
    if (has_java) {
        rjd3jars::reload_dictionaries()
    }
}
