#' @include utils.R
NULL

.onAttach <- function(libname, pkgname) {
    if (rjd3toolkit::get_java_version() < rjd3toolkit::minimal_java_version) {
        packageStartupMessage(sprintf("Your java version is %s. %s or higher is needed.",
                                      rjd3toolkit::get_java_version(), rjd3toolkit::minimal_java_version))
    }
}

.onLoad <- function(libname, pkgname) {

    jar_dir <- system.file("java", package = "rjd3nowcasting")
    jars <- list.files(jar_dir, pattern = "\\.jar$", full.names = TRUE)
    rJava::.jaddClassPath(jars)

    result <- rJava::.jpackage(pkgname, lib.loc = libname)
    if (!result) stop("Loading java packages failed", call. = FALSE)

    if (rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version) {
        rjd3toolkit::reload_dictionaries()
    }

    #  proto.dir <- system.file("proto", package = pkgname)
    #  RProtoBuf::readProtoFiles2(protoPath = proto.dir)
}
