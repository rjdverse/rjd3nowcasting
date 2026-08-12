#' @include utils.R
#' @importFrom stats ts.plot start
#' @importFrom graphics axis legend par barplot points lines
#' @importFrom grDevices palette.colors colorRampPalette
NULL

#' @export
print.JD3_DFMESTIMATES <- function(x, ...) {
    cat("DFM Estimates\n")
    cat(strrep("-", 20), "\n", sep = "")
    cat(sprintf("Has converged: %s\n", toupper(x$has_converged)))
    cat(sprintf("Log-likelihood: %.4f\n", x$log_likelihood))

    invisible(x)
}

#' @export
print.JD3_DFMRESULTS <- function(x, ...){

    sample_mean_stdev <- round(x$preprocessing$sample_mean_stdev, 3)
    param_mcoeff <- round(x$parameters$measurement_coefficients, 3)
    param_mvar <- round(x$parameters$measurement_errors_variance, 3)

    loadings <- cbind(sample_mean_stdev, param_mcoeff, param_mvar)
    colnames(loadings) <- c(
        "Sample mean",
        "Stdev",
        paste0("Coeff. ", colnames(param_mcoeff)),
        "Idiosyncratic variance"
    )

    var_coef <- round(x$parameters$var_coefficients, 3)
    var_err_variance <- round(x$parameters$var_errors_variance, 3)

    cat("DFM Results\n")
    cat(strrep("-", 20), "\n", sep = "")

    cat("\nLoadings:\n")
    print(loadings)
    cat("\nVAR Model:\n")
    print(var_coef)
    cat("\nInnovative Variance:\n")
    print(var_err_variance)

    invisible(x)
}

#' @export
print.JD3_DFMFORECASTS <- function(x, ...) {
    cat("DFM Forecasts (forecasts only)\n")
    cat(strrep("-", 30), "\n", sep = "")
    print(x$forecasts_only)

    invisible(x)
}

#' @export
print.JD3_DFMNEWS <- function(x, ...){
    cat("DFM NEWS\n")
    cat(strrep("-", 20), "\n", sep = "")

    summary(x)

    cat("\nForecasts:\n")
    print(round(x$forecasts, 3))

    invisible(x)
}

#' @export
summary.JD3_DFMNEWS <- function(object, ...) {
    x <- object

    nf <- ncol(x$forecasts)

    base <- cbind(
        x$impacts[, 1:2],
        round(as.matrix(sapply(x$impacts[, 3:5], as.numeric)), 3)
    )

    # Extract weights and impacts
    weights <- as.matrix(sapply(x$weights[, 6:ncol(x$weights), drop = FALSE], as.numeric))
    weights <- rbind(weights, colSums(weights))
    impacts <- as.matrix(sapply(x$impacts[, 6:ncol(x$impacts), drop = FALSE], as.numeric))

    # Interleave weights and impacts
    wi_mat <- matrix(NA, nrow = nrow(weights), ncol = 2 * nf)
    wi_mat[, seq(1, 2 * nf, by = 2)] <- round(weights, 3)
    wi_mat[, seq(2, 2 * nf, by = 2)] <- round(impacts, 3)

    colnames(wi_mat) <- as.vector(rbind(
        colnames(x$weights)[6:(5 + nf)],
        colnames(x$impacts)[6:(5 + nf)]
    ))

    # Final summary
    summary_news <- cbind(base, wi_mat)

    cat("\nTarget series:\n")
    cat(x$target_series, "\n")
    cat("\nNews analysis:\n")
    print(summary_news, row.names = FALSE)

    invisible(x)
}

#' @export
plot.JD3_DFMFORECASTS <- function(x, series_name = NULL, ...){

    fcst <- x$forecasts
    fcst_stdev <- x$forecasts_stdev
    fcst_only <- x$forecasts_only

    if (is.null(series_name)) {
        series_name <- colnames(fcst)[1]
    }

    if (!(series_name %in% colnames(fcst))) {
        stop("series name not found!")
    }

    s <- fcst[, series_name]
    s_lb <- s - 1.28 * fcst_stdev[, series_name]
    s_ub <- s + 1.28 * fcst_stdev[, series_name]
    sf <- fcst_only[, series_name]

    stats::ts.plot(
        s_lb,
        s_ub,
        s,
        sf,
        gpars = list(
            main = series_name,
            sub = "Forecasts with a 80% prediction interval",
            xlab = "",
            ylab = "",
            lty = c(3, 3, 1, 1),
            xaxt = "n",
            type = "o",
            pch = 20,
            cex = 0.8,
            las = 2,
            col = c("orange", "orange", "black", "red")
        )
    )
    graphics::axis(1,
                   at = seq(stats::start(s)[1], stats::end(s)[1], by = 1),
                   las = 2)
    graphics::legend(
        "topleft",
        legend = c("series", "forecasts", "80% PI"),
        col = c("black", "red", "orange"),
        lty = c(1, 1, 3),
        cex = 0.8
    )

    invisible(x)
}

#' @export
plot.JD3_DFMNEWS <- function(x, ...) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    start_col <- 6
    impacts <- x$impacts
    fcsts <- x$forecasts

    nr <- nrow(impacts)
    nf <- ncol(fcsts)

    # Extract impacts (excluding total row)
    impacts_mat <- as.matrix(
        sapply(impacts[-nr, start_col:(start_col + nf - 1)], as.numeric)
    )
    colnames(impacts_mat) <- colnames(fcsts)
    rownames(impacts_mat) <- impacts[-nr, 1]

    # Total impacts
    total_impacts <- as.numeric(
        impacts[nr, start_col:(start_col + nf - 1)]
    )

    # Colors
    plot_colors <- if (nr <= 10) {
        grDevices::palette.colors(nr - 1)
    } else {
        grDevices::colorRampPalette(grDevices::palette.colors(8))(nr - 1)
    }

    # Legend text
    legend_text <- paste0(impacts[-nr, 1], " (", impacts[-nr, 2], ")")
    legend_text <- c(legend_text, "Total impact")

    graphics::par(mar = c(5.1, 4.1, 4.1, 16.1), xpd = TRUE)

    # Handle positive and negative contributions separately
    impacts_positive <- impacts_mat
    impacts_negative <- impacts_mat

    impacts_positive[impacts_positive < 0] <- 0
    impacts_negative[impacts_negative > 0] <- 0

    # Plotting range
    impacts_range <- c(
        min(colSums(impacts_negative)),
        max(colSums(impacts_positive))
    )

    # Barplot (negative first)
    bar_positions <- graphics::barplot(
        impacts_negative,
        ylim = impacts_range,
        main = paste0("Impacts of news on forecast revisions (", x$target_series, ")"),
        col = plot_colors,
        las = 1,
        ...
    )

    # Overlay positive values
    graphics::barplot(
        impacts_positive,
        add = TRUE,
        ylim = rev(impacts_range),
        col =  plot_colors,
        las = 1,
    )

    # Overlay total impact
    graphics::points(bar_positions, total_impacts, pch = 3, col = "red", cex = 1)
    graphics::lines(bar_positions, total_impacts, col = "red")

    # Legend
    graphics::legend(
        "bottomright",
        legend = legend_text,
        inset = c(-0.6, 0),
        col = c(plot_colors, "red"),
        pch = c(rep(15, nr - 1), 3)
    )

    invisible(x)
}
