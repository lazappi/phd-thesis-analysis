#' Outlier histogram
#'
#' Histogram with lines showing median and MADs from median
#'
#' @param data data.frame to plot
#' @param x Name of the column to plot
#' @param mads Vector of mads to show
#' @param bins Number of bins for histogram
#' @param show_zero Whether to show lines that have more extreme values
#'
#' @return ggplot object
outlierHistogram <- function(data, x, mads = 1:3, bins = 30,
                             show_zero = FALSE) {
    x_data <- data[, x]
    med    <- median(x_data)
    MAD    <- mad(x_data, center = med, na.rm = TRUE)

    outs <- lapply(mads, function(n) {
        lower  <- med - n * MAD
        higher <- med + n * MAD
        n_low  <- sum(x_data < lower)
        n_high <- sum(x_data > higher)

        list(n = n, lower = lower, higher = higher,
             n_low = n_low, n_high = n_high)
    })

    gg <- ggplot(data, aes(x = !!ensym(x))) +
        geom_histogram(bins = bins) +
        geom_vline(xintercept = med, colour = "blue") +
        annotate("text", x = med, y = 0, label = "Median", colour = "blue",
                 angle = 90, hjust = -0.1, vjust = -0.5) +
        annotate("text", x = med, y = Inf, label = round(med, 2), colour = "blue",
                 angle = 90, hjust = 1, vjust = -0.5) +
        theme_minimal() +
        theme(legend.position = "none")

    cols <- scales::brewer_pal(palette = "Reds",
                               direction = -1)(length(mads) + 1)

    for (i in seq_along(outs)) {
        out <- outs[[i]]

        if (show_zero) {
            show_low  <- TRUE
            show_high <- TRUE
        } else {
            show_low  <- out$n_low > 0
            show_high <- out$n_high > 0
        }

        if (show_low) {
            gg <- gg +
                geom_vline(xintercept = out$lower, colour = cols[i]) +
                annotate("text", x = out$lower, y = 0,
                         label = paste0("-", out$n, " MADs"), colour = cols[i],
                         angle = 90, hjust = -0.1, vjust = -0.5) +
                annotate("text", x = out$lower, y = Inf,
                         label = paste0(round(out$lower, 2),
                                        " (", out$n_low, " lower)"),
                         colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
        }

        if (show_high) {
            gg <- gg +
                geom_vline(xintercept = out$higher, colour = cols[i]) +
                annotate("text", x = out$higher, y = 0,
                         label = paste0("+", out$n, " MADs"), colour = cols[i],
                         angle = 90, hjust = -0.1, vjust = -0.5) +
                annotate("text", x = out$higher, y = Inf,
                         label = paste0(round(out$higher, 2),
                                        " (", out$n_high, " higher)"),
                         colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
        }
    }

    return(gg)
}
