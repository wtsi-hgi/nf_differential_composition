#!/usr/bin/env Rscript

SCRIPT_NAME <- "make_cell_metric_plots.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tidyverse"))

set.seed(0)


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--results_files"),
            type = "character",
            default = "",
            help = paste0(
              "Tab-delimited results file.",
              " Must contain the following columns:",
              " cell_label,condition_effect,pvalue,method.",
              " Optional columns:",
              " formula."
            )
        ),

        optparse::make_option(c("--out_file"),
            type = "character",
            default = "differential_composition",
            help = paste0(
                "Name (and possibly path) of output file. Will have tsv.gz",
                " appended to it.",
                " [default: %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Compares results from differential composition tests."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    # Set input files
    files_string <- param[["results_files"]]
    base <- param[["out_file"]]

    # Merge all of the results
    files <- unlist(strsplit(files_string, ","))
    df_list <- list()
    i <- 1
    for (file in files) {
        df_list[[i]] <- data.table::fread(
            cmd = paste("gunzip -c", file),
            sep = "\t",
            header = T,
            stringsAsFactors = F
        )
        i <- i + 1
    }
    df_results <- data.table::rbindlist(df_list, use.names = TRUE, fill = TRUE)

    # Need to unique results because could have some duplicate results in
    # forward model selection.
    df_results <- unique(df_results)
    # Check that we don't have duplicates.
    df_results$check_var <- paste(
        df_results$variable,
        df_results$variable_effect, # needed if multiple factors e.g.: v1,v2,v3
        df_results$variable_target,
        df_results$formula,
        df_results$method,
        sep = "-"
    )
    if (length(unique(table(df_results$check_var))) > 1) {
        stop("ERROR: not unique results in merge.")
    }
    df_results$check_var <- paste(
        df_results$variable,
        df_results$variable_effect,
        df_results$formula,
        df_results$method,
        sep = "-"
    )
    if (length(unique(table(df_results$check_var))) > 1) {
        warning(paste0(
            "WARNING: Same model fit multiple times with different target",
            " variable. This could happen due to forward model selection."
        ))
    }

    # Save the output file
    df_results <- df_results[
        with(df_results, order(pvalue, abs(estimate))),
    ]
    gzfh <- gzfile(
        paste0(base, "-all_results.tsv.gz"),
        "w",
        compression = 9
    )
    write.table(
        df_results,
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
    close(gzfh)

    df_results <- data.frame(df_results)
    if ("variable_effect" %in% names(df_results)) {
        if (sum(df_results$variable_effect == "FALSE", na.rm=T) != 0) {
            stop("FALSE in variable_effect (misinterpreted F?).")
        }
    }

    # For development
    #df_results <- subset(df_results, condition != "hallmark_inflammatory_response")

    # Plot the results
    df_results <- subset(df_results, variable == variable_target)
    df_results$method_model <- paste(
        df_results$method,
        df_results$formula,
        sep = ":"
    )
    df_results$method_model <- unlist(lapply(
        df_results$method_model, FUN=function(x) {
            gsub(" ", "", gsub(":NA$", "", gsub(":$", "", x)))
        }
    ))
    df_results$variable_and_effect <- paste(
        df_results$variable_target,
        df_results$variable_effect,
        sep = ":"
    )
    filt <- df_results$variable_target == df_results$variable_effect
    df_results$variable_and_effect[filt] <- df_results$variable_target[filt]
    df_results$variable_and_effect <- unlist(lapply(
        df_results$variable_and_effect, FUN=function(x) {
            gsub(" ", "", gsub(":NA$", "", gsub(":$", "", x)))
        }
    ))

    # For each model, plot the qvalues across cell types
    if ("qvalue_bh" %in% colnames(df_results)) {
        # make a list of results grouped by method and formula
        df_results_list <- list(NULL)
        method_model <- unique(df_results$method_model)
        j <- 1
        for (i in method_model) {
            df_results_list[[j]] <- subset(df_results, method_model == i)
            j <- j + 1
        }

        for (df_results_i in df_results_list) {
            # Sort the dataframe
            df_results_i$hit <- df_results_i$qvalue_bh < 0.05
            df_results_i <- df_results_i[
                with(df_results_i, order(pvalue, abs(estimate))),
            ]
            df_results_i$cell_type <- factor(
                df_results_i$cell_type,
                levels = unique(df_results_i$cell_type)
            )
            plt <- ggplot2::ggplot(df_results_i, ggplot2::aes(
                x = cell_type,
                y = -log10(pvalue),
                fill = hit,
                label = variable_and_effect
            ))
            plt <- plt + ggplot2::theme_bw(base_size = 12)
            plt <- plt + ggplot2::geom_bar(stat = "identity")
            # plt <- plt + ggplot2::geom_hline(
            #     yintercept = -log10(0.05),
            #     linetype = "dashed"
            # )
            plt <- plt + ggplot2::labs(
                x = "Cell type",
                y = "-log10(p-value)",
                fill = "FDR < 5%",
                title = paste0(
                    unique(df_results_i$method_model),
                    "\n",
                    unique(df_results_i$variable_and_effect)
                )
            )
            plt <- plt + ggplot2::theme(
                # legend.position = "none",
                axis.text.x = ggplot2::element_text(angle = -90, hjust = 0)
            )
            plt <- plt + ggplot2::scale_fill_brewer(palette = "Dark2")
            nrow <- 1.5
            if (length(unique(df_results_i$variable_effect)) > 2) {
                plt <- plt + ggplot2::facet_wrap(~ variable_effect, ncol = 1)
                nrow <- 1.5 * length(unique(df_results_i$variable_effect)) / 2
            }
            cowplot::save_plot(
                filename = paste(
                    base,
                    unique(df_results_i$method_model),
                    "hits.png",
                    sep = "-"
                ),
                plot = plt,
                nrow = nrow,
                ncol = 2
            )
        }
    }


    # x axis = cluster, y = -log10p, color = trait, facet = model
    # pdf(
    #     file = paste(base, "all_results-wide.pdf", sep = "-"),
    #     height = 15,
    #     width = 25
    # )
    plt <- ggplot2::ggplot(df_results, ggplot2::aes(
        x = cell_type,
        y = -log10(pvalue),
        color = variable_and_effect,
        label = variable_and_effect
    ))
    plt <- plt + ggplot2::geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    )
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::geom_text(alpha = 0.75)
    # plt <- plt + ggplot2::scale_y_continuous(trans="-log10")
    # plt <- plt + ggplot2::theme(legend.position = "none")
    plt <- plt + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
    )
    plt <- plt + ggplot2::labs(
        x = "Cell type",
        y = "-log10(p-value)",
        fill = NULL
    )
    plt <- plt + ggplot2::facet_wrap(~ method, ncol = 1, scales = "free")
    cowplot::save_plot(
        filename = paste(
            base, "all_results-wide.pdf", sep = "-"
        ),
        plot = plt,
        nrow = 3,
        ncol = 2
    )
    # print(plt)
    # dev.off()
    # pdf(
    #     file = paste(base, "all_results-wide.pdf", sep = "-"),
    #     height = 15,
    #     width = 10
    # )
    plt <- plt + ggplot2::coord_flip()
    plt <- plt + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0)
    )
    cowplot::save_plot(
        filename = paste(
            base, "all_results-long.pdf", sep = "-"
        ),
        plot = plt,
        nrow = 2,
        ncol = 3
    )
    # print(plt)
    # dev.off()
    plt2 <- ggplot2::ggplot(df_results, ggplot2::aes(
        x = method,
        y = -log10(pvalue),
        color = variable_and_effect,
        label = variable_and_effect
    ))
    plt2 <- plt2 + ggplot2::geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    )
    plt2 <- plt2 + ggplot2::facet_wrap(
        ~ cell_type,
        scales = "free"
    )
    plt2 <- plt2 + ggplot2::theme_bw(base_size = 12)
    plt2 <- plt2 + ggplot2::geom_text(alpha = 0.75)
    # plt <- plt + ggplot2::theme(legend.position = "none")
    plt2 <- plt2 + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
    )
    plt2 <- plt2 + ggplot2::labs(
        x = "Method",
        y = "-log10(p-value)",
        fill = NULL
    )
    cowplot::save_plot(
        filename = paste(
            base, "all_results-facet.pdf", sep = "-"
        ),
        plot = plt2,
        nrow = 5,
        ncol = 5
    )

    # Per variable_and_effect, look at the correlation of -log10 pvalues
    # across methods
    pdf(
        file = paste(base, "model_comparison.pdf", sep = "-"),
        height = 8,
        width = 8
    )
    print(min(df_results$pvalue, na.rm = T))
    if (min(df_results$pvalue, na.rm = T) == 0) {
        filt <- df_results$pvalue == 0
        warning(paste0("Found", sum(filt), "cases where p-value == 0."))
        df_results$pvalue <- min(df_results$pvalue[~filt], na.rm = T)
    }
    df_results$pvalue_neg_log10 <- -log10(df_results$pvalue)
    for (var in unique(df_results$variable_and_effect)) {
        df_plt <- subset(df_results, variable_and_effect == var)
        ax_max <- max(df_plt$pvalue_neg_log10, na.rm = T) + 0.5
        df_plt <- reshape2::dcast(
            subset(df_results, variable_and_effect == var),
            cell_type ~ method_model,
            value.var="pvalue_neg_log10"
        )
        if (ncol(df_plt) != 2) {
            print(pairs(
                df_plt[, seq(ncol(df_plt)-1)+1],
                panel = function(x, y, ...) {
                    # points(x, y, ...);
                    text(x, y, df_plt[,"cell_type"])
                },
                main = paste(var, "-log10(p)")
            ))
            print(pairs(
                df_plt[, seq(ncol(df_plt)-1)+1],
                panel = function(x, y, ...) {
                    # points(x, y, ...);
                    text(x, y, df_plt[,"cell_type"])
                },
                main = paste(var, "-log10(p)"),
                xlim = c(0, ax_max),
                ylim = c(0, ax_max)
            ))
        }
    }
    dev.off()

    # For some reason cowplot::save_plot produces these plots
    if (file.exists("Rplots.pdf")) {
        file.remove("Rplots.pdf")
    }

    return(0)
}


main <- function() {
    # Run analysis
    run_time <- system.time(df_results <- command_line_interface())
    execution_summary <- paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours.\n"
    )
    cat(execution_summary)
    #return(execution_summary)
    #return(0) # For nextflow, more meaningful to return execution_summary
}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
