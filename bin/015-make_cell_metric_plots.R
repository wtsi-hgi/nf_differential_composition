#!/usr/bin/env Rscript

SCRIPT_NAME <- "make_cell_metric_plots.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyverse"))

set.seed(0)


#' Casts columns in a dataframe.
#'
#' \code{cast_variables} casts columns in a dataframe.
#'
#' @param df Data.frame.
#'     Row names of dataframe should correspond to sample names of used to
#'     identify each file (i.e., all(names(files) %in% rownames(metadata_df)))
#' @param verbose Logical.
#'     Write extra output to std.out. Default = TRUE.
#'
#' @return Data.frame.
#'     Updated Data.frame
#'
#' @export
cast_variables <- function(
        df,
        cols,
        cast_func,
        cast_func_description,
        verbose = T
    ) {

    if (verbose) {
        print(sprintf("Casting columns to be %s...", cast_func_description))
    }
    for (col in cols) {
        df[col] <- cast_func(df[[col]])
    }

    return(df)
}


make_plot_dataframe <- function(
    df_cell_count_metric,
    df_meta,
    variable_target,
    df_pvalues=NULL
) {
    # Make a long dataframe of with the following columns:
    # group
    # cell_type
    # value [counts or proportions]
    #
    # Assumes rownames of df_cell_count_metric and df_meta are experiment_id
    df_plot <- cbind(
        df_cell_count_metric,
        df_meta[rownames(df_cell_count_metric), variable_target]
    )
    # fix bug if only one variable
    v1 <- colnames(df_plot)
    colnames(df_plot) <- relist(
        replace(
            v1,
            v1=='df_meta[rownames(df_cell_count_metric), variable_target]',
            "group"
        ),
        skeleton = colnames(df_plot)
    )
    df_plot$experiment_id <- rownames(df_plot)
    df_plot <- reshape2::melt(
        df_plot,
        id.vars = c("experiment_id", "group"),
        variable.name = "cell_type",
        value.name = "value"
    )
    if (!is.null(df_pvalues)) {
        # Subset results dataframe down to variables we are interested in
        filt <- df_pvalues[["variable"]] == variable_target
        df_pvalues <- df_pvalues[filt,]
        if (nrow(df_pvalues) == 0) {
            warning(paste0(
                "Skipping pvalues for: ",
                variable_target,
                ". Likely due to this term being included from a forward model."
            ))
        } else {
            if (is.numeric(df_plot$group)) {
                df_plot <- merge(
                    df_plot,
                    df_pvalues,
                    by = c("cell_type"),
                    all.x = T,
                    all.y = F
                )
            } else {
                df_plot$key <- paste0(df_plot$cell_type, "::", df_plot$group)
                df_pvalues$key <- paste0(
                    df_pvalues$cell_type,
                    "::",
                    df_pvalues$variable_effect
                )
                df_pvalues$cell_type <- NULL
                df_plot <- merge(
                    df_plot,
                    df_pvalues,
                    by = c("key"),
                    all.x = T,
                    all.y = F
                )
                df_plot$key <- NULL
            }
        }
    }

    return(df_plot)
}


# Modified from SmillieCS-31348891
# https://github.com/cssmillie/ulcerative_colitis
make_barplot <- function(
        df_long,
        value="mean",
        error="se",
        xlab="",
        ylab="",
        title="",
        pvalue_column=NULL,
        coord_flip=FALSE,
        pos='dodge'
    ) {
    # At a minimum df_long needs the columns:  group  cell_type value

    # Value function
    if (value == 'mean') {
        vf <- mean
    } else if(value == 'median') {
        vf <- median
    } else {
        stop()
    }

    # Error function
    se <- function(x, na.rm=T) {
        sd(x, na.rm=na.rm) / sqrt(length(x))
    }
    if (error == 'sd') {
        ef <- sd
    } else if(error == 'se'){
        ef <- se
    } else {
        ef <- function(x, ...){0}
    }

    # Now aggregate that data and calculate the mean and sd across samples
    df_plot <- df_long %>%
        dplyr::group_by(group, cell_type) %>%
        dplyr::mutate(
            u = vf(value),
            s = ef(value)
        ) %>%
        dplyr::select(-c(experiment_id, value)) %>%
        unique() %>%
        dplyr::group_by(group, cell_type) %>%
        dplyr::mutate(
            n = n()
        ) %>%
        as.data.frame()
    if (max(df_plot$n) > 1) {
        stop("ERROR")
    }

    # Add in optional p-value labels
    if (!is.null(pvalue_column)) {
        if (!(pvalue_column %in% colnames(df_plot))) {
            stop("Missing pvalue column")
        }
        df_plot$lab1 <- ifelse(
            df_plot[[pvalue_column]] <= 0.01, '**',
            ifelse(df_plot[[pvalue_column]] <= 0.05, '*', '')
        )
    }

    # Plot data
    border <- NA
    if (pos == 'stack') {
        p <- ggplot2::ggplot(df_plot)
        p <- p + ggplot2::geom_bar(
            ggplot2::aes(x=cell_type, y=u, fill=group),
            colour=border,
            size=0.25,
            stat='identity'
        )
        if (error %in% c('sd', 'se')) {
            p <- p + ggplot2::geom_errorbar(
                ggplot2::aes(x=cell_type, ymin = u-s, ymax = u+s, fill=group),
                stat='identity',
                width=0.25
            )
        }
    } else {
        pos <- ggplot2::position_dodge(.9)
        p <- ggplot2::ggplot(df_plot)
        p <- p + ggplot2::geom_bar(
            ggplot2::aes(x=cell_type, y=u, fill=group),
            colour=border,
            size=.25,
            stat='identity',
            position=pos
        )
        if (error %in% c('sd', 'se')) {
            p <- p + ggplot2::geom_errorbar(
                ggplot2::aes(x=cell_type, ymin=u-s, ymax=u+s, fill=group),
                stat='identity',
                position=pos,
                width=.25
            )
        }
    }

    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::scale_fill_brewer(palette = "Dark2", name = NULL)
    p <- p + ggplot2::labs(
        x = xlab,
        y = ylab,
        title = title
    )

    # Facet wrap
    # if (do.facet == TRUE) {
    #     p <- p + ggplot2::facet_grid(group ~ ., scales='free')
    # }

    dy <- max(df_plot$u + df_plot$s, na.rm=T) * 0.01
    if (coord_flip == FALSE) {
        p <- p + ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
        )
        if(!is.null(pvalue_column)) {
            p <- p + ggplot2::geom_text(
                ggplot2::aes(x=cell_type, y=u+s+dy, label=lab1, group=group),
                hjust='center',
                vjust=0,
                size=5,
                angle=0,
                position=pos
            )
        }
    } else {
        p <- p + ggplot2::coord_flip()
        if (!is.null(pvalue_column)) {
            p <- p + ggplot2::geom_text(
                ggplot2::aes(x=cell_type, y=u+s+dy, label=lab1, group=group),
                hjust='center',
                vjust=1,
                size=5,
                angle=90,
                position=pos
            )
        }
    }

    return(p)
}


make_cell_count_plot <- function(df_plot, var) {
    plt <- ggplot2::ggplot(df_plot, ggplot2::aes(
        x = group,
        y = value
    ))
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::labs(
        x = var,
        y = "Number of cells",
        fill = NULL
    )
    plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
    if (is.numeric(df_plot$group)) {
        plt <- plt + ggplot2::geom_point(alpha = 0.5)
        plt <- plt + ggplot2::geom_smooth(
            #method = "lm",
            se = FALSE
        )
    } else {
        plt <- plt + ggplot2::geom_boxplot(outlier.shape = NA,
            ggplot2::aes(fill = group)
        )
        plt <- plt + ggplot2::geom_jitter(alpha = 0.5)
        plt <- plt + ggplot2::scale_fill_brewer(palette = "Dark2")
    }
    plt <- plt + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
    )
    return(plt)
}


make_cell_prop_plot <- function(df_plot, var) {
    plt <- ggplot2::ggplot(df_plot, ggplot2::aes(
        x = group,
        y = value
    ))
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::labs(
        x = var,
        y = "Proportion of cells",
        fill = NULL
    )
    plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
    if (is.numeric(df_plot$group)) {
        plt <- plt + ggplot2::geom_point(alpha = 0.5)
        plt <- plt + ggplot2::geom_smooth(
            #method = "lm",
            se = FALSE
        )
    } else {
        plt <- plt + ggplot2::geom_boxplot(outlier.shape = NA,
            ggplot2::aes(fill = group)
        )
        plt <- plt + ggplot2::geom_jitter(alpha = 0.5)
        plt <- plt + ggplot2::scale_fill_brewer(palette = "Dark2")
    }
    plt <- plt + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
    )
}


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--counts_file"),
            type = "character",
            help = paste0(
                "Tab-delimited input file cell type count per sample.",
                " First column = index (i.e., sample / experiment_id)",
                " all other columns = cell types."
            )
        ),

        optparse::make_option(c("--metadata_file"),
            type = "character",
            help = paste0(
              "Tab-delimited metadata file.",
              " First column = index (i.e., sample / experiment_id)",
              " all other columns = cell types."
            )
        ),

        optparse::make_option(c("--results_file"),
            type = "character",
            default = "",
            help = paste0(
              "Tab-delimited results file. All variables in the formula will ",
              " be plotted.",
              " Must contain the following columns:",
              " cell_label,variable_effect,pvalue,formula."
            )
        ),

        optparse::make_option(c("--variables_to_plot"),
            type = "character",
            default = "",
            help = paste0(
                "Comma separated list of variables to plot (optional).",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--continuous_variables"),
            type = "character",
            default = "",
            help = paste0(
                "Comma seperated list of variables to cast as numeric.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--discrete_variables"),
            type = "character",
            default = "",
            help = paste0(
                "Comma seperated list of variables to cast as character.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--discrete_levels"),
            type = "character",
            default = "",
            help = paste0(
                "Order of levels to cast discrete variables. The first",
                " variable will be reference. Example:",
                " smoking_status::no,yes,ex-smoker;;sex::F,M.",
                " In the example, no will be refence for smoking_status",
                " variable F will be reference for sex variable.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--model_plots_only"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Only print plots where p-values from model are incorporated.",
                " [default: %default]"
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

        # optparse::make_option(c("--verbose"),
        #     type = "logical",
        #     action = "store_true",
        #     default = FALSE,
        #     help = paste0(
        #         "Verbose mode (write extra info to std.err).",
        #         " [default: %default]"
        #     )
        # )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Plots results from differential celltype composition analysis."
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
    f_counts <- param[["counts_file"]]
    f_meta <- param[["metadata_file"]]
    f_results <- param[["results_file"]]
    model_plots_only <- param[["model_plots_only"]]
    variables_plot_str <- param[["variables_to_plot"]]
    continuous_variables_str <- param[["continuous_variables"]]
    discrete_variables_str <- param[["discrete_variables"]]
    levels_str <- param[["discrete_levels"]]
    base <- param[["out_file"]]

    # f_counts <- "diff_comp-counts.tsv.gz"
    # f_meta <- "diff_comp-metadata.tsv.gz"
    # formula_rhs <- "~ disease_status"
    # variable_target <- "disease_status"
    # forward_selection_variables_str <- "age,sex,I(age^2),smoker_at_time_of_biopsy"
    # continuous_variables_str <- "age,I(age^2)"
    # discrete_variables_str <- "disease_status,sex,smoker_at_time_of_biopsy"
    # levels_str <- "disease_status::healthy,cd"

    # Update the variables
    variables_plot <- c()
    if (variables_plot_str != "") {
        variables_plot <- strsplit(
            x = variables_plot_str,
            split = ",",
            fixed = TRUE
        )[[1]]
    }

    # Get variables to have their levels set.
    variables_and_levels <- c()
    if (levels_str != "") {
        vars_and_levels_str <- unlist(strsplit(levels_str, ";;"))
        for (i in vars_and_levels_str) {
            var <- strsplit(i, "::")[[1]][1]
            lvls <- unlist(strsplit(strsplit(i, "::")[[1]][2], ","))
            variables_and_levels[[var]] <- lvls
        }
    }

    # Read in the counts_file.
    df_counts <- read.csv(
        f_counts,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
    )
    # mtx_counts <- as.matrix(df_counts)

    # Read in the metadata file.
    df_meta <- read.csv(
        f_meta,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
        #row.names = 1
    )
    rownames(df_meta) <- df_meta[[colnames(df_meta)[1]]]
    df_meta <- df_meta[rownames(df_counts),]

    # Cast discrete_vars first so levels are not overwritten
    if (discrete_variables_str != "") {
        discrete_vars <- strsplit(
            x = discrete_variables_str,
            split = ",",
            fixed = TRUE
        )[[1]]
        df_meta <- cast_variables(
            df_meta,
            discrete_vars,
            as.character,
            "characters",
            verbose = FALSE
        )
    }

    # Get variables to set as factors
    variables_and_levels <- c()
    if (levels_str != "") {
        vars_and_levels_str <- unlist(strsplit(levels_str, ";;"))
        for (i in vars_and_levels_str) {
            var <- strsplit(i, "::")[[1]][1]
            lvls <- unlist(strsplit(strsplit(i, "::")[[1]][2], ","))
            variables_and_levels[[var]] <- lvls
        }
    }

    # Set the factor levels for variables
    for (var in names(variables_and_levels)) {
        # Check to make sure these arrays are identical
        if (!identical(
            sort(variables_and_levels[[var]]),
            sort(unique(df_meta[[var]]))
        )) {
            stop(cat(
                "Missing variable level.\nvariables_and_levels:\t",
                paste(sort(variables_and_levels[[var]]), collapse = ","),
                "\ndf_meta:\t",
                paste(sort(unique(df_meta[[var]])), collapse = ","),
                "\n"
            ))
        }
        df_meta[[var]] <- factor(
            as.character(df_meta[[var]]),
            levels = variables_and_levels[[var]],
            labels = variables_and_levels[[var]]
        )
    }

    # Cast any variables that need to be cast
    # Do this after setting levels for anything so one can control how
    # a factor is cast to numeric.
    if (continuous_variables_str != "") {
        continuous_vars <- strsplit(
            x = continuous_variables_str,
            split = ",",
            fixed = TRUE
        )[[1]]
        df_meta <- cast_variables(
            df_meta,
            continuous_vars,
            as.numeric,
            "numeric",
            verbose = FALSE
        )
        # Mean impute missing variables
        for (i in continuous_vars) {
            filt <- is.na(df_meta[[i]])
            if (any(filt)) {
                cat(paste0("Mean imputing:\t", i, "\n"))
                df_meta[filt, i] <- mean(df_meta[[i]], na.rm = TRUE)
            }
        }
    }

    # Read in the model results
    if (f_results != "") {
        df_results <- read.csv(
            f_results,
            sep = "\t",
            header = TRUE,
            stringsAsFactors = FALSE
            #row.names = 1
        )
        # Now use the formula column to get variables
        # NOTE: could be more than one fomula because of step forward
        if ("formula" %in% colnames(df_results)) {
            formulas <- unique(df_results$formula)
            for (i in formulas) {
                cat("Getting terms for plotting: ", i, "\n")
                variables_plot <- unique(c(
                    attr(terms(as.formula(i)), "term.labels"),
                    variables_plot
                ))
            }
        }
        # Use variable column to also get any values (e.g., if this results)
        # df did not have a formula column
        if ("variable" %in% colnames(df_results)) {
            vars <- unique(df_results$variable)
            filt <- vars != "intercept"
            variables_plot <- unique(c(vars[filt], variables_plot))
        }
    } else {
        df_results <- NULL
    }

    # Evaluate any terms that need to be evaluated on original scale
    for (i in variables_plot) {
        # If the term is one that we need to evaluate and generate, then do so
        if (grepl("I(", i, fixed = TRUE)) {
            cat(paste0("Evaluating:\t", i, "\n"))
            df_meta[[i]] <- eval(
                parse(text = i),
                df_meta
            )
        }
    }

    # Prep
    df_results_list <- list(NULL) # list of results grouped by formula
    if (!is.null(df_results)) {
        if ("formula" %in% colnames(df_results)) {
            formulas <- unique(df_results$formula)
            j <- 1
            for (i in formulas) {
                df_results_list[[j]] <- subset(df_results, formula == i)
                j <- j + 1
            }
        }
    }
    df_proportions <- as.matrix(df_counts)
    df_proportions <- as.data.frame(
        df_proportions / rowSums(df_proportions)
    )

    for (df_results_i in df_results_list) {
        for (var in variables_plot) {
            print(var)
            var_clean <- var
            var_clean <- gsub("I\\(", "", var_clean)
            var_clean <- gsub("\\)", "", var_clean)
            var_clean <- gsub("\\^", "power", var_clean)

            formula_str <- unique(df_results_i$formula)
            formula_str <- gsub("\\+", "", formula_str)
            formula_str <- gsub("\\~", "", formula_str)
            formula_str_short <- substr(formula_str, 1, 50)

            # Make plots of counts #############################################
            df_plot <- make_plot_dataframe(
                df_cell_count_metric = df_counts,
                df_meta = df_meta,
                variable_target = var,
                df_pvalues = df_results_i
            )

            pvalue_column <- NULL
            if ("qvalue_bh" %in% colnames(df_plot)) {
                pvalue_column <- "qvalue_bh"
            } else if("pvalue" %in% colnames(df_plot)) {
                pvalue_column <- "pvalue"
            }

            title <- ""
            if ("formula" %in% colnames(df_plot)) {
                title <- unique(df_plot$formula)
            }

            if (!is.numeric(df_plot$group)) {
                plt <- make_barplot(
                    df_long=df_plot,
                    value="mean",
                    error="se",
                    xlab="",
                    ylab="Number of cells",
                    title=title,
                    pvalue_column=pvalue_column,
                    coord_flip=FALSE
                )
                cowplot::save_plot(
                    plot = plt,
                    filename = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_count_barplot-mean.pdf",
                        sep = "-"
                    ),
                    nrow = 1,
                    ncol = 2.5
                )
                plt <- make_barplot(
                    df_long=df_plot,
                    value="median",
                    error="se",
                    xlab="",
                    ylab="Number of cells",
                    title=title,
                    pvalue_column=pvalue_column,
                    coord_flip=FALSE
                )
                cowplot::save_plot(
                    plot = plt,
                    filename = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_count_barplot-median.pdf",
                        sep = "-"
                    ),
                    nrow = 1,
                    ncol = 2.5
                )
            }

            if (!model_plots_only) {
                pdf(
                    file = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_count_plots.pdf",
                        sep = "-"
                    ),
                    height = 10,
                    width = 10
                )
                plt <- make_cell_count_plot(df_plot, var)
                print(plt)
                plt <- plt + ggplot2::ylim(
                    0,
                    max(df_plot$value, na.rm = TRUE) + 5
                )
                print(plt)
                dev.off()

                # Make the same plot but only for hits
                if (!is.null(pvalue_column)) {
                    df_tmp <- subset(df_plot, df_plot[[pvalue_column]] < 0.05)
                    if (nrow(df_tmp) > 0) {
                        pdf(
                            file = paste(
                                base,
                                formula_str_short,
                                var_clean,
                                paste0(pvalue_column, "_lessthan0pt05"),
                                "cell_count_plots.pdf",
                                sep = "-"
                            ),
                            height = 10,
                            width = 10
                        )
                        plt <- make_cell_count_plot(df_tmp, var)
                        print(plt)
                        plt <- plt + ggplot2::ylim(
                            0,
                            max(df_plot$value, na.rm = TRUE) + 5
                        )
                        print(plt)
                        dev.off()
                    }
                }
            }
            ####################################################################

            # Make plots of propotions #########################################
            df_plot <- make_plot_dataframe(
                df_cell_count_metric = df_proportions,
                df_meta = df_meta,
                variable_target = var,
                df_pvalues = df_results_i
            )
            if (!is.numeric(df_plot$group)) {
                plt <- make_barplot(
                    df_long = df_plot,
                    value = "mean",
                    error = "se",
                    xlab = "",
                    ylab = "Proportion of cells",
                    title=title,
                    pvalue_column = pvalue_column,
                    coord_flip = FALSE
                )
                cowplot::save_plot(
                    plot = plt,
                    filename = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_proportion_barplot-mean.pdf",
                        sep = "-"
                    ),
                    nrow = 1,
                    ncol = 2.5
                )
                plt <- make_barplot(
                    df_long = df_plot,
                    value = "median",
                    error = "se",
                    xlab = "",
                    ylab = "Proportion of cells",
                    title=title,
                    pvalue_column = pvalue_column,
                    coord_flip = FALSE
                )
                cowplot::save_plot(
                    plot = plt,
                    filename = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_proportion_barplot-median.pdf",
                        sep = "-"
                    ),
                    nrow = 1,
                    ncol = 2.5
                )
            }

            if (!model_plots_only) {
                pdf(
                    file = paste(
                        base,
                        formula_str_short,
                        var_clean,
                        "cell_proportion_plots.pdf",
                        sep = "-"
                    ),
                    height = 10,
                    width = 10
                )
                plt <- make_cell_prop_plot(df_plot, var)
                print(plt)
                plt <- plt + ggplot2::ylim(
                    0,
                    max(df_plot$value, na.rm = TRUE)
                )
                print(plt)
                dev.off()

                # Make the same plot but only for hits
                if (!is.null(pvalue_column)) {
                    df_tmp <- subset(df_plot, df_plot[[pvalue_column]] < 0.05)
                    if (nrow(df_tmp) > 0) {
                        pdf(
                            file = paste(
                                base,
                                formula_str_short,
                                var_clean,
                                paste0(pvalue_column, "_lessthan0pt05"),
                                "cell_proportion_plots.pdf",
                                sep = "-"
                            ),
                            height = 10,
                            width = 10
                        )
                        plt <- make_cell_prop_plot(df_tmp, var)
                        print(plt)
                        plt <- plt + ggplot2::ylim(
                            0,
                            max(df_plot$value, na.rm = TRUE)
                        )
                        print(plt)
                        dev.off()
                    }
                }
            }
            ###################################################################
        }
    }

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
