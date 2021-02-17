#!/usr/bin/env Rscript

SCRIPT_NAME <- "run_mannwhitneyu_test.R"

suppressPackageStartupMessages(library("optparse"))

set.seed(0)


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

        optparse::make_option(c("--condition"),
            type = "character",
            default = "",
            help = paste0(
                "Condition column. If categorical >2 levels, first is",
                " selected as a reference.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--levels"),
            type = "character",
            default = "",
            help = paste0(
                "Order of levels to cast catagorical variables. The first",
                " variable will be reference. Example:",
                " smoking_status::no,yes,ex-smoker;;sex::F,M.",
                " In the example, no will be refence for smoking_status",
                " variable F will be reference for sex variable.",
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
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Performs differential celltype composition analysis."
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
    condition <- param[["condition"]]
    levels_str <- param[["levels"]]
    base <- param[["out_file"]]

    # Get variables to cast
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
    for (var in condition) {
        if (var %in% names(variables_and_levels)) {
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
                df_meta[[var]],
                levels = variables_and_levels[[var]],
                labels = variables_and_levels[[var]]
            )
        }
    }
    df_meta$condition_dummy_var_032352 <- df_meta[[condition]]


    # Get all possible condition pairing using the first variable as a
    # reference
    if (!is.numeric(df_meta[["condition_dummy_var_032352"]])) {
        tmp_list <- levels(df_meta[["condition_dummy_var_032352"]])
        if (length(tmp_list) == 0) {
            tmp_list <- unique(df_meta[["condition_dummy_var_032352"]])
        }
        condition_parings <- data.frame("condition" = tmp_list[-1])
        condition_parings$reference <- tmp_list[1]
    } else {
        print("WARNING condition is numeric, not running MWU test.")
        quit() # Don't use stop because raise error
    }


    # Calculate the proportion of cells for each sample
    df_proportions <- as.matrix(df_counts)
    df_proportions <- as.data.frame(df_proportions / rowSums(df_proportions))

    list_of_results <- list()
    i <- 1
    for (group_name in colnames(df_counts)) {
        for (row in 1:nrow(condition_parings)) {
            x_var <- condition_parings[row, "condition"]
            y_var <- condition_parings[row, "reference"]

            # Mann Whiney U test on proportions
            df_tmp <- cbind(
                df_meta[rownames(df_counts),],
                df_proportions[group_name]
            )
            df_tmp$gn <- df_tmp[[group_name]]

            # Init null results in case test fails because not enough case
            # controls
            list_of_results[[i]] <- data.frame(
                "cell_type" = group_name,
                "statistic" = NA,
                "pvalue" = NA,
                "estimate" = NA,
                "conf_int_95" = NA,
                "method" = "mann_whitney_u",
                "data_type" = "fraction_of_cells",
                "condition" = condition,
                "condition_effect" = x_var
            )

            foo <- tryCatch({
                # The unpaired two-samples Wilcoxon test (also known as Wilcoxon
                # rank sum test or Mann-Whitney test)
                res <- wilcox.test(
                    subset(df_tmp, condition_dummy_var_032352 == x_var)$gn,
                    subset(df_tmp, condition_dummy_var_032352 == y_var)$gn,
                    alternative = "two.sided",
                    exact = TRUE,
                    conf.int = TRUE
                )
                # NOTE: Does not estimate the difference in medians (a common
                # misconception) but rather the median of the difference between
                # a sample from x and a sample from y.
                # NOTE: condition_effect = condition that is higher if estimate
                # pos.
                list_of_results[[i]] <- data.frame(
                    "cell_type" = group_name,
                    "statistic" = res$statistic,
                    "pvalue" = res$p.value,
                    "estimate" = res$estimate,
                    "conf_int_95" = paste(res$conf.int, collapse = ","),
                    "method" = "mann_whitney_u",
                    #"data_type" = "proportion_of_cells",
                    "condition" = condition,
                    "condition_effect" = x_var
                )
            }, error = function(e) {
                print(e)
            })
            i <- i + 1
        }
    }
    df_tmp <- as.data.frame(do.call(rbind, list_of_results))
    df_tmp$qvalue_bh <- p.adjust(
        df_tmp$pvalue,
        method = "BH"
    )

    # Make the final results dataframe and sort it
    df_res_univariate <- df_tmp
    df_res_univariate <- df_res_univariate[
        with(df_res_univariate, order(method, pvalue, cell_type)),
    ]

    # Save the results
    gzfh <- gzfile(
        paste0(base, "-mannwhitneyu_results.tsv.gz"),
        "w",
        compression = 9
    )
    write.table(
        df_res_univariate,
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
    close(gzfh)


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
