#!/usr/bin/env Rscript

SCRIPT_NAME <- "run_dirichlet_regression.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DirichletReg"))
#suppressPackageStartupMessages(library("rgl")) # OPTOINAL: DirichletReg plots
# Another possible package is Compositional::diri.reg
# suppressPackageStartupMessages(library("Compositional"))

set.seed(0)

# set setting, but re-enable upon exit
old <- options(
    stringsAsFactors = FALSE,
    try.outFile = stdout(),
    show.error.messages = TRUE
)
on.exit(options(old), add = TRUE)


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


extractAIC.DirichletRegModel <- function(fit, scale, k = 2, ...) {
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
    # n <- nobs(fit)
    # npar <- fit$npar
    # dev <- -2*fit$logLik
    # c(npar, dev + k * npar)
}


#' Performs stepwise forward model selection.
#'
#' \code{step_forward_dirichreg} iteratively adds terms to DirichReg model
#'
#' @param model DirichReg.model.
#'     Initial DirichReg.model
#' @param other_variables List.
#'     List of other variables to consider adding to model.
#' @param data DirichletReg::DR_data.
#'     Model data.
#' @param optimization_method Character.
#'     What optimization method to use? LRT = likelihood-ratio test. Currently
#'     the only supported option is LRT.
#' @param optimization_var_threshold Numeric.
#'     If optimization value is >= this variable, then the model is done.
#' @param verbose Logical.
#'     Write extra output to std.out. Default = TRUE.
#'
#' @return DirichReg.model.
#'     Final model.
#'
#' @importFrom DirichletReg DirichReg
#' @export
step_forward_dirichreg <- function(
        model,
        other_variables,
        data,
        optimization_method = "LRT",
        optimization_var_threshold = 0.05,
        verbose = FALSE
    ) {

    if (verbose) {
        cat("\n\ncalling step_forward_dirichreg\n")
    }

    # TODO: add option for other methods: AIC(model), BIC(model)
    if (optimization_method != "LRT") {
        stop("Invalid optimization_method.")
    }

    # Get the variables in our current model
    variables_cur <- attr(terms(model$formula), "term.labels")

    # Init ll for less complex model to compare to new model
    ll_fit_h0 <- logLik(model)

    # Fit new models
    new_model <- model
    optimization_var <- 1
    for (var in other_variables) {
        if (verbose) {
            cat("[", var,
                "]: testing this term.\n"
            )
        }
        # If the variable is already in model, then skip it
        if (var %in% variables_cur) {
            if (verbose) {
                cat("[", var,
                    "]: skipping this term because already in model.\n"
                )
            }
            next
        }
        # Fit the model adding var to our list of covariates
        variables_i <- c(variables_cur, var)
        formula_rhs <- paste0(variables_i, collapse = " + ")
        model_i <- do.call(DirichletReg::DirichReg, args = list(
            formula = as.formula(paste("Y", formula_rhs, sep = "~")),
            data = data,
            model = "common",
            verbosity = 0,
            control = list(maxit = 50000, iterlim = 50000)
        ))
        encountered_error <- TRUE
        # This catches collinearity issues.
        # Strange though that after the catch, none of the std out of
        # subsequent code is written.
        # NOTE: This does not catch the case of perfect collinarity
        try(
            expr = {
                sink("/dev/null")
                tmp <- summary(model_i)
                sink()
                encountered_error <- FALSE
            },
            silent = TRUE
        ) # end try catch
        if (encountered_error) {
            #if (verbose) {
                cat("[", var,
                    "]: adding this term causes model not to converge.\n"
                )
            #}
            next
        }
        # Calculate our optimization_var, comparing this model to original
        # model
        optimization_var_i <- NULL
        if (optimization_method == "LRT") {
            # Our h0 model is nested in this model since we added var to
            # variables_cur
            ll_fit_h1 <- logLik(model_i)
            teststat <- -2 * (as.numeric(ll_fit_h0) - as.numeric(ll_fit_h1))
            optimization_var_i <- pchisq(
                teststat,
                df = (attr(ll_fit_h1, "df") - attr(ll_fit_h0, "df")),
                lower.tail = FALSE
            )
        } # end optimization_method == "LRT"
        if (optimization_var_i < optimization_var) {
            if (verbose) {
                cat("[", var, "]: adding this term improves the model.\n")
                cat("\t[ formula ]:\tY ~ ",
                    formula_rhs,
                    "\n"
                )
                cat("\t[ optimization_var ]:\t",
                    optimization_var_i,
                    "(optimization_var_i) <",
                    optimization_var,
                    "(optimization_var)",
                    "\n"
                )
            }
            new_model <- model_i
            optimization_var <- optimization_var_i
            #new_variables <- variables_i
        } else {
            if (verbose) {
                cat("[", var,
                    "]: adding this term does not improve the model.\n"
                )
                cat("\t[ formula ]:\tY ~",
                    formula_rhs,
                    "\n"
                )
                cat("\t[ optimization_var ]:\t",
                    optimization_var_i,
                    "(optimization_var_i) >=",
                    optimization_var,
                    "(optimization_var)",
                    "\n"
                )
            }
        } # end optimization_var_i < optimization_var
    } # end model for loop

    if (optimization_var < optimization_var_threshold) {
        # Get the variables in the new model
        new_variables <- attr(terms(new_model$formula), "term.labels")
        # Calculate the variables for the next call
        other_variables <- other_variables[
            !(other_variables %in% new_variables)
        ]
        if (verbose) {
            cat("UPDATED model:",
                "\n\t[ formula ]:\t",
                paste0("Y ~ ", paste0(new_variables, collapse = " + ")),
                "\n\t[ remaining_vars_to_add ]:\t",
                paste0(other_variables, collapse = ","),
                "\n\t[ optimization_var ]:\t",
                optimization_var,
                "\n"
            )
        } # end verbose
        new_model <- step_forward_dirichreg(
            new_model,
            other_variables,
            data,
            optimization_method,
            optimization_var_threshold,
            verbose
        )
    } else {
        if (verbose) {
            cat("\n\nNo model passed threshold. Returning original model.\n")
            cat("FINAL model:",
                "\n\t[ formula ]:\t",
                paste0("Y ~ ", paste0(variables_cur, collapse = " + ")),
                "\n\t[ remaining_vars_not_added ]:\t",
                paste0(other_variables, collapse = ","),
                "\n\t[ optimization_var ]:\t",
                optimization_var,
                "(optimization_var_i) >",
                optimization_var_threshold,
                "(optimization_var_threshold)",
                "\n"
            )
        } # end verbose
        return(model)
    } # end optimization_var < optimization_var_threshold
}


#' Performs Dirichlet regression.
#'
#' \code{dirichlet_regression} performs Dirichlet regression on single cell
#' count data, similar to SmillieCS-31348891
#' (https://github.com/cssmillie/ulcerative_colitis). Note that for modelling,
#' the cell count data is transformed to proportions.
#'
#' @param df_counts Data.frame.
#'     Data.frame of counts where row names are samples and column names
#'     are celltypes.
#' @param df_independent_vars Data.frame.
#'     Data.frame of variables that maybe included in the model where
#'     row names are samples and columns are variables. NOTE: row order should
#'     match that of \code{df_counts}.
#' @param formula_rhs Character.
#'     Right hand side of a formula declaration. Example: "~ sex + age".
#' @param variable_target Character.
#'     Target variable that we want to test. If provided, a reduced model
#'     where this variable is removed will be generated and tested against
#'     the full model using the likelihood-ratio test. Example: "sex".
#' @param fit_ols Logical.
#'     Fit OLS as proposed by Aitchison where \code{y = log(y/y_reference)}.
#' @param forward_selection_variables List.
#'     If provided, will perform stepwise forward model selection attempting
#'     to add the terms listed here to the model specifed in
#'     \code{formula_rhs}.
#' @param verbose Logical.
#'     Write extra output to std.out. Default = TRUE.
#'
#' @return List.
#'     A list with three entries:
#'          1. model:       The final DirichReg.model.
#'          2. df_results:  A clean dataframe of the results.
#'          3. dr_data:     The DirichletReg::DR_data object.
#'
#' @importFrom DirichletReg DR_data
#' @importFrom DirichletReg DirichReg
#' @export
dirichlet_regression <- function(
        df_counts,
        df_independent_vars,
        formula_rhs,
        variable_target = NULL,
        fit_ols = FALSE,
        forward_selection_variables = c(),
        optimization_method = "LRT",
        optimization_var_threshold = 0.05,
        verbose = TRUE
    ){

    # Ensure the sample ids across rows match
    df_independent_vars <- df_independent_vars[rownames(df_counts),]

    # Get formula values
    formula_rhs <- gsub("[[:space:]]", "", formula_rhs) # strip whitespace
    formula_rhs <- gsub("~", "", formula_rhs)
    formula_str <- paste("Y", formula_rhs, sep = "~")
    formula_obj <- as.formula(formula_str)  # , env = environment()
    independent_vars <- attr(terms(formula_obj), "term.labels")
    cell_types <- colnames(df_counts)

    # Calculate regression
    if (class(df_counts) != "data.frame") {
        df_counts <- as.data.frame(df_counts)
    }
    # If the alternative parametrization will be used, the base (i.e., omitted)
    # variable, which is the first by default, can be set to be any other of
    # the variables. The transformed variables must be re-integrated into the
    # original data frame with a new, unique name
    # (e.g., DF$Y <- DR_data(DF[, 1:3]); see examples below) which is then
    # passed on the the fitting function.
    sink("/dev/null")
    df_counts$Y <- DirichletReg::DR_data(df_counts) # Y is now propotions
    sink()
    data <- cbind(df_counts, df_independent_vars)

    # NOTE: There are two parametrizations as described below
    # Maier, M. DirichletReg: Dirichlet Regression for Compositional Data in R.
    #     Research Report Series 125, (2014).
    # (a) A common approach: the Dirichlet distribution’s alpha parameters are
    #     directly modeled by covariates using a log-link (i.e., all alpha
    #     parameters are modeled independently).
    # (b) An alternative approach: where the Dirichlet distribution is
    #     reparametrized to allow for modeling ‘means’ and ‘dispersion’.
    #     Heteroscedasticity is therefore taken in account implicitly in the
    #     common parametrization and explicitly in the alternative
    #     parametrization (i.e., means and precision are modeled).
    #
    #     Like in betareg, the alternative parametrization is defined using
    #     two blocks of variables. The first gives the dependent variable and
    #     predictors for the means while the latter specifies the precision
    #     model to set up the according predictors. Therefore
    #     y ~ x1 * x2 | x1 + x2 sets up a model where the main effects and
    #     the interaction of x1 and x2 are used to predict the mean,
    #     but the precision is modeled merely by their main effects.
    #
    # NOTE: weights parameter can be used to weight observations (by
    # probability or frequency).
    #
    # NOTE: If alternative approach, one of the Y variables must be made a
    # reference category, meaning no results for that cell type.
    model_0 <- do.call(DirichletReg::DirichReg, args = list(
        formula = formula_obj,
        data = data,
        model = "common",
        verbosity = 0,
        control = list(maxit = 50000, iterlim = 50000)
    ))
    if (length(forward_selection_variables) > 0) {
        # # An aproach using step in R
        # tmp_vars <- unique(c(independent_vars, forward_selection_variables))
        # tmp_formula <- paste0("Y ~ ", paste0(tmp_vars, collapse = " + "))
        # # NOTE: In theory this should work but it fails due to error in silly
        # # Formula parser:
        # # Error in Formula(new) : inherits(object, "formula") is not TRUE
        # model_full <- do.call(MASS::stepAIC, args = list(
        #     object = model_0,
        #     scope = list(
        #         upper = as.formula(
        #             paste0("~", paste0(tmp_vars, collapse = "+"))
        #         ),
        #         lower = as.formula(
        #             paste0("~", paste0(independent_vars, collapse = "+"))
        #         )
        #     ),
        #     direction = 'forward'
        # ))
        # A custom approach using LRT
        model_full <- step_forward_dirichreg(
            model_0,
            forward_selection_variables,
            data,
            optimization_method = optimization_method,
            optimization_var_threshold = optimization_var_threshold,
            verbose = verbose
        )
        # Update variables based on the new full model:
        independent_vars <- attr(terms(model_full$formula), "term.labels")
        formula_rhs <- paste(independent_vars, collapse = " + ")
    } else {
        model_full <- model_0
    }
    # As described in Maier white paper (see above), we can also use the
    # OLS-regression applied to log-ratio transformed data. This is likely to
    # show substantial effects due to heteroscedasticity.
    if (fit_ols) {
        data_ols <- df_independent_vars
        data_ols$Y <- log(df_counts$Y / (1 - df_counts$Y))
        model_full_ols <- lm(
            as.formula(paste("Y", formula_rhs, sep = "~")),
            data_ols
        )
        model_full$model_ols_logratios <- model_full_ols
    }
    # If the variable_target != NULL, we will compare the formula with and
    # without the variable_target term to see if adding this term better
    # models the data.
    if (!is.null(variable_target)) {
        independent_vars_reduced <- independent_vars[
            independent_vars != variable_target
        ]
        formula_rhs_reduced <- paste0(
            independent_vars_reduced,
            collapse = " + "
        )
        if (formula_rhs_reduced == "") {
            formula_rhs_reduced = "1"
        }
        cat("Reduced formula:\t Y ~", formula_rhs_reduced, "\n")
        tryCatch({
            model_reduced <- do.call(DirichletReg::DirichReg, args = list(
                formula = as.formula(
                    paste("Y", formula_rhs_reduced, sep = "~")
                ),
                data = data,
                model = "common",
                verbosity = 0,
                control = list(maxit = 50000, iterlim = 50000)
            ))
            model_full$model_reduced <- model_reduced

            # Compare our two models using anova as described in:
            # Maier, M. DirichletReg: Dirichlet Regression for Compositional
            #     Data in R. Research Report Series 125, (2014).
            model_full$model_reduced_anova <- anova(model_reduced, model_full)
            if (
                min(
                    model_full$model_reduced_anova[['Pr(>Chi)']], na.rm = T
                ) > 0.05
            ) {
                warning(paste0(
                    "Full model not better fit than reduced: Y ~ ",
                    formula_rhs_reduced
                ))
            }
            # NOTE: I verified that in some cases the LRT p-value may be 1 and
            # anova very small. Not sure why this happens.
            #
            # We can also do this manually with the likelihood-ratio-test
            # (e.g., via lmtest::lrtest)
            # Here is a good intro https://api.rpubs.com/tomanderson_34/lrt
            #
            # Fnd the loglikelihoods of each model.
            ll_fit_h0 <- logLik(model_reduced)
            ll_fit_h1 <- logLik(model_full)
            # Calculate the test statistic
            teststat <- -2 * (as.numeric(ll_fit_h0)-as.numeric(ll_fit_h1))
            # The test statistic follows a chi-squared distribution with
            # degrees of freedom equal to the difference in the number of
            # free parameters between the reduced model and the full model
            model_full$model_reduced_lrt_pvalue <- pchisq(
                teststat,
                df = (attr(ll_fit_h1, "df") - attr(ll_fit_h0, "df")),
                lower.tail = FALSE
            )
            if (model_full$model_reduced_lrt_pvalue > 0.05) {
                warning(paste0(
                    "Full model not better fit than reduced: Y ~ ",
                    formula_rhs_reduced
                ))
            }
        }, error = function(e) {
            warning(paste0(
                "Reduced model did not converge: Y ~ ",
                formula_rhs_reduced
            ))
        }) # end try catch
    }

    # Compile results of the full model
    sink("/dev/null")
    model_full_results <- summary(model_full)
    sink()
    # Clean up the columns of the coefficient matrix
    res_sum <- list(
        "Estimate" = "estimate",
        "Std. Error" = "std_err",
        "z value" = "zscore",
        "Pr(>|z|)" = "pvalue"
    )
    for (col in colnames(model_full_results$coef.mat)) {
        if (!(col %in% names(res_sum))) {
            res_sum[[col]] <- col
        }
    }
    df_results <- data.frame(model_full_results$coef.mat)
    rownames(df_results) <- NULL
    colnames(df_results) <- unlist(res_sum[
        colnames(model_full_results$coef.mat)
    ])
    df_results$variable <- rownames(model_full_results$coef.mat)
    # Get cell_type:variable
    df_results$cell_type <- names(model_full_results$coefficients)
    # Strip variable from cell type column
    df_results$cell_type <- unlist(apply(
        df_results[, c("variable", "cell_type")],
        1,
        FUN = function(row) {
            new_str <- gsub(
                paste0(":", row["variable"]),
                "",
                row["cell_type"],
                fixed = TRUE
            )
            new_str <- gsub(":\\(Intercept\\)", "", new_str)
            return(new_str)
        }
    ))
    # Clean up the variable names
    df_results$variable <- unlist(lapply(
        df_results$variable,
        FUN = function(x) {
            return(gsub("\\(Intercept\\)", "intercept", x))
        }
    ))
    # Get the effect variable for categorical variables
    df_results$variable_effect <- unlist(lapply(
        df_results$variable,
        FUN = function(x) {
            new_str <- x
            for (i in independent_vars) {
                new_str <- sub(i, "", new_str)
                if (new_str == "") { # fix continuous var case
                    new_str <- i
                }
            }
            return(new_str)
        }
    ))
    # Strip the effect variable (e.g., sexM -> sex)
    df_results$variable <- unlist(apply(
        df_results[, c("variable", "variable_effect")],
        1,
        FUN = function(row) {
            if (row["variable_effect"] == row["variable"]) {
                return(row["variable"])
            } else {
                return(gsub(row["variable_effect"], "", row["variable"]))
            }
        }
    ))
    # Add info on the target variable
    if (!is.null(variable_target)) {
        df_results$variable_target <- variable_target
        df_results$model_reduced_anova_pvalue <- model_full$model_reduced_anova[[
            'Pr(>Chi)'
        ]][2]
        df_results$model_reduced_lrt_pvalue <- model_full[[
            "model_reduced_lrt_pvalue"
        ]]
    }
    # Add the formula of the model we are returning
    df_results$formula <- paste0("Y~", paste(independent_vars, collapse = "+"))

    # One could also add the AIC and BIC of the models
    # df_results$aic <- AIC(model_full)
    # df_results$bic <- BIC(model_full)

    # Add cell type list to the model
    model_full$cell_types <- cell_types

    return(list(
        "model" = model_full,
        "df_results" = df_results,
        "data_dr" = data
    ))
}


#' Plot results from \code{dirichlet_regression}.
#'
#' \code{plot_dirichlet_regression} generates plots to evaluate a
#' \code{dirichlet_regression} model.
#'
#' @param data_dr Data.frame.
#'     Model data from \code{dirichlet_regression}. data_dr$Y should be
#'     DirichletReg::DR_data. data_dr$experiment_id shoud correspond to
#'     samples. For numeric variables included in \code{model_dirichreg},
#'     each variable should have a corresponding '_unscaled' variable.
#' @param model_dirichreg DirichReg.model.
#'     DirichReg.model
#'
#' @return List.
#'     A list of plots
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
plot_dirichlet_regression <- function(
        data_dr,
        model_dirichreg
    ) {

    # Get the variables in model_dirichreg
    independent_vars <- attr(terms(model_dirichreg$formula), "term.labels")
    cell_types <- model_dirichreg$cell_types
    formula_str <- paste0("Y ~ ", paste0(independent_vars, collapse = " + "))

    # Intialize our return info
    plts_return <- list()
    plt_i <- 1

    # NOTE: details on model prediction:
    # * fitted() only returns values where we have observed data.
    # * predict() can be used for new data. Be sure to use type="response"
    #             flag.
    # See more at link below:
    # https://stackoverflow.com/questions/12201439/is-there-a-difference-between-the-r-functions-fitted-and-predict

    # Plot the predicted fit for each value ################################
    for (i in independent_vars) {
        # Make a dataframe of predicted data
        # df_plt <- data.frame(predict(
        #     model_dirichreg,
        #     data_dr,
        #     type = "response"
        # ))
        df_plt <- data.frame(fitted(model_dirichreg))
        colnames(df_plt) <- cell_types
        df_plt$experiment_id <- data_dr$experiment_id
        df_plt$x <- data_dr[[i]]
        if (is.numeric(data_dr[[i]])) {
            df_plt$x <- data_dr[[paste0(i, "_unscaled")]]
        }
        df_plt <- reshape2::melt(
            df_plt,
            id.vars = c("experiment_id", "x"),
            variable.name = "cell_type",
            value.name = "proportion_predicted"
        )

        # Make a dataframe of observed dadta
        sink("/dev/null")
        df_tmp <- data.frame(print(as.matrix(data_dr$Y)))
        sink()
        colnames(df_tmp) <- colnames(data_dr$Y)
        df_tmp$experiment_id <- rownames(data_dr)
        df_tmp <- reshape2::melt(
            df_tmp,
            id.vars = c("experiment_id"),
            variable.name = "cell_type",
            value.name = "proportion_observed"
        )

        # Combine the data frames
        df_plt$key <- as.character(paste(
            df_plt$experiment_id, df_plt$cell_type, sep = "-"
        ))
        df_tmp$key <- as.character(paste(
            df_tmp$experiment_id, df_tmp$cell_type, sep = "-"
        ))
        df_tmp$experiment_id <- NULL
        df_tmp$cell_type <- NULL
        df_plt <- merge(df_plt, df_tmp, by = c("key"), all.x = T, all.y = T)

        # Make the plot
        if (is.numeric(data_dr[[i]])) {
            plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                x = x,
                y = proportion_observed,
                color = cell_type
            ))
            plt <- plt + ggplot2::theme_bw(base_size = 12)
            plt <- plt + ggplot2::geom_point(alpha = 0.75)
            # plt <- plt + ggplot2::geom_line(data = df_plt, ggplot2::aes(
            #     x = x,
            #     y = proportion_predicted,
            #     color = cell_type
            # ))
            plt <- plt + ggplot2::geom_smooth(data = df_plt, ggplot2::aes(
                x = x,
                y = proportion_predicted,
                color = cell_type
            ))
            plt <- plt + ggplot2::labs(
                x = i,
                y = "Proportion",
                title = paste0(formula_str, ":\nModel fit to data")
            )
            plt <- plt + ggplot2::theme(legend.position = "none")
        } else {
            df_plt <- reshape2::melt(
                df_plt,
                id.vars = c("key", "experiment_id", "x", "cell_type"),
                variable.name = "proportion_type",
                value.name = "proportion"
            )
            plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                x = x,
                y = proportion,
                color = proportion_type
            ))
            plt <- plt + ggplot2::theme_bw(base_size = 12)
            plt <- plt + ggplot2::geom_boxplot(alpha = 0.75)
            plt <- plt + ggplot2::labs(
                x = i,
                y = "Proportion",
                color = "",
                fill = "",
                title = paste0(formula_str, ":\nModel fit to data")
            )
            # plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
            #     x = proportion_observed,
            #     y = proportion_predicted,
            #     color = x
            # ))
            # plt <- plt + ggplot2::geom_point(alpha = 0.75)
            # plt <- plt + ggplot2::geom_smooth(method = "lm")
        }
        plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
        plt <- plt + ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
        )
        plts_return[[plt_i]] <- plt
        plt_i <- plt_i + 1
    }

    # QQ plot of the residuals #############################################
    # There are three types of resduals:
    # 1. raw
    # 2. Pearson residuals (standardized)
    # 3. sum of the squared Pearson residuals (composite)
    # Because of the skewness and heteroscedasticity, raw residuals are
    # only of limited importance.
    sink("/dev/null")
    df_tmp <- data.frame(
        print(residuals(model_dirichreg, type = "standardized"))
    )
    sink()
    df_tmp$experiment_id <- rownames(df_tmp)
    df_plt <- reshape2::melt(
        df_tmp,
        id.vars = c("experiment_id"),
        variable.name = "cell_type",
        value.name = "residuals"
    )
    plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
        sample = residuals
    ))
    plt <- plt + ggplot2::stat_qq(
        alpha = 0.75,
        distribution = stats::qnorm
    )
    plt <- plt + ggplot2::stat_qq_line(distribution = stats::qnorm)
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::labs(
        x = "Theoretical quantiles",
        y = "Pearson residuals (standardized residuals)",
        title = paste0(formula_str, ":\nQQ plot of residuals")
    )
    plts_return[[plt_i]] <- plt
    plt_i <- plt_i + 1
    plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
    plts_return[[plt_i]] <- plt
    plt_i <- plt_i + 1

    # Plot residuals across all independent variables. #####################
    # Plot the predicted fit for each value
    for (i in independent_vars) {
        # Make a dataframe of residuals data
        sink("/dev/null")
        df_plt <- data.frame(
            print(residuals(model_dirichreg, type = "standardized"))
        )
        sink()
        df_plt$experiment_id <- rownames(df_plt)
        df_plt$x <- data_dr[[i]]
        if (is.numeric(data_dr[[i]])) {
            df_plt$x <- data_dr[[paste0(i, "_unscaled")]]
        }
        df_plt <- reshape2::melt(
            df_plt,
            id.vars = c("experiment_id", "x"),
            variable.name = "cell_type",
            value.name = "residuals"
        )
        if (is.numeric(data_dr[[i]])) {
            plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                x = x,
                y = residuals,
                color = cell_type
            ))
            plt <- plt + ggplot2::geom_point(alpha = 0.75)
            plt <- plt + ggplot2::geom_smooth()
        } else {
            plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                x = x,
                y = residuals,
                fill = x
            ))
            plt <- plt + ggplot2::geom_boxplot(outlier.shape = NA)
            plt <- plt + ggplot2::geom_jitter(alpha = 0.75)
        }
        plt <- plt + ggplot2::theme_bw(base_size = 12)
        plt <- plt + ggplot2::geom_hline(yintercept = 0)
        plt <- plt + ggplot2::theme(legend.position = "none")
        plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
        plt <- plt + ggplot2::labs(
            x = i,
            y = "Pearson residuals (standardized residuals)",
            title = paste0(formula_str, ":\nVariable by residuals")
        )
        plt <- plt + ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
        )
        plts_return[[plt_i]] <- plt
        plt_i <- plt_i + 1
    }

    # Plot the residuals against the fitted values. ########################
    # Make a dataframe of residuals
    sink("/dev/null")
    df_plt <- data.frame(
        print(residuals(model_dirichreg, type = "standardized"))
    )
    sink()
    df_plt$experiment_id <- rownames(df_plt)
    df_plt <- reshape2::melt(
        df_plt,
        id.vars = c("experiment_id"),
        variable.name = "cell_type",
        value.name = "residuals"
    )

    # Make a dataframe of predicted data
    df_tmp <- data.frame(fitted(model_dirichreg))
    colnames(df_tmp) <- cell_types
    df_tmp$experiment_id <- data_dr$experiment_id
    df_tmp <- reshape2::melt(
        df_tmp,
        id.vars = c("experiment_id"),
        variable.name = "cell_type",
        value.name = "proportion_predicted"
    )

    # Combine the data frames
    df_plt$key <- paste(
        df_plt$experiment_id, df_plt$cell_type, sep = "-"
    )
    df_tmp$key <- paste(
        df_tmp$experiment_id, df_tmp$cell_type, sep = "-"
    )
    df_tmp$experiment_id <- NULL
    df_tmp$cell_type <- NULL
    df_plt <- merge(df_plt, df_tmp, by = c("key"))

    # Make the plot
    plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
        x = proportion_predicted,
        y = residuals,
        color = cell_type
    ))
    plt <- plt + ggplot2::geom_point(alpha = 0.75)
    plt <- plt + ggplot2::geom_smooth()
    plt <- plt + ggplot2::geom_hline(yintercept = 0)
    plt <- plt + ggplot2::theme_bw(base_size = 12)
    plt <- plt + ggplot2::theme(legend.position = "none")
    plt <- plt + ggplot2::facet_wrap(~ cell_type, scales = "free")
    plt <- plt + ggplot2::labs(
        x = "Fitted values",
        y = "Pearson residuals (standardized residuals)",
        title = paste0(formula_str, ":\nFitted values by residuals")
    )
    plt <- plt + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0)
    )
    plts_return[[plt_i]] <- plt
    plt_i <- plt_i + 1

    return(plts_return)

} # end  make plots


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @import magrittr
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

        optparse::make_option(c("--formula_rhs"),
            type = "character",
            help = paste0(
                "Right hand side of formula for dirichlet model, starting with",
                " ~. Y is automatically set to propotions.",
                " Example 1: ~ disease_status.",
                " Example 2: ~ disease_status + age + sex.",
                " Note: Independent variable(s) can be continuous or ",
                " categorical.",
                " If categorical and 3 levels, then 1st level will be used as",
                " reference.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--variable_target"),
            type = "character",
            default = "",
            help = paste0(
                "Target variable in formula_rhs.",
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

        optparse::make_option(c("--forward_selection_variables"),
            type = "character",
            default = "",
            help = paste0(
                "Comma seperated list of variables to consider for",
                " stepwise forward model selection.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--optimization_var_threshold"),
            type = "numeric",
            default = 0.05,
            help = paste0(
                "Threshold for including a term in",
                " stepwise forward model selection.",
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
    formula_rhs <- param[["formula_rhs"]]
    variable_target <- param[["variable_target"]]
    continuous_variables_str <- param[["continuous_variables"]]
    discrete_variables_str <- param[["discrete_variables"]]
    levels_str <- param[["discrete_levels"]]
    forward_selection_variables_str <- param[["forward_selection_variables"]]
    optimization_var_threshold <- param[["optimization_var_threshold"]]
    base <- param[["out_file"]]

    # For dev
    # f_counts <- "diff_comp-counts.tsv.gz"
    # f_meta <- "diff_comp-metadata.tsv.gz"
    # formula_rhs <- "~ disease_status"
    # variable_target <- "disease_status"
    # forward_selection_variables_str <- "age,sex,I(age^2),smoker_at_time_of_biopsy"
    # continuous_variables_str <- "age,I(age^2)"
    # discrete_variables_str <- "disease_status,sex,smoker_at_time_of_biopsy"
    # levels_str <- "disease_status::healthy,cd"

    # Parse any forward selection variables
    forward_selection_variables <- c()
    if (forward_selection_variables_str != "") {
        forward_selection_variables <- strsplit(
            x = forward_selection_variables_str,
            split = ",",
            fixed = TRUE
        )[[1]]
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

    # Scale continuous variables included in the model or that have the
    # potential to be included in the model.
    #
    # First get a list of all terms in the model or that may be added to
    # the model
    formula_obj <- as.formula(
        paste("Y", gsub("~", "", formula_rhs), sep = "~")
    )
    possible_terms <- c(
        attr(terms(formula_obj), "term.labels"),
        forward_selection_variables
    )
    # First evaluate any terms that need to be evaluated on original scale
    for (i in possible_terms) {
        # If the term is one that we need to evaluate and generate, then do so
        if (grepl("I(", i, fixed = TRUE)) {
            cat(paste0("Evaluating:\t", i, "\n"))
            df_meta[[i]] <- eval(
                parse(text = i),
                df_meta
            )
        }
    }
    # Now all terms are set, scale any terms that need to be scaled
    for (i in possible_terms) {
        # If the term is numeric, then scale it
        if (is.numeric(df_meta[[i]])) {
            cat(paste0("Scaling:\t", i, "\n"))
            # Save the unscaled variable as _unscaled
            df_meta[[paste0(i, "_unscaled")]] <- df_meta[[i]]
            df_meta[[i]] <- scale(
                df_meta[[paste0(i, "_unscaled")]],
                center = TRUE,
                scale = TRUE
            )
        }
    }


    # Fit the starting model.
    results <- dirichlet_regression(
        df_counts = df_counts,
        df_independent_vars = df_meta,
        formula_rhs = formula_rhs,
        variable_target = variable_target,
        fit_ols = FALSE,
        forward_selection_variables = c()
    )

    # Save the R object of the results
    saveRDS(
        results$model,
        file = paste0(base, "-dirichletreg_model.Rds.gz"),
        compress = TRUE
    )

    # Set a useful dataframe of the results
    df_results <- results$df_results
    df_results$method <- "dirichlet_regression"

    # Make plots
    plts <- plot_dirichlet_regression(
        results$data_dr,
        results$model
    )
    if (length(plts) > 0) {
        pdf(
            file = paste(base, "dirichletreg_model.pdf", sep = "-"),
            height = 15,
            width = 15
        )
            # # Plot the DR data
            # # NOTE: requires rgl package
            # print(plot(results$data_dr$Y))
            for (i in seq(length(plts))) {
                print(plts[[i]])
            }
        dev.off()
    }


    # Optionally run the forward selection model and:
    # 1. Save the model
    # 2. Plot the results
    # 3. Add in the results to the previous model
    if (length(forward_selection_variables) > 0) {
        results_forward <- dirichlet_regression(
            df_counts = df_counts,
            df_independent_vars = df_meta,
            formula_rhs = formula_rhs,
            variable_target = variable_target,
            fit_ols = FALSE,
            forward_selection_variables = forward_selection_variables,
            optimization_var_threshold = optimization_var_threshold
        )

        # Save the R object of the results
        saveRDS(
            results_forward$model,
            file = paste0(base, "-dirichletreg_forwardmodel.Rds.gz"),
            compress = TRUE
        )

        # Make plots
        plts <- plot_dirichlet_regression(
            results_forward$data_dr,
            results_forward$model
        )
        if (length(plts) > 0) {
            pdf(
                file = paste(base, "dirichletreg_forwardmodel.pdf", sep = "-"),
                height = 15,
                width = 15
            )
                # # Plot the DR data
                # # NOTE: requires rgl package
                # print(plot(results$data_dr$Y))
                for (i in seq(length(plts))) {
                    print(plts[[i]])
                }
            dev.off()
        }

        # Add in results to previous model
        results_forward$df_results$method <- "dirichlet_regression-step_forward"
        df_results <- rbind(
            df_results,
            results_forward$df_results
        )
    }


    # Save the model results
    df_results <- df_results[
        with(df_results, order(pvalue, abs(estimate))),
    ]
    # Perform FDR for each term independently
    df_results <- df_results %>%
        dplyr::group_by(variable, variable_effect, formula, method) %>%
        dplyr::mutate(
            n = n(),  # Check to make sure == to number cell types
            qvalue_bh = p.adjust(pvalue, method = "BH")
        ) %>%
        as.data.frame
    if (length(unique(subset(df_results, variable != "intercept")$n)) != 1) {
        print(subset(df_results, n > 1 & variable != "intercept"))
        stop("ERROR: length(unique(df_results$n)) != 1")
    } else {
        df_results$n <- NULL
    }
    gzfh <- gzfile(
        paste0(base, "-dirichletreg_results.tsv.gz"),
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
