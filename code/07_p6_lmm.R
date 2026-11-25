# Frozen SAP S10 Aim 1 and G5 6.1a: visit-only random effects with Satterthwaite df.
arguments <- commandArgs(trailingOnly = TRUE)
tryCatch({
    stopifnot(length(arguments) == 8L)
    local_library <- arguments[[1L]]
    .libPaths(c(local_library, .Library), include.site = FALSE)
    stopifnot(requireNamespace("lme4", lib.loc = local_library, quietly = TRUE),
              requireNamespace("lmerTest", lib.loc = local_library, quietly = TRUE))
    source(arguments[[2L]], local = TRUE)
    control <- list(output = arguments[[3L]], hours = as.numeric(arguments[4:5]),
                    df = as.integer(arguments[[6L]]), confidence = as.numeric(arguments[[7L]]))
    set.seed(as.integer(arguments[[8L]]))
    selected <- pre_dose_frames(read_input(), control)$before_dose
    complete <- selected[stats::complete.cases(selected[, c("prior_htn", "age_years", "sex")]), ]
    basis <- new_basis(complete, control)
    complete <- add_basis(complete, basis)
    spline_terms <- paste0("spline", seq_len(control$df), collapse = " + ")
    formula <- stats::as.formula(paste("sbp ~ (", spline_terms,
                            ") * prior_htn + age_years + sex + (1 + time_hours | visit_id)"))
    result <- fit_attempt(lmerTest::lmer(formula, data = complete, REML = TRUE, na.action = stats::na.fail,
                        control = lme4::lmerControl(optimizer = "bobyqa")))
    stopifnot(!is.null(result$model))
    model <- result$model
    coefficients <- summary(model, ddf = "Satterthwaite")$coefficients
    critical <- stats::qt((1 + control$confidence) / 2, coefficients[, "df"])
    estimates <- data.frame(term = rownames(coefficients), estimate = coefficients[, "Estimate"],
                            std_error = coefficients[, "Std. Error"], df = coefficients[, "df"],
                            ci_low = coefficients[, "Estimate"] - critical * coefficients[, "Std. Error"],
                            ci_high = coefficients[, "Estimate"] + critical * coefficients[, "Std. Error"],
                            p_value = coefficients[, "Pr(>|t|)"])
    term_names <- names(lme4::fixef(model))
    indices <- which(grepl(":", term_names, fixed = TRUE))
    contrast <- diag(length(term_names))[indices, , drop = FALSE]
    joint <- lmerTest::contestMD(model, contrast, ddf = "Satterthwaite")
    info <- model@optinfo$conv
    covariance <- as.matrix(lme4::VarCorr(model)$visit_id)
    diagnostics <- data.frame(
        as.list(subset_counts(complete, "modeled")),
        spline_knot_1_hours = basis$knots[[1L]], spline_knot_2_hours = basis$knots[[2L]],
        converged = isTRUE(info$opt == 0L) && !length(info$lme4$messages),
        singular = lme4::isSingular(model), optimizer = "bobyqa", warning_count = result$warnings,
        reml_loglik = as.numeric(stats::logLik(model)), aic = stats::AIC(model), residual_sd = stats::sigma(model),
        visit_intercept_sd = sqrt(covariance[1L, 1L]), visit_slope_sd_mmhg_per_hour = sqrt(covariance[2L, 2L]),
        visit_intercept_slope_corr = stats::cov2cor(covariance)[1L, 2L],
        joint_interaction_f = joint[, "F value"], joint_interaction_numerator_df = joint[, "NumDF"],
        joint_interaction_denominator_df = joint[, "DenDF"], joint_interaction_p_value = joint[, "Pr(>F)"]
    )
    write_aggregate(estimates, control$output, "coefficients")
    write_aggregate(diagnostics, control$output, "diagnostics")
    cat("P6_LMM_COMPLETED\n")
}, error = function(e) {
    cat("P6_LMM_FAILED\n")
    quit(save = "no", status = 1L, runLast = FALSE)
})
