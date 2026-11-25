# SAP v0.3 S10 and S5.6. Read the outcome-free modeling input on stdin.
options(warn = 1)

write_aggregate <- function(data, path, name) {
    utils::write.csv(data, file.path(path, paste0(name, ".csv")), row.names = FALSE, na = "")
}

parse_flag <- function(values) {
    values <- tolower(trimws(as.character(values)))
    ifelse(values %in% c("true", "t", "1"), 1L,
           ifelse(values %in% c("false", "f", "0"), 0L, NA_integer_))
}

read_input <- function() {
    connection <- file("stdin", open = "rt", encoding = "UTF-8")
    on.exit(close(connection))
    data <- utils::read.csv(connection, colClasses = "character", na.strings = c("", "NA"))
    columns <- c("visit_id", "patient_id", "minutes_from_arrival", "sbp", "prior_htn", "age_years",
                 "sex", "has_medication_list", "ed_acute_aht", "ed_first_acute_aht_minutes")
    stopifnot(identical(names(data), columns), nrow(data) > 0L)
    for (column in c("minutes_from_arrival", "sbp", "age_years", "ed_first_acute_aht_minutes")) {
        data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
    }
    for (column in c("prior_htn", "has_medication_list", "ed_acute_aht")) {
        data[[column]] <- parse_flag(data[[column]])
    }
    sex <- toupper(trimws(data$sex))
    data$sex <- factor(ifelse(sex %in% c("F", "FEMALE"), "female",
                       ifelse(sex %in% c("M", "MALE"), "male", NA_character_)),
                       levels = c("female", "male"))
    stopifnot(!anyNA(data$visit_id), !anyNA(data$patient_id))
    data$time_hours <- data$minutes_from_arrival / 60
    data[order(data$patient_id, data$visit_id, data$time_hours, data$sbp), , drop = FALSE]
}

new_basis <- function(data, control) {
    knots <- stats::quantile(data$time_hours, probs = seq_len(control$df - 1L) / control$df, names = FALSE)
    stopifnot(length(unique(knots)) == length(knots),
              all(knots > control$hours[1L]), all(knots < control$hours[2L]))
    list(knots = knots, boundary = control$hours)
}

add_basis <- function(data, basis) {
    spline <- splines::ns(data$time_hours, knots = basis$knots, Boundary.knots = basis$boundary)
    colnames(spline) <- paste0("spline", seq_len(ncol(spline)))
    cbind(data, spline)
}

fit_attempt <- function(expression) {
    warnings <- 0L
    model <- tryCatch(withCallingHandlers(expression, warning = function(w) {
        warnings <<- warnings + 1L
        invokeRestart("muffleWarning")
    }), error = function(e) NULL)
    list(model = model, warnings = warnings)
}

counts <- function(data) {
    c(reading_n = nrow(data), visit_n = length(unique(data$visit_id)),
      patient_n = length(unique(data$patient_id)))
}

subset_counts <- function(data, prefix) {
    values <- counts(data)
    setNames(values, paste0(prefix, "_", names(values)))
}

fit_lmm <- function(data, fixed) {
    for (hierarchy in c("patient_intercept_visit_intercept_slope", "visit_intercept_slope")) {
        random <- if (hierarchy == "patient_intercept_visit_intercept_slope") {
            list(patient_id = ~ 1, visit_id = ~ time_hours)
        } else {
            list(visit_id = ~ time_hours)
        }
        for (optimizer in c("nlminb", "optim")) {
            fit <- fit_attempt(nlme::lme(
                fixed = fixed, random = random, data = data, method = "REML",
                na.action = stats::na.fail,
                control = nlme::lmeControl(opt = optimizer, maxIter = 150, msMaxIter = 200,
                                           msMaxEval = 400, allow.n.lt.q = TRUE)
            ))
            if (!is.null(fit$model)) {
                return(c(fit, list(hierarchy = hierarchy, optimizer = optimizer)))
            }
        }
    }
    list(model = NULL, warnings = NA_integer_, hierarchy = NA_character_, optimizer = NA_character_)
}

wald_curve <- function(fit, data, fixed, curve_times, basis, control, specification) {
    covariates <- data[!duplicated(data$visit_id), c("age_years", "sex")]
    design_terms <- stats::delete.response(stats::terms(fixed))
    calculate <- function(group, hour) {
        new_data <- add_basis(data.frame(covariates, prior_htn = group, time_hours = hour), basis)
        colMeans(stats::model.matrix(design_terms, data = new_data))
    }
    vif <- stats::vcov(fit)
    beta <- nlme::fixef(fit)
    critical <- stats::qnorm((1 + control$confidence) / 2)
    effect <- function(design) {
        design <- design[names(beta)]
        estimate <- sum(design * beta)
        std_error <- sqrt(as.numeric(t(design) %*% vif %*% design))
        c(estimate = estimate, ci_low = estimate - critical * std_error, ci_high = estimate + critical * std_error)
    }
    output <- lapply(0:1, function(group) {
        baseline <- calculate(group, curve_times[1L])
        do.call(rbind, lapply(curve_times, function(hour) {
            design <- calculate(group, hour)
            mean_values <- effect(design)
            change <- effect(design - baseline)
            data.frame(specification = specification, prior_htn = group, time_hours = hour,
                       estimate_mmhg = mean_values[[1L]], ci_low = mean_values[[2L]], ci_high = mean_values[[3L]],
                       change_from_arrival_mmhg = change[[1L]], change_ci_low = change[[2L]], change_ci_high = change[[3L]])
        }))
    })
    do.call(rbind, output)
}

run_lmm <- function(window, before_dose, control) {
    coefficients <- list()
    diagnostics <- list()
    curves <- list()
    fixed <- stats::as.formula(paste("sbp ~ (", paste0("spline", seq_len(control$df), collapse = " + "),
                                     ") * prior_htn + age_years + sex"))
    primary_complete <- before_dose[stats::complete.cases(before_dose[, c("prior_htn", "age_years", "sex")]), ]
    basis <- new_basis(primary_complete, control)
    curve_times <- seq(control$hours[1L], control$hours[2L], length.out = control$grid_points)
    specs <- c("pre_acute_dose", "medication_list_only", "exclude_acute_treated")
    for (specification in specs) {
        eligible <- switch(specification,
            pre_acute_dose = window,
            medication_list_only = window[which(window$has_medication_list == 1L), ],
            exclude_acute_treated = window[which(window$ed_acute_aht == 0L), ]
        )
        selected <- before_dose[before_dose$visit_id %in% eligible$visit_id, , drop = FALSE]
        complete <- primary_complete[primary_complete$visit_id %in% selected$visit_id, , drop = FALSE]
        complete <- add_basis(complete, basis)
        fit_result <- fit_lmm(complete, fixed)
        fit <- fit_result$model
        diagnostic <- data.frame(
            specification = specification,
            as.list(subset_counts(eligible, "window")), as.list(subset_counts(selected, "pre_dose")),
            as.list(subset_counts(complete, "modeled")),
            unknown_acute_time_reading_n = sum(eligible$ed_acute_aht == 1L & !is.finite(eligible$ed_first_acute_aht_minutes), na.rm = TRUE),
            unknown_acute_flag_reading_n = sum(is.na(eligible$ed_acute_aht)),
            missing_covariate_reading_n = nrow(selected) - nrow(complete),
            spline_knot_1_hours = basis$knots[[1L]], spline_knot_2_hours = basis$knots[[2L]],
            status = if (is.null(fit)) "failed" else "converged", random_effects = fit_result$hierarchy,
            optimizer = fit_result$optimizer, warning_count = fit_result$warnings,
            reml_loglik = NA_real_, aic = NA_real_, residual_sd = NA_real_,
            patient_intercept_sd = NA_real_, visit_intercept_sd = NA_real_,
            visit_slope_sd_mmhg_per_hour = NA_real_, visit_intercept_slope_corr = NA_real_,
            standardized_residual_q025 = NA_real_, standardized_residual_q50 = NA_real_,
            standardized_residual_q975 = NA_real_, joint_interaction_f = NA_real_,
            joint_interaction_numerator_df = NA_real_, joint_interaction_denominator_df = NA_real_,
            joint_interaction_p_value = NA_real_
        )
        if (!is.null(fit)) {
            tab <- summary(fit)$tTable
            critical <- stats::qt((1 + control$confidence) / 2, df = tab[, "DF"])
            coefficients[[specification]] <- data.frame(
                specification = specification, term = rownames(tab), estimate = tab[, "Value"],
                std_error = tab[, "Std.Error"], df = tab[, "DF"],
                ci_low = tab[, "Value"] - critical * tab[, "Std.Error"],
                ci_high = tab[, "Value"] + critical * tab[, "Std.Error"], p_value = tab[, "p-value"]
            )
            diagnostic$reml_loglik <- as.numeric(stats::logLik(fit))
            diagnostic$aic <- stats::AIC(fit)
            diagnostic$residual_sd <- fit$sigma
            if (!is.null(fit$modelStruct$reStruct$patient_id)) {
                patient_cov <- as.matrix(fit$modelStruct$reStruct$patient_id) * fit$sigma ^ 2
                diagnostic$patient_intercept_sd <- sqrt(patient_cov[1L, 1L])
            }
            visit_cov <- as.matrix(fit$modelStruct$reStruct$visit_id) * fit$sigma ^ 2
            diagnostic$visit_intercept_sd <- sqrt(visit_cov[1L, 1L])
            diagnostic$visit_slope_sd_mmhg_per_hour <- sqrt(visit_cov[2L, 2L])
            diagnostic$visit_intercept_slope_corr <- stats::cov2cor(visit_cov)[1L, 2L]
            residual_quantiles <- stats::quantile(stats::residuals(fit, type = "normalized"),
                                                  probs = c((1 - control$confidence) / 2, 0.5, (1 + control$confidence) / 2))
            diagnostic$standardized_residual_q025 <- residual_quantiles[[1L]]
            diagnostic$standardized_residual_q50 <- residual_quantiles[[2L]]
            diagnostic$standardized_residual_q975 <- residual_quantiles[[3L]]
            term_names <- names(nlme::fixef(fit))
            interaction_terms <- which(grepl(":", term_names, fixed = TRUE))
            contrast <- diag(length(term_names))[interaction_terms, , drop = FALSE]
            colnames(contrast) <- term_names
            joint <- stats::anova(fit, L = contrast)
            diagnostic$joint_interaction_f <- joint[, "F-value"]
            diagnostic$joint_interaction_numerator_df <- joint[, "numDF"]
            diagnostic$joint_interaction_denominator_df <- joint[, "denDF"]
            diagnostic$joint_interaction_p_value <- joint[, "p-value"]
            curves[[specification]] <- wald_curve(fit, complete, fixed, curve_times, basis, control, specification)
        }
        diagnostics[[specification]] <- diagnostic
    }
    write_aggregate(do.call(rbind, diagnostics), control$output, "p2_mcmed_lmm_diagnostics")
    if (length(coefficients)) {
        write_aggregate(do.call(rbind, coefficients), control$output, "p2_mcmed_lmm_coefficients")
        write_aggregate(do.call(rbind, curves), control$output, "p2_mcmed_lmm_adjusted_curves")
    }
}

lcmm_class_summary <- function(model) {
    posterior <- model$pprob
    total <- nrow(posterior)
    do.call(rbind, lapply(seq_len(model$ng), function(group) {
        assigned <- !is.na(posterior$class) & posterior$class == group
        assigned_count <- sum(assigned)
        data.frame(k = model$ng, class = group, visits = assigned_count, proportion = assigned_count / total,
                   mean_assigned_posterior = if (assigned_count) mean(posterior[assigned, paste0("prob", group)]) else NA_real_)
    }))
}

run_lcmm <- function(before_dose, control) {
    stopifnot(requireNamespace("lcmm", quietly = TRUE))
    suppressWarnings(suppressPackageStartupMessages(library(lcmm)))
    basis <- new_basis(before_dose, control)
    data <- add_basis(before_dose, basis)
    # hlme requires a numeric subject identifier. This mapping stays in memory.
    data$visit_number <- match(data$visit_id, unique(data$visit_id))
    spline_terms <- paste0("spline", seq_len(control$df), collapse = " + ")
    fixed <- stats::as.formula(paste("sbp ~", spline_terms))
    mixture <- stats::as.formula(paste("~", spline_terms))
    models <- list()
    summaries <- list()
    diagnostics <- list()
    initial <- fit_attempt(hlme(fixed, random = ~ time_hours, subject = "visit_number", data = data,
                                 maxiter = control$final_iterations))
    for (k in control$class_counts) {
        set.seed(control$seed + k)
        method <- "single_class"
        if (k == 1L) {
            result <- initial
        } else if (!is.null(initial$model) && initial$model$conv == 1L) {
            method <- "gridsearch"
            result <- fit_attempt(gridsearch(
                hlme(fixed, mixture = mixture, random = ~ time_hours, subject = "visit_number", ng = k,
                     data = data, maxiter = control$final_iterations),
                rep = control$starts, maxiter = control$grid_iterations, minit = initial$model
            ))
            if (is.null(result$model)) {
                method <- "single_start_fallback"
                result <- fit_attempt(hlme(fixed, mixture = mixture, random = ~ time_hours, subject = "visit_number", ng = k,
                                           data = data, B = initial$model, maxiter = control$final_iterations))
            }
        } else {
            method <- "single_class_fit_unavailable"
            result <- list(model = NULL, warnings = NA_integer_)
        }
        model <- result$model
        key <- as.character(k)
        valid_posteriors <- !is.null(model) && !is.null(model$pprob) && nrow(model$pprob) == length(unique(data$visit_id))
        groups <- if (valid_posteriors) lcmm_class_summary(model) else NULL
        min_proportion <- if (is.null(groups)) NA_real_ else min(groups$proportion)
        min_posterior <- if (is.null(groups)) NA_real_ else min(groups$mean_assigned_posterior)
        eligible <- !is.null(model) && model$conv == 1L && is.finite(model$BIC) &&
            is.finite(min_proportion) && min_proportion >= control$min_proportion &&
            is.finite(min_posterior) && min_posterior >= control$min_posterior
        if (!is.null(model)) {
            models[[key]] <- model
        }
        if (!is.null(groups)) {
            summaries[[key]] <- groups
        }
        diagnostics[[key]] <- data.frame(
            k = k, as.list(counts(data)),
            spline_knot_1_hours = basis$knots[[1L]], spline_knot_2_hours = basis$knots[[2L]],
            start_method = method, grid_starts = if (method == "gridsearch") control$starts else NA_integer_,
            grid_max_iterations = if (method == "gridsearch") control$grid_iterations else NA_integer_,
            final_max_iterations = control$final_iterations,
            converged = !is.null(model) && model$conv == 1L,
            convergence_code = if (is.null(model)) NA_integer_ else model$conv,
            iterations = if (is.null(model)) NA_integer_ else model$niter,
            warning_count = result$warnings,
            log_likelihood = if (is.null(model)) NA_real_ else model$loglik,
            aic = if (is.null(model)) NA_real_ else model$AIC,
            bic = if (is.null(model)) NA_real_ else model$BIC,
            min_class_proportion = min_proportion, min_mean_assigned_posterior = min_posterior,
            eligible = eligible
        )
        write_aggregate(do.call(rbind, diagnostics), control$output, "p2_mcmed_lcmm_models")
    }
    selection <- do.call(rbind, diagnostics)
    accepted <- selection[selection$eligible, , drop = FALSE]
    selected_k <- if (nrow(accepted)) accepted$k[which.min(accepted$bic)] else NA_integer_
    selection$selected <- !is.na(selected_k) & selection$k == selected_k
    write_aggregate(selection, control$output, "p2_mcmed_lcmm_models")
    if (length(summaries)) {
        class_data <- do.call(rbind, summaries)
        class_data$selected <- !is.na(selected_k) & class_data$k == selected_k
        write_aggregate(class_data, control$output, "p2_mcmed_lcmm_classes")
    }
    if (is.na(selected_k)) {
        return(invisible(NULL))
    }
    grid <- add_basis(data.frame(time_hours = seq(control$hours[1L], control$hours[2L], length.out = control$grid_points)), basis)
    chosen <- models[[as.character(selected_k)]]
    # predictY inspects the formula expressions stored in the call.
    chosen$call$fixed <- fixed
    if (selected_k > 1L) {
        chosen$call$mixture <- mixture
    }
    curves <- lcmm::predictY(chosen, grid, var.time = "time_hours", draws = TRUE)$pred
    trajectory <- do.call(rbind, lapply(seq_len(selected_k), function(group) {
        prefix <- if (selected_k == 1L) "Ypred" else paste0("Ypred_class", group)
        data.frame(
            k = selected_k, class = group, time_hours = grid$time_hours,
            predicted_sbp_mmhg = curves[, prefix], ci_low = curves[, paste0("lower.", prefix)],
            ci_high = curves[, paste0("upper.", prefix)]
        )
    }))
    write_aggregate(trajectory, control$output, "p2_mcmed_lcmm_trajectories")
}

pre_dose_frames <- function(data, control) {
    eligible <- is.finite(data$sbp) & is.finite(data$time_hours) &
        data$time_hours >= control$hours[1L] & data$time_hours <= control$hours[2L]
    window <- data[eligible, , drop = FALSE]
    # A reading tied to the first dose does not establish a pre-dose value.
    before <- window$ed_acute_aht == 0L | (window$ed_acute_aht == 1L &
        is.finite(window$ed_first_acute_aht_minutes) & window$minutes_from_arrival < window$ed_first_acute_aht_minutes)
    before_dose <- window[which(before), , drop = FALSE]
    list(window = window, before_dose = before_dose)
}

main <- function(control) {
    frames <- pre_dose_frames(read_input(), control)
    stage <<- "lmm"
    run_lmm(frames$window, frames$before_dose, control)
    stage <<- "lcmm"
    run_lcmm(frames$before_dose, control)
}

p2_cli <- function(arguments) {
    stage <<- "setup"
    tryCatch({
    stopifnot(length(arguments) == 14L)
    .libPaths(c(arguments[[1L]], .Library), include.site = FALSE)
    stopifnot(requireNamespace("nlme", quietly = TRUE))
    control <- list(
        output = arguments[[2L]], seed = as.integer(arguments[[3L]]), hours = as.numeric(arguments[4:5]),
        df = as.integer(arguments[[6L]]), grid_points = as.integer(arguments[[7L]]),
        confidence = as.numeric(arguments[[8L]]), class_counts = as.integer(strsplit(arguments[[9L]], ",", fixed = TRUE)[[1L]]),
        min_proportion = as.numeric(arguments[[10L]]), min_posterior = as.numeric(arguments[[11L]]),
        starts = as.integer(arguments[[12L]]), grid_iterations = as.integer(arguments[[13L]]),
        final_iterations = as.integer(arguments[[14L]])
    )
    stage <<- "input"
    main(control)
    cat("P2_MODELS_COMPLETED\n")
}, error = function(e) {
    cat(paste0("P2_MODELS_FAILED_", stage, "\n"))
    quit(save = "no", status = 1L, runLast = FALSE)
    })
}

if (sys.nframe() == 0L) {
    p2_cli(commandArgs(trailingOnly = TRUE))
}
