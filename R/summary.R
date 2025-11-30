#' Compute and print a summary while returning a list you can access with $
summarise_obs_data <- function(observed_data, print = TRUE) {
  need_cols <- c("status_ill", "status_dead", "V_healthy", "T_obs")
  missing_cols <- setdiff(need_cols, names(observed_data))
  if (length(missing_cols)) stop("Missing columns in datasets$obs: ", paste(missing_cols, collapse = ", "))

  # --- vectors
  n <- nrow(observed_data)
  observed_illness <- as.logical(observed_data$status_ill)
  observed_death   <- as.logical(observed_data$status_dead)
  unknown_transition <- !observed_illness & (observed_data$V_healthy < observed_data$T_obs)

  # --- counts
  n_any_death                 <- sum(observed_death, na.rm = TRUE)
  n_death_via_illness         <- sum(observed_death & observed_illness, na.rm = TRUE)
  n_death_exact               <- sum(observed_death & !unknown_transition & !observed_illness, na.rm = TRUE)
  n_death_missing_transition  <- sum(observed_death &  unknown_transition & !observed_illness, na.rm = TRUE)

  n_censoring_exact               <- sum(!observed_death & !unknown_transition, na.rm = TRUE)
  n_censoring_missing_transition  <- sum(!observed_death &  unknown_transition & !observed_illness, na.rm = TRUE)

  # --- tables
  table_illness <- table(Illness = observed_illness, Death = observed_death, useNA = "ifany")
  table_death   <- table(
    Death = observed_death[!observed_illness],
    LastStageKnown = !unknown_transition[!observed_illness],
    useNA = "ifany"
  )

  n_missing = n_censoring_missing_transition + n_death_missing_transition
  n_ill = sum(observed_data$status_ill)
  n_no_ill_obs = n - sum(observed_data$status_ill)

  out <- list(
    n = n,
    counts = list(
      n_any_death = n_any_death,
      n_death_via_illness = n_death_via_illness,
      n_death_exact = n_death_exact,
      n_death_missing_transition = n_death_missing_transition,
      n_censoring_exact = n_censoring_exact,
      n_censoring_missing_transition = n_censoring_missing_transition,
      n_missing = n_missing,
      n_ill = n_ill,
      n_no_ill_obs = n_no_ill_obs,
      p_missing = n_missing/n_no_ill_obs,
      p_missing_any = n_missing/n
    ),
    tables = list(
      table_illness = table_illness,
      table_death = table_death
    ),
    flags = list(
      observed_illness = observed_illness,
      observed_death = observed_death,
      unknown_transition = unknown_transition
    )
  )
  class(out) <- c("summary_data", class(out))

  if (isTRUE(print)) print(out)
  out
}

true_path_for_unknown <- function(exact, obs) {
  need_cols <- c("status_ill", "status_dead", "V_healthy", "T_obs")
  missing_cols <- setdiff(need_cols, names(obs))
  if (length(missing_cols)) stop("Missing columns in datasets$obs: ", paste(missing_cols, collapse = ", "))

  # --- vectors
  n <- nrow(obs)
  observed_illness <- as.logical(obs$status_ill)
  observed_death   <- as.logical(obs$status_dead)
  unknown_transition <- !observed_illness & (obs$V_healthy < obs$T_obs)

  illness_happened_before_T_obs <- exact$T_ill <= obs$T_obs


  exact$path_w_censoring <- ifelse(!observed_death & !observed_illness,
                                   ifelse(illness_happened_before_T_obs, "1->2->?", "1->?"),
                                   NA)

  table_path_unknown_transition <- NULL
  if (!is.null(exact$path)) {
    # assumes exact has 1:1 rows with obs; adapt if your key is different
    idx <- !observed_illness & unknown_transition
    idx <- idx & seq_len(n) <= length(exact$path)  # guard if lengths differ

    exact$path_w_censoring[observed_death] <- as.character(exact$path[observed_death])
    exact$path_w_censoring[observed_illness & !observed_death] <- "1->2->?"

    idx = ifelse(idx, "unknown", "known")
    table_path_unknown_transition <- table(Path = exact$path_w_censoring, "Transistion" = idx, useNA = "ifany")
  }
  table_path_unknown_transition
}




#' Pretty printing for summary_data objects
print.summary_data <- function(x, ...) {
  cat("Summary\n")
  cat("-------\n")
  cat("N observations:", x$n, "\n\n")

  cat("Counts\n")
  cat("  Any death:                          ", x$counts$n_any_death, "\n")
  cat("  Death via illness:                  ", x$counts$n_death_via_illness, "\n")
  cat("  Death (exact, no prior illness):    ", x$counts$n_death_exact, "\n")
  cat("  Death (missing transition):         ", x$counts$n_death_missing_transition, "\n")
  cat("  Censoring (exact):                  ", x$counts$n_censoring_exact, "\n")
  cat("  Censoring (missing transition):     ", x$counts$n_censoring_missing_transition, "\n\n")

  cat("Tables\n")
  cat("$tables$table_illness\n"); print(addmargins(x$tables$table_illness))
  cat("\n$tables$table_death\n");   print(addmargins(x$tables$table_death))
  if (!is.null(x$tables$table_path_unknown_transition)) {
    cat("\n$tables$table_path_unknown_transition\n")
    print(addmargins(x$tables$table_path_unknown_transition))
  }
  invisible(x)
}
