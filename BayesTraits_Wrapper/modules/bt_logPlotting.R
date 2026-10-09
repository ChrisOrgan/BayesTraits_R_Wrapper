# Unified BayesTraits Log Plotting

# Load libraries
library(ggplot2)
library(stringr)
library(bayestestR)
library(coda)
library(dplyr)
library(tibble)
library(rcartocolor)

options(scipen = 999)


# UTILITY FUNCTIONS ----------------------------------------------------
# Read a single BayesTraits log file into a data.frame
read_bt_log <- function(file_path) {
  tbl <- readLines(file_path)
  start_line <- grep("Iteration.*Lh|Lh.*Iteration", tbl)[1]
  df <- read.delim(file_path, header = TRUE, skip = start_line - 1, na.strings = c("N\\A", "NA", "-999", ""))
  df <- df[, colSums(!is.na(df)) > 0]
  df$Tree.No <- NULL
  return(df)
}

keywords <- c("Lh", "Beta", "Alpha", "Sigma", "q", "R.")

# Extract relevant columns by keyword
get_cols <- function(df) {
  res <- list()
  for (kw in keywords) {
    if (kw == "Sigma") {
      # Match Sigma or Var, case-insensitive
      matched <- grep("^(Sigma|Var)", names(df), value = TRUE, ignore.case = TRUE)
    } else {
      matched <- grep(paste0("^", kw), names(df), value = TRUE, ignore.case = TRUE)
    }
    if (length(matched) > 0) res[[kw]] <- matched
  }
  return(res)
}

# Extract signature for grouping multiple runs
extract_relevant <- function(file) {
  lines <- readLines(file, warn = FALSE)
  get_same_line <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) > 0) trimws(line[1]) else NA
  }
  tree_file  <- get_same_line("Tree File Name")
  data_file  <- get_same_line("Data File Name")
  sample_per <- get_same_line("Sample Period")
  iterations <- get_same_line("Iterations")
  burnin     <- get_same_line("Burn In")

  prior_start <- grep("Prior Information", lines)
  tree_start  <- grep("Tree Information", lines)

  prior_block <- NA
  if (length(prior_start) > 0 && length(tree_start) > 0) {
    prior_block <- paste(trimws(lines[(prior_start + 1):(tree_start - 1)]), collapse = "\n")
  }

  paste(tree_file, data_file, sample_per, iterations, burnin, prior_block, sep = "||")
}

make_results_table <- function(df, param_cols, diagnostics = "single") {

  results <- data.frame(
    Mean    = sapply(param_cols, function(nm) mean(df[[nm]])),
    Q2.5    = sapply(param_cols, function(nm) quantile(df[[nm]], 0.025)),
    Q50     = sapply(param_cols, function(nm) quantile(df[[nm]], 0.5)),
    Q97.5   = sapply(param_cols, function(nm) quantile(df[[nm]], 0.975)),
    HPD2.5  = sapply(param_cols, function(nm) unlist(ci(df[[nm]], method = "HDI"))[2]),
    HPD97.5 = sapply(param_cols, function(nm) unlist(ci(df[[nm]], method = "HDI"))[3]),
    Above0  = sapply(param_cols, function(nm) mean(df[[nm]] > 0))
  )

  rownames(results) <- param_cols

  if (diagnostics == "single") {
    results$ESS <- sapply(param_cols, function(nm) {
      if (length(unique(df[[nm]])) > 1) effectiveSize(df[[nm]]) else NA
    })

    results$Geweke <- sapply(param_cols, function(nm) {
      if (length(unique(df[[nm]])) > 1) geweke.diag(df[[nm]])$z else NA
    })

    results$Converged <- ifelse(
      is.na(results$ESS), NA,
      ifelse(results$ESS > 200 & abs(results$Geweke) < 2, "Yes", "Maybe")
    )
  }

  if (diagnostics == "multi") {

    save_runID <- df$num_runID
    combined_numeric <- df %>% select(where(is.numeric))
    vars <- sapply(combined_numeric, var)
    combined_var <- combined_numeric[, which(vars > 0.0001)]

    combined_var$num_runID <- save_runID
    if (ncol(combined_var) == 0) {
      print("Parameter and Lh Variances too low to do Gelman diagnostic")

    } else {

      df_list <- split(combined_var, combined_var$num_runID)
      df_list <- lapply(df_list, function(x) {
        rownames(x) <- NULL
        x %>%
          column_to_rownames("Iteration") %>%
          select(matches("^(Lh|beta|alpha|sigma|var)", ignore.case = TRUE))
      })

      dfs_mcmc <- lapply(df_list, as.mcmc)
      dfs_mcmclist <- as.mcmc.list(dfs_mcmc)
      gelman_cols <- param_cols[!grepl("^R\\.", param_cols)]
      gelman_cols <- gelman_cols[which(gelman_cols %in% colnames(combined_var))]
      gelman_result <- gelman.diag(dfs_mcmclist)[[1]][gelman_cols, 1]
      results$Gelman.psrf <- gelman_result[rownames(results)]
    }
  }

  # Round numeric values
  results[] <- lapply(results, function(x) if (is.numeric(x)) signif(x, 4) else x)
  return(results)
}

print_aligned_table <- function(df, digits = 4) {

  # Format numeric columns
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) {
      format(signif(x, digits), justify = "right")
    } else {
      format(x, justify = "left")
    }
  })

  # Ensure rownames are included and aligned
  df_out <- cbind(Parameter = rownames(df), df)
  rownames(df_out) <- NULL

  print(df_out, row.names = FALSE, right = TRUE)
  df_out
}

# Read the log marginal likelihood from a BayesTraits .Stones.txt file (stepping-stone
# sampler output). Returns NA if the file doesn't exist or wasn't found (e.g. the run
# didn't use the "Stones" command).
read_log_marginal_likelihood <- function(stones_file) {
  if (!file.exists(stones_file)) return(NA_real_)
  lines <- readLines(stones_file, warn = FALSE)
  ml_line <- grep("^Log marginal likelihood:", lines, value = TRUE)
  if (length(ml_line) == 0) return(NA_real_)
  as.numeric(trimws(sub("^Log marginal likelihood:", "", ml_line[length(ml_line)])))
}

# Finds base name of logfiles to name summary plot folder
longest_common_substring <- function(vec) {
  # Start with the first string
  s1 <- vec[1]
  len <- nchar(s1)

  # Try all substrings of s1 from longest to shortest
  for (l in len:1) {
    for (start in 1:(len - l + 1)) {
      sub <- substr(s1, start, start + l - 1)
      if (all(grepl(sub, vec, fixed = TRUE))) {
        return(sub)
      }
    }
  }
  return("")  # no common substring
}


# SINGLE LOG FILE PLOTS ----------------------------------------------------
plot_single_log <- function(file_path, homedir) {
  print(paste("plotting", file_path), quote = FALSE)
  df <- read_bt_log(file_path)
  outName <- sub("\\..*", "", basename(file_path))
  outdir <- file.path(homedir, paste0(outName, "_results"))
  dir.create(outdir, showWarnings = FALSE)

  cols <- get_cols(df)
  col_colors <- c(Lh = "grey40", Beta = "firebrick4", Alpha = "royalblue4", Sigma = "seagreen4", q = "mediumpurple3", R. = "orange3")

  for (kw in names(cols)) {
    param_cols <- cols[kw]
    col_color <- col_colors[kw]

    for (param in param_cols) {
      # Trace plot
      p_trace <- ggplot(df, aes(x = df[, 1], y = .data[[param]])) +
        geom_line(color = col_color, linewidth = 0.1) +
        theme_light() +
        labs(x = "Iteration", y = param, title = paste("Trace plot for", param))
      ggsave(file.path(outdir, paste0(outName, param, "_Trace.png")), p_trace, width = 8, height = 6)

      # Density plot
      p_density <- ggplot(df, aes(x = .data[[param]])) +
        geom_density(fill = col_color, alpha = 0.3) +
        theme_light() +
        labs(x = param, title = paste("Density plot for", param))

      # Sigma xlim adjustment
      if (kw == "Sigma") {
        ci_vals <- unlist(ci(df[[param]]))
        xl <- c(ci_vals[2] - (2 * sd(unlist(df[[param]]))), ci_vals[3] + (2 * sd(unlist(df[[param]]))))
        p_density <- p_density + xlim(xl)
      }
      ggsave(file.path(outdir, paste0(outName, param, "_Density.png")), p_density, width = 8, height = 6)
    }
  }

  # Build parameter list (same as original logic)
  all_cols <- unlist(cols[c("Lh", "Beta", "Alpha", "Sigma", "R.", "q")], use.names = FALSE)

  # Create results table
  results <- make_results_table(df, all_cols, diagnostics = "single")

  # Prepend the log marginal likelihood from the stepping-stone sampler output
  # (.Stones.txt), if the run produced one
  stones_file <- file.path(dirname(file_path), paste0(outName, ".Stones.txt"))
  log_marg_lik <- read_log_marginal_likelihood(stones_file)
  if (!is.na(log_marg_lik)) {
    lml_row <- results[1, ]
    lml_row[] <- NA
    rownames(lml_row) <- "LnL_stonesEST"
    lml_row$Mean <- signif(log_marg_lik, 4)
    results <- rbind(lml_row, results)
  }

  # Write to file (row names -> explicit "Model_Output" column, since write.table
  # otherwise writes row names as an unlabeled first column)
  results_out <- cbind(Model_Output = rownames(results), results)
  write.table(results_out,
              file = file.path(outdir, paste0(outName, "ResultTbl.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
}


# MULTI-RUN SUMMARY PLOTS ----------------------------------------------------
plot_multi_run_summary <- function(log_files, homedir) {

  col_colors <- c(Lh = "grey40", Beta = "firebrick4", Alpha = "royalblue4", Sigma = "seagreen4", q = "mediumpurple3", R. = "orange3")

  # Group files by signature
  signatures <- sapply(log_files, extract_relevant)
  matches <- split(log_files, signatures)
  matches <- matches[sapply(matches, length) > 1]

  cat("Found", length(matches), "Logfile groups with identical input, parameters, & settings.\n")

  for (group_idx in seq_along(matches)) {
    cat("Processing match group:", group_idx, "\n")

    group_files <- matches[[group_idx]]

    if (length(group_files) == 1) {

      plot_single_log(group_files[1], homedir)
      print(paste("only one logfile found in match", group_idx, "\n "))

    } else {
      outName_base <- longest_common_substring(matches[[group_idx]])
      outName <- paste0("AllRuns_", outName_base)
      outdir <- file.path(homedir, paste0(outName, "results"))
      dir.create(outdir, showWarnings = FALSE)

      # Read all runs
      all_runs <- lapply(seq_along(group_files), function(idx) {
        df <- read_bt_log(group_files[idx])
        df$num_runID <- factor(idx)
        df
      })
      combined <- do.call(rbind, all_runs)

      cols <- get_cols(combined)

      # Define individual run colors
      ind_run_colors <- c(
        Lh = "grey80",
        Beta = "indianred2",
        Alpha = "steelblue2",
        Sigma = "palegreen2",
        q = "darkorchid1",
        R. = "sienna1"
      )

      for (kw in keywords) {
        param_cols <- cols[[kw]]
        if (length(param_cols) == 0) next
        col_color <- col_colors[[kw]]  # summary density color

        for (param in param_cols) {
          # Trace plot
          p_trace <- ggplot(combined, aes(x = combined[, 1], y = .data[[param]], color = num_runID, group = num_runID)) +
            geom_line(linewidth = 0.15, alpha = 0.5) +
            scale_color_carto_d(palette = "ag_Sunset") +
            theme_light() +
            labs(x = "Iteration", y = param, title = paste("Trace plot for", param))
          ggsave(file.path(outdir, paste0(outName, param, "_Trace.png")), p_trace, width = 10, height = 6)

          # Density plot: individual + summary
          p_density <- ggplot() + theme_light() +
            labs(x = param, title = paste("Density plot for", param))

          # Individual run densities with custom colors
          for (run_idx in levels(combined$num_runID)) {
            run_data <- combined[combined$num_runID == run_idx, param, drop = FALSE]  # keep as data frame
            p_density <- p_density +
              geom_density(data = run_data,
                           aes(x = .data[[param]]),
                           color = scales::alpha(col_colors[[kw]], 0.7),
                           linewidth = 0.25,
                           alpha = 0.1,
                           fill = ind_run_colors[[kw]]
              )
          }

          # Summary density over all runs (semi-opaque)
          summary_vals <- combined[[param]]
          p_density <- p_density +
            geom_density(aes(x = summary_vals),
                         fill = col_color,
                         alpha = 0.4,
                         color = "grey25",
                         linewidth = 0.5)

          # Sigma xlim adjustment
          if (kw == "Sigma") {
            ci_vals <- unlist(ci(summary_vals))
            xl <- c(ci_vals[2] - (2 * sd(unlist(combined[[param]]))), ci_vals[3] + (2 * sd(unlist(combined[[param]]))))
            p_density <- p_density + xlim(xl)
          }

          ggsave(file.path(outdir, paste0(outName, param, "_Density.png")), p_density, width = 10, height = 6)
        }
      }

      # Combine all parameter columns
      all_cols <- unlist(cols[c("Lh", "Beta", "Alpha", "Sigma", "R.", "q")], use.names = FALSE)

      # Results based on combined chains
      results <- make_results_table(combined, all_cols, diagnostics = "multi")

      cat("\n===== Multi-run summary:", outName, "=====\n")
      results_align <- print_aligned_table(results)
      cat("========================================\n\n")

      # Write to file
      write.table(results_align,
                  file = file.path(outdir, paste0(outName, "ResultTbl.txt")),
                  quote = FALSE, sep = "\t")

    } # closes if/else testing if only one logfile was located

  } # closes loop over logfile match sets in WD

} # closes multirun plot function


# MAIN EXECUTION (CLI) ----------------------------------------------------
# Only runs when this file is executed directly (e.g. `Rscript bt_logPlotting.R`),
# not when it's source()'d (e.g. from main.R) just to load the functions above.
# sys.nframe() == 0 only at the true top level of a directly-run script.
if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)
  homedir <- getwd()
  plotsingles <- TRUE

  if (length(args) > 0) {
    if ("plotsingles" %in% args) {
      plotsingles <- TRUE
      args <- args[!args %in% "plotsingles"]
    }
    if (length(args) > 0) {
      all_files <- args
    } else {
      all_files <- list.files()
    }
  } else {
    all_files <- list.files()
  }

  all_log_files <- all_files[grep("\\.Log\\.txt$", all_files)]

  # Single log files
  if (plotsingles == TRUE) {
    for (file_path in all_log_files) {
      plot_single_log(file_path, homedir)
    }
  }

  # Multi-run summary (only if multiple runs exist)
  plot_multi_run_summary(all_log_files, homedir)

}
