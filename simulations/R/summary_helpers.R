# sequence of deltas and labels for plotting
delta_shift <- seq(-2, 2, by = 0.5)
label_maker <- c(
  "-2.0" = "Shift down by 2.0 SD units",
  "-1.5" = "Shift down by 1.5 SD units",
  "-1.0" = "Shift down by 1.0 SD units",
  "-0.5" = "Shift down by 0.5 SD units",
  "0" = "No shift (0 std. units)",
  "0.5" = "Shift up by 0.5 SD units",
  "1.0" = "Shift up by 1.0 SD units",
  "1.5" = "Shift up by 1.5 SD units",
  "2.0" = "Shift up by 2.0 SD units",
  "ipcw_eif" = "Augmented EIF (efficient with ML)",
  "ipcw_loss" = "Augmented Loss (efficient with GLMs)"
)


# helper function to get 5 most recent files for efficient estimator simulations
get_recent_files <- function(file_names, num_files = 10) {
  file_dates <- as.Date(str_extract(file_names, "2019......"))
  recent_files <- sort(file_dates, decreasing = TRUE)[seq_len(num_files)]
  recent_files_regex <- paste(recent_files, collapse = '|')
  which_recent <- str_detect(file_names, recent_files_regex)
  return(file_names[which_recent])
}
