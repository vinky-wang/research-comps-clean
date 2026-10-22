# Run from research-comps-clean/extension_test using the existing R 4.2 launcher
source('test_helpers.R')
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) sub('^--phase=', '', args[1]) else 'all'
stopifnot(length(args) <= 1L, phase %in% c('all', 'pilot', 'stability', 'report'))
verify_snapshot()
tryCatch({
  if (phase %in% c('all','pilot')) run_pilot()
  if (phase %in% c('all','stability')) run_stability()
}, finally = write_reports())
