#' Run the Full MethREfinder Pipeline
#'
#' This function executes the complete MethREfinder analysis pipeline, including:
#' 1. Fetching methylation sensitivity data from REBASE.
#' 2. Extracting recognition sequence windows from the provided sequence.
#' 3. Matching enzymes against the provided target sequence based on methylation sensitivity.
#'
#' @param target_seq Character. DNA sequence to analyze and match enzymes against.
#' @param mod_type Character. Modification type (e.g., \"6mA\", \"5mC\").
#' @param mod_position Integer. Position of the modified base within the recognition sequence.
#' @param direction Character. Direction for sequence extraction (\"forward\" or \"reverse\").
#' @param window_sizes Integer vector. Window sizes for sequence extraction.
#'
#' @return A list containing:
#' \describe{
#'   \item{methylation_table}{Data frame with methylation sensitivity data from REBASE.}
#'   \item{windows}{Extracted sequence windows.}
#'   \item{matched_enzymes}{Matched enzymes for the given target sequence.}
#' }
#' @examples
#' \dontrun{
#' result <- run_MethREfinder(
#'   target_seq = \"CGAANNNNNNNTARC\",
#'   mod_type = \"6mA\",
#'   mod_position = 4,
#'   direction = \"forward\",
#'   window_sizes = c(4, 5, 6)
#' )
#' print(result$methylation_table)
#' print(result$windows)
#' print(result$matched_enzymes)
#' }
#' @export
run_MethREfinder <- function(target_seq, mod_type, mod_position, direction, window_sizes) {
  # 1. Fetch REBASE methylation sensitivity data
  cat("Step 1. Fetching REBASE methylation sensitivity data...\n")
  methylation_table <- fetch_rebase_methylation()
  cat("Step 1. Complete!\n\n")
  
  # 2. Extract windows from sequence
  cat("Step 2. Extracting sequence windows...\n")
  windows <- extract_all_windows_RecSeq(
    sequence = target_seq,
    mod_position = mod_position,
    mod_type = mod_type,
    direction = direction,
    window_sizes = window_sizes
  )
  cat("Step 2. Complete!\n\n")
  
  # 3. Match enzymes against target sequence
  cat("Step 3. Matching enzymes against target sequence...\n")
  matched_enzymes <- match_enzyme_sequences(
    mod_type = mod_type,
    target_seq = target_seq
  )
  cat("Step 3. Complete!\n\n")
  
  # 결과를 리스트로 반환
  cat("All steps completed successfully!\n")
  return(list(
    methylation_table = methylation_table,
    windows = windows,
    matched_enzymes = matched_enzymes
  ))
}