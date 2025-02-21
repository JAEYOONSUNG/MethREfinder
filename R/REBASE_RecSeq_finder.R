#' Extract All Subwindows Clamped (No Skipping)
#'
#' @description
#' Enumerates all subwindows of specified sizes around a modified position, with clamping to ensure valid positions.
#'
#' @param sequence A character string representing the DNA sequence (A,T,G,C,N), uppercase.
#' @param mod_position Integer. The 1-based index of the modified base.
#' @param mod_type Character. A label representing the modification type (e.g., "6mA" or c("4mC","5mC","6mA")). Default is "manual".
#' @param direction Character. Either "forward" or "reverse". Determines whether to reverse complement the window.
#' @param window_sizes Integer vector. Specifies the window sizes to extract (e.g., c(4,5,6)).
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{window_size}: Length of the window.
#'   \item \code{start}, \code{end}: Start and end positions of the window.
#'   \item \code{mod_containing_seq}: Extracted sequence window (padded with 'N' if needed).
#'   \item \code{mod_pos}: 1-based position of the modification within the window.
#'   \item \code{direction}: Orientation ("forward" or "reverse").
#'   \item \code{mod_type}: Modification type label.
#' }
#'
#' @examples
#' df_result <- extract_all_windows_RecSeq(
#'   sequence = "CGAANNNNNNNTARC",
#'   mod_position = 4,
#'   mod_type = "6mA",
#'   direction = "forward",
#'   window_sizes = c(4,5,6)
#' )
#' print(df_result)
#'
#' @export

extract_all_windows_RecSeq <- function(
    sequence,
    mod_position,
    mod_type     = "manual",
    direction    = c("forward","reverse"),
    window_sizes = c(4,5,6)
) {
  direction <- match.arg(direction)
  nseq      <- nchar(sequence)

  rev_comp <- function(dna) {
    # IUPAC 코드 변환 테이블 (1:1 치환)
    iupac_complement <- c(
      "A" = "T", "T" = "A", "G" = "C", "C" = "G",
      "R" = "Y", "Y" = "R", "S" = "S", "W" = "W",
      "K" = "M", "M" = "K", "B" = "V", "V" = "B",
      "D" = "H", "H" = "D", "N" = "N"
    )

    # 서열 뒤집기
    rev_dna <- rev(strsplit(dna, "")[[1]])

    # Reverse Complement 변환 적용
    dna_rc <- sapply(rev_dna, function(base) iupac_complement[base])

    # 문자열로 반환
    return(paste0(dna_rc, collapse=""))
  }

  extract_exact_window <- function(seq, start_idx, win_size) {
    end_idx <- start_idx + win_size - 1
    prefixN <- 0
    suffixN <- 0

    if (start_idx < 1) {
      prefixN   <- 1 - start_idx
      start_idx <- 1
    }
    if (end_idx > nseq) {
      suffixN <- end_idx - nseq
      end_idx <- nseq
    }
    core_sub <- ""
    if (start_idx <= end_idx && start_idx <= nseq && end_idx >= 1) {
      core_sub <- substr(seq, start_idx, end_idx)
    }
    left_part  <- paste0(rep("N", prefixN), collapse="")
    right_part <- paste0(rep("N", suffixN), collapse="")
    raw_window <- paste0(left_part, core_sub, right_part)

    lw <- nchar(raw_window)
    if (lw < win_size) {
      needed <- win_size - lw
      raw_window <- paste0(raw_window, paste0(rep("N", needed), collapse=""))
    } else if (lw > win_size) {
      raw_window <- substr(raw_window, 1, win_size)
    }
    return(list(
      window_seq     = raw_window,
      prefixN_count  = prefixN
    ))
  }

  results <- list()

  for (ws in window_sizes) {
    s_min <- mod_position - (ws - 1)
    s_max <- mod_position
    starts_vec <- seq(s_min, s_max)

    row_list <- lapply(starts_vec, function(s0) {
      e0 <- s0 + ws - 1
      winfo <- extract_exact_window(sequence, s0, ws)
      wseq  <- winfo$window_seq
      pN    <- winfo$prefixN_count

      offset_fwd <- pN + (mod_position - s0) + 1

      if (direction == "reverse") {
        wseq <- rev_comp(wseq)
        offset_final <- ws - offset_fwd + 1
      } else {
        offset_final <- offset_fwd
      }

      if (offset_final < 1)  offset_final <- 1
      if (offset_final > ws) offset_final <- ws

      data.frame(
        window_size        = ws,
        start              = s0,
        end                = e0,
        mod_containing_seq = wseq,
        mod_pos            = offset_final,
        stringsAsFactors   = FALSE
      )
    })

    df_ws <- do.call(rbind, row_list)
    results[[as.character(ws)]] <- df_ws
  }

  methylasensitive_window <- do.call(rbind, results)
  rownames(methylasensitive_window) <- NULL
  methylasensitive_window$direction <- direction
  methylasensitive_window$mod_type  <- mod_type

  assign("methylasensitive_window", methylasensitive_window, envir = .GlobalEnv)

  results <- list() # Reset results after assigning the final result

  return(methylasensitive_window)
}
