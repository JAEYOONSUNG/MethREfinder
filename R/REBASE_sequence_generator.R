#' Ambiguous DNA values for nucleotide matching
#' @export
ambiguous_dna_values <- list(
  A = "A", C = "C", G = "G", T = "T",
  M = "AC", R = "AG", W = "AT", S = "CG",
  Y = "CT", K = "GT", V = "ACG", H = "ACT",
  D = "AGT", B = "CGT", N = "ACGT"
)

#' Null coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Check ambiguous match between two sequences
#' @export
ambiguous_match <- function(seq1, seq2) {
  if (nchar(seq1) != nchar(seq2)) return(FALSE)
  for (i in seq_len(nchar(seq1))) {
    char1 <- substr(seq1, i, i)
    char2 <- substr(seq2, i, i)
    allowed1 <- ambiguous_dna_values[[char1]] %||% char1
    allowed2 <- ambiguous_dna_values[[char2]] %||% char2
    if (!any(strsplit(allowed1, "")[[1]] %in% strsplit(allowed2, "")[[1]])) {
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Reverse complement of a DNA sequence
#' @export
rev_comp <- function(dna) {
  iupac_complement <- c(
    "A" = "T", "T" = "A", "G" = "C", "C" = "G",
    "R" = "Y", "Y" = "R", "S" = "S", "W" = "W",
    "K" = "M", "M" = "K", "B" = "V", "V" = "B",
    "D" = "H", "H" = "D", "N" = "N"
  )
  rev_dna <- rev(strsplit(dna, "")[[1]])
  dna_rc <- sapply(rev_dna, function(base) iupac_complement[base])
  return(paste0(dna_rc, collapse = ""))
}

# ============================================
# FORWARD OVERLAP: window_suffix + RecSeq_prefix
# ============================================

#' Find longest FORWARD overlap (window suffix matches RecSeq prefix)
#' @export
longest_overlap_forward <- function(window, recseq, mod_pos) {
  window <- trimws(window)
  recseq <- trimws(recseq)
  if (is.na(window) || is.na(recseq) || window == "" || recseq == "" || is.na(mod_pos)) {
    return(0)
  }
  
  # Remove trailing N from window
  w_chars <- strsplit(window, "")[[1]]
  trailing_N <- 0
  for (i in rev(seq_along(w_chars))) {
    if (w_chars[i] == "N") trailing_N <- trailing_N + 1
    else break
  }
  w_clean <- substr(window, 1, nchar(window) - trailing_N)
  if (nchar(w_clean) == 0) return(0)
  
  mod_pos_adj <- min(mod_pos, nchar(w_clean))
  max_len <- 0
  limit <- min(nchar(w_clean), nchar(recseq))
  
  for (k in seq_len(limit)) {
    w_suffix <- substr(w_clean, nchar(w_clean) - k + 1, nchar(w_clean))
    r_prefix <- substr(recseq, 1, k)
    
    if (ambiguous_match(w_suffix, r_prefix)) {
      overlap_start <- nchar(w_clean) - k + 1
      if (overlap_start <= mod_pos_adj) {
        max_len <- k
      }
    }
  }
  return(max_len)
}

#' Concatenate for FORWARD overlap
#' @export
concat_forward <- function(window, recseq, mod_pos) {
  k <- longest_overlap_forward(window, recseq, mod_pos)
  if (k == 0) return(NA_character_)
  
  # Remove trailing N
  w_chars <- strsplit(window, "")[[1]]
  trailing_N <- 0
  for (i in rev(seq_along(w_chars))) {
    if (w_chars[i] == "N") trailing_N <- trailing_N + 1
    else break
  }
  w_clean <- substr(window, 1, nchar(window) - trailing_N)
  
  # Resolve overlap (prefer more specific base)
  w_suffix <- substr(w_clean, nchar(w_clean) - k + 1, nchar(w_clean))
  r_prefix <- substr(recseq, 1, k)
  
  resolved <- ""
  for (i in 1:k) {
    cw <- substr(w_suffix, i, i)
    cr <- substr(r_prefix, i, i)
    allowed_w <- ambiguous_dna_values[[cw]] %||% cw
    if (cr %in% strsplit(allowed_w, "")[[1]]) {
      resolved <- paste0(resolved, cr)
    } else {
      resolved <- paste0(resolved, cw)
    }
  }
  
  pre <- substr(w_clean, 1, nchar(w_clean) - k)
  post <- substr(recseq, k + 1, nchar(recseq))
  return(paste0(pre, resolved, post))
}

# ============================================
# BACKWARD OVERLAP: RecSeq_suffix + window_prefix
# ============================================

#' Find longest BACKWARD overlap (RecSeq suffix matches window prefix)
#' @export
longest_overlap_backward <- function(window, recseq, mod_pos) {
  window <- trimws(window)
  recseq <- trimws(recseq)
  if (is.na(window) || is.na(recseq) || window == "" || recseq == "" || is.na(mod_pos)) {
    return(0)
  }
  
  # Remove leading N from window
  w_chars <- strsplit(window, "")[[1]]
  leading_N <- 0
  for (i in seq_along(w_chars)) {
    if (w_chars[i] == "N") leading_N <- leading_N + 1
    else break
  }
  w_clean <- substr(window, leading_N + 1, nchar(window))
  if (nchar(w_clean) == 0) return(0)
  
  # Adjust mod_pos for leading N removal
  mod_pos_in_clean <- mod_pos - leading_N
  if (mod_pos_in_clean < 1) return(0)
  
  max_len <- 0
  limit <- min(nchar(w_clean), nchar(recseq))
  
  for (k in seq_len(limit)) {
    r_suffix <- substr(recseq, nchar(recseq) - k + 1, nchar(recseq))
    w_prefix <- substr(w_clean, 1, k)
    
    if (ambiguous_match(r_suffix, w_prefix)) {
      # mod_pos must be within or after the overlap region
      if (mod_pos_in_clean <= k) {
        max_len <- k
      }
    }
  }
  return(max_len)
}

#' Concatenate for BACKWARD overlap
#' @export
concat_backward <- function(window, recseq, mod_pos) {
  k <- longest_overlap_backward(window, recseq, mod_pos)
  if (k == 0) return(NA_character_)
  
  # Remove leading N
  w_chars <- strsplit(window, "")[[1]]
  leading_N <- 0
  for (i in seq_along(w_chars)) {
    if (w_chars[i] == "N") leading_N <- leading_N + 1
    else break
  }
  w_clean <- substr(window, leading_N + 1, nchar(window))
  
  # Resolve overlap
  r_suffix <- substr(recseq, nchar(recseq) - k + 1, nchar(recseq))
  w_prefix <- substr(w_clean, 1, k)
  
  resolved <- ""
  for (i in 1:k) {
    cr <- substr(r_suffix, i, i)
    cw <- substr(w_prefix, i, i)
    allowed_w <- ambiguous_dna_values[[cw]] %||% cw
    allowed_r <- ambiguous_dna_values[[cr]] %||% cr
    # Prefer more specific
    if (nchar(allowed_r) <= nchar(allowed_w)) {
      resolved <- paste0(resolved, cr)
    } else {
      resolved <- paste0(resolved, cw)
    }
  }
  
  pre <- substr(recseq, 1, nchar(recseq) - k)
  post <- substr(w_clean, k + 1, nchar(w_clean))
  return(paste0(pre, resolved, post))
}

# ============================================
# HELPER FUNCTIONS
# ============================================

#' Filter methylation sensitivity data for blocked enzymes
#' @export
filter_meth_blocked <- function(methylation_sensitivity) {
  methylation_sensitivity %>%
    dplyr::filter(cleavage == "blocked")
}

#' Filter methylation sensitivity data by modification type
#' @export
filter_meth_by_mod_type <- function(methylation_sensitivity, mod_type) {
  methylation_sensitivity %>%
    dplyr::filter(modification == mod_type)
}

#' Cross join and filter NA values
#' @export
cross_join_and_filter_na <- function(df_filtered, df_meth_blocked) {
  tidyr::crossing(df_filtered, df_meth_blocked) %>%
    dplyr::filter(!is.na(mod_containing_seq) & !is.na(forward_RecSeq))
}

#' Check if extended sequence is compatible with target sequence
#' @export
check_extended_compatible <- function(extended_seq, target_seq) {
  L_ext <- nchar(extended_seq)
  L_tar <- nchar(target_seq)
  
  if (L_ext < L_tar) return(FALSE)
  if (L_ext == L_tar) return(ambiguous_match(extended_seq, target_seq))
  
  # Check if target can be found within extended
  for (start in 1:(L_ext - L_tar + 1)) {
    candidate <- substr(extended_seq, start, start + L_tar - 1)
    if (ambiguous_match(candidate, target_seq)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Extend concatenated sequence to match target length (FORWARD direction)
#' @export
extend_to_target_forward <- function(concat_seq, target_seq) {
  if (is.na(concat_seq)) return(NA_character_)
  
  L_concat <- nchar(concat_seq)
  L_target <- nchar(target_seq)
  
  # Find overlap: concat_suffix matches target_prefix
  max_overlap <- 0
  for (k in seq_len(min(L_concat, L_target))) {
    c_suffix <- substr(concat_seq, L_concat - k + 1, L_concat)
    t_prefix <- substr(target_seq, 1, k)
    if (ambiguous_match(c_suffix, t_prefix)) {
      max_overlap <- k
    }
  }
  
  if (max_overlap == 0) return(NA_character_)
  
  extended <- paste0(concat_seq, substr(target_seq, max_overlap + 1, L_target))
  return(extended)
}

#' Extend concatenated sequence to match target length (BACKWARD direction)
#' @export
extend_to_target_backward <- function(concat_seq, target_seq) {
  if (is.na(concat_seq)) return(NA_character_)
  
  L_concat <- nchar(concat_seq)
  L_target <- nchar(target_seq)
  
  # Find overlap: concat_prefix matches target_suffix
  max_overlap <- 0
  for (k in seq_len(min(L_concat, L_target))) {
    c_prefix <- substr(concat_seq, 1, k)
    t_suffix <- substr(target_seq, L_target - k + 1, L_target)
    if (ambiguous_match(c_prefix, t_suffix)) {
      max_overlap <- k
    }
  }
  
  if (max_overlap == 0) return(NA_character_)
  
  extended <- paste0(substr(target_seq, 1, L_target - max_overlap), concat_seq)
  return(extended)
}

#' Validate mod_position: user's methylated position must match REBASE mod_position
#' @export
validate_mod_position <- function(extended_seq, target_seq, recseq, 
                                   user_mod_position, rebase_mod_position) {
  if (is.na(extended_seq) || is.na(recseq) || 
      is.na(user_mod_position) || is.na(rebase_mod_position)) {
    return(FALSE)
  }
  
  L_ext <- nchar(extended_seq)
  L_rec <- nchar(recseq)
  L_tar <- nchar(target_seq)
  
  # Step 1: Find where target aligns in extended
  target_start <- NA
  for (s in 1:(L_ext - L_tar + 1)) {
    if (ambiguous_match(substr(extended_seq, s, s + L_tar - 1), target_seq)) {
      target_start <- s
      break
    }
  }
  if (is.na(target_start)) return(FALSE)
  
  # Step 2: User's mod position in extended coordinates
  user_pos_in_ext <- target_start + user_mod_position - 1
  
  # Step 3: Find where RecSeq aligns in extended
  for (r in 1:(L_ext - L_rec + 1)) {
    if (ambiguous_match(substr(extended_seq, r, r + L_rec - 1), recseq)) {
      # Step 4: Check if user's position falls within RecSeq
      if (user_pos_in_ext >= r && user_pos_in_ext <= r + L_rec - 1) {
        # Position within RecSeq
        pos_in_recseq <- user_pos_in_ext - r + 1
        # Step 5: Must match REBASE mod_position
        if (pos_in_recseq == rebase_mod_position) {
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

# ============================================
# MAIN FUNCTION
# ============================================

#' Match enzyme sequences with target sequence
#'
#' @description
#' Matches enzymes against a target sequence using both FORWARD and BACKWARD overlap.
#' - FORWARD: window_suffix + RecSeq_prefix (RecSeq extends after window)
#' - BACKWARD: RecSeq_suffix + window_prefix (RecSeq extends before window)
#'
#' @param mod_type Character string (e.g., "4mC", "5mC", "6mA").
#' @param target_seq Character string (target DNA sequence).
#' @return Data frame with matched enzyme sequences.
#' @export
match_enzyme_sequences <- function(mod_type, target_seq) {
  df_meth_blocked <- filter_meth_blocked(methylation_sensitivity)
  df_meth_filtered <- filter_meth_by_mod_type(df_meth_blocked, mod_type)
  df_joined <- cross_join_and_filter_na(methylasensitive_window, df_meth_filtered)
  
  # Get user's original mod_position
  user_mod_position <- methylasensitive_window$start[1] + methylasensitive_window$mod_pos[1] - 1
  
  results <- list()
  
  # ============================================
  # METHOD 1: FORWARD overlap
  # ============================================
  df_fwd <- df_joined %>%
    dplyr::mutate(
      overlap_fwd = purrr::pmap_int(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        longest_overlap_forward
      ),
      concat_fwd = purrr::pmap_chr(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        concat_forward
      )
    ) %>%
    dplyr::filter(overlap_fwd > 0, !is.na(concat_fwd))
  
  if (nrow(df_fwd) > 0) {
    df_fwd <- df_fwd %>%
      dplyr::mutate(
        extended_seq = purrr::map_chr(concat_fwd, function(seq) {
          ext <- extend_to_target_forward(seq, target_seq)
          if (!is.na(ext) && check_extended_compatible(ext, target_seq)) ext
          else NA_character_
        }),
        match_type = "forward"
      ) %>%
      dplyr::filter(!is.na(extended_seq))
    
    if (nrow(df_fwd) > 0) {
      df_fwd <- df_fwd %>%
        dplyr::mutate(
          valid = purrr::pmap_lgl(
            list(extended_seq, forward_RecSeq, mod_position),
            function(ext, rec, rebase_pos) {
              validate_mod_position(ext, target_seq, rec, user_mod_position, rebase_pos)
            }
          )
        ) %>%
        dplyr::filter(valid) %>%
        dplyr::select(-valid, -overlap_fwd, -concat_fwd)
      
      results$forward <- df_fwd
    }
  }
  
  # ============================================
  # METHOD 2: BACKWARD overlap
  # ============================================
  df_bwd <- df_joined %>%
    dplyr::mutate(
      overlap_bwd = purrr::pmap_int(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        longest_overlap_backward
      ),
      concat_bwd = purrr::pmap_chr(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        concat_backward
      )
    ) %>%
    dplyr::filter(overlap_bwd > 0, !is.na(concat_bwd))
  
  if (nrow(df_bwd) > 0) {
    df_bwd <- df_bwd %>%
      dplyr::mutate(
        extended_seq = purrr::map_chr(concat_bwd, function(seq) {
          ext <- extend_to_target_backward(seq, target_seq)
          if (!is.na(ext) && check_extended_compatible(ext, target_seq)) ext
          else NA_character_
        }),
        match_type = "backward"
      ) %>%
      dplyr::filter(!is.na(extended_seq))
    
    if (nrow(df_bwd) > 0) {
      df_bwd <- df_bwd %>%
        dplyr::mutate(
          valid = purrr::pmap_lgl(
            list(extended_seq, forward_RecSeq, mod_position),
            function(ext, rec, rebase_pos) {
              validate_mod_position(ext, target_seq, rec, user_mod_position, rebase_pos)
            }
          )
        ) %>%
        dplyr::filter(valid) %>%
        dplyr::select(-valid, -overlap_bwd, -concat_bwd)
      
      results$backward <- df_bwd
    }
  }
  
  # ============================================
  # Combine results
  # ============================================
  if (length(results) > 0) {
    result_df <- dplyr::bind_rows(results) %>%
      dplyr::distinct(Enzyme, forward_RecSeq, extended_seq, .keep_all = TRUE)
    rownames(result_df) <- NULL
  } else {
    result_df <- data.frame(
      Enzyme = character(),
      forward_RecSeq = character(),
      extended_seq = character(),
      match_type = character(),
      stringsAsFactors = FALSE
    )
  }
  
  result_df <<- result_df
  return(result_df)
}
