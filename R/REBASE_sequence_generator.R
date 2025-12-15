#' Ambiguous DNA values for nucleotide matching
#'
#' @description
#' A list mapping each nucleotide (or ambiguous IUPAC code) to its possible base representations.
#'
#' @return A named list where each element is a string of acceptable bases.
#'
#' @examples
#' ambiguous_dna_values$N  # returns "ACGT"
#'
#' @export
ambiguous_dna_values <- list(
  A = "A", C = "C", G = "G", T = "T",
  M = "AC", R = "AG", W = "AT", S = "CG",
  Y = "CT", K = "GT", V = "ACG", H = "ACT",
  D = "AGT", B = "CGT", N = "ACGT"
)

#' Null coalescing operator
#' @param x First value
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Check ambiguous match between two sequences
#'
#' @description
#' This function checks whether two sequences match ambiguously based on IUPAC nucleotide codes.
#'
#' @param seq1 A character string representing the first DNA sequence.
#' @param seq2 A character string representing the second DNA sequence.
#'
#' @return Logical value (`TRUE` or `FALSE`) indicating if the sequences match ambiguously.
#'
#' @examples
#' ambiguous_match("A", "M") # TRUE (M represents A or C)
#' ambiguous_match("T", "G") # FALSE
#'
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
#'
#' @description
#' Computes the reverse complement of a DNA sequence considering IUPAC nucleotide codes.
#'
#' @param dna A character string representing a DNA sequence.
#'
#' @return A character string representing the reverse complement of the input sequence.
#'
#' @examples
#' rev_comp("ATGC") # returns "GCAT"
#'
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

#' Find the longest overlap between two strings with ambiguous matching
#'
#' @param a A character string (mod_containing_seq).
#' @param b A character string (forward_RecSeq).
#' @param mod_pos Integer. The 1-based position of the modification in window `a`.
#' @return An integer representing the length of the longest overlap.
#' @examples
#' longest_overlap("CGAAN", "GAAGA", 4)  # returns 2 (AA와 GA 오버랩)
#' @importFrom stringr str_detect fixed
#' @export
longest_overlap <- function(a, b, mod_pos) {
  a <- trimws(a)
  b <- trimws(b)
  if (is.na(a) || is.na(b) || a == "" || b == "" || is.na(mod_pos)) {
    return(0)
  }

  # 끝의 연속적인 N 제외
  a_chars <- strsplit(a, "")[[1]]
  trailing_N_count <- 0
  for (i in rev(seq_along(a_chars))) {
    if (a_chars[i] == "N") {
      trailing_N_count <- trailing_N_count + 1
    } else {
      break
    }
  }
  a_no_trailing_N <- substr(a, 1, nchar(a) - trailing_N_count)
  mod_pos_adjusted <- min(mod_pos, nchar(a_no_trailing_N))

  max_len <- 0
  limit <- min(nchar(a_no_trailing_N), nchar(b))

  for (k in seq_len(limit)) {
    overlap_match <- TRUE
    for (i in seq_len(k)) {
      char_a <- substr(a_no_trailing_N, nchar(a_no_trailing_N) - k + i, nchar(a_no_trailing_N) - k + i)
      char_b <- substr(b, i, i)
      char_a <- toupper(char_a)
      char_b <- toupper(char_b)

      allowed_a <- ambiguous_dna_values[[char_a]]
      allowed_b <- ambiguous_dna_values[[char_b]]
      if (is.null(allowed_a)) allowed_a <- char_a
      if (is.null(allowed_b)) allowed_b <- char_b

      # Check if there's any overlap between allowed bases
      set_a <- strsplit(allowed_a, "")[[1]]
      set_b <- strsplit(allowed_b, "")[[1]]
      if (length(intersect(set_a, set_b)) == 0) {
        overlap_match <- FALSE
        break
      }
    }
    if (overlap_match) {
      overlap_start_in_a <- nchar(a_no_trailing_N) - k + 1
      if (overlap_start_in_a <= mod_pos_adjusted) {
        max_len <- k
      }
    }
  }
  return(max_len)
}

#' Concatenate two strings by removing duplicated ambiguous overlap
#'
#' @param a Character string (mod_containing_seq).
#' @param b Character string (forward_RecSeq).
#' @param mod_pos Integer. The 1-based position of the modification in window `a`.
#' @return A concatenated character string.
#' @examples
#' concat_with_overlap("CGAAN", "GAAGA", 4)  # returns "CGAAGA"
#' @export
concat_with_overlap <- function(a, b, mod_pos) {
  k <- longest_overlap(a, b, mod_pos)
  if (k == 0) {
    return(paste0(a, b))
  }

  # 끝의 N 제외한 a의 유효 부분
  a_chars <- strsplit(a, "")[[1]]
  trailing_N_count <- 0
  for (i in rev(seq_along(a_chars))) {
    if (a_chars[i] == "N") {
      trailing_N_count <- trailing_N_count + 1
    } else {
      break
    }
  }
  a_no_trailing_N <- substr(a, 1, nchar(a) - trailing_N_count)

  # 오버랩 부분 추출
  overlap_a <- substr(a_no_trailing_N, nchar(a_no_trailing_N) - k + 1, nchar(a_no_trailing_N))
  overlap_b <- substr(b, 1, k)

  resolved_overlap <- ""
  for (i in 1:k) {
    char_a <- substr(overlap_a, i, i)
    char_b <- substr(overlap_b, i, i)
    allowed_a <- ambiguous_dna_values[[char_a]]
    if (is.null(allowed_a)) allowed_a <- char_a

    if (char_b %in% strsplit(allowed_a, "")[[1]]) {
      resolved_overlap <- paste0(resolved_overlap, char_b)
    } else {
      resolved_overlap <- paste0(resolved_overlap, char_a)
    }
  }

  pre_overlap <- substr(a_no_trailing_N, 1, nchar(a_no_trailing_N) - k)
  post_overlap <- substr(b, k + 1, nchar(b))
  return(paste0(pre_overlap, resolved_overlap, post_overlap))
}

#' Filter methylation sensitivity data for blocked enzymes
#'
#' @description
#' Filters methylation sensitivity data to retain only enzymes whose cleavage is blocked.
#'
#' @param methylation_sensitivity A data frame from REBASE containing enzyme information.
#'
#' @return A filtered data frame with only blocked enzymes.
#'
#' @examples
#' df_blocked <- filter_meth_blocked(methylation_sensitivity)
#'
#' @export
filter_meth_blocked <- function(methylation_sensitivity) {
  methylation_sensitivity %>%
    dplyr::filter(cleavage == "blocked")
}

#' Filter methylation sensitivity data by modification type
#'
#' @description
#' Filters methylation sensitivity data based on the specified modification type.
#'
#' @param methylation_sensitivity A data frame from REBASE containing enzyme information.
#' @param mod_type A character string specifying the modification type (e.g., "4mC").
#'
#' @return A filtered data frame containing only enzymes with the specified modification type.
#'
#' @examples
#' df_filtered <- filter_meth_by_mod_type(methylation_sensitivity, "4mC")
#'
#' @export
filter_meth_by_mod_type <- function(methylation_sensitivity, mod_type) {
  methylation_sensitivity %>%
    dplyr::filter(modification == mod_type)
}

#' Cross join and filter NA values
#'
#' @description
#' Performs a cross join between two data frames and removes rows with NA values in key columns.
#'
#' @param df_filtered A data frame containing sequence windows.
#' @param df_meth_blocked A data frame containing blocked enzyme data.
#'
#' @return A data frame with cross-joined data and removed NA values.
#'
#' @examples
#' df_joined <- cross_join_and_filter_na(df_filtered, df_meth_blocked)
#'
#' @importFrom tidyr crossing
#' @importFrom dplyr filter
#' @export
cross_join_and_filter_na <- function(df_filtered, df_meth_blocked) {
  tidyr::crossing(df_filtered, df_meth_blocked) %>%
    dplyr::filter(!is.na(mod_containing_seq) & !is.na(forward_RecSeq))
}

#' Check if modification position is within the overlap region
#'
#' @description
#' Determines whether the modification position falls within the overlap region between two sequences.
#'
#' @param context A character string representing the modification-containing sequence.
#' @param forward_RecSeq A character string representing the recognition sequence.
#' @param mod_pos An integer representing the modification position in the window `context`.
#' @param overlap_len An integer indicating the length of the overlap.
#'
#' @return A logical value (`TRUE` or `FALSE`) indicating whether the modification position is within the overlap region.
#'
#' @examples
#' check_mod_position_in_overlap("CCWGG", "CCAGG", 2, 3)
#'
#' @export
check_mod_position_in_overlap <- function(context, forward_RecSeq, mod_pos, overlap_len) {
  if (is.na(overlap_len) || is.na(mod_pos) || overlap_len == 0) {
    return(FALSE)
  }
  overlap_start_in_context <- nchar(context) - overlap_len + 1
  return(mod_pos >= overlap_start_in_context && mod_pos <= nchar(context))
}

#' Find intersection between concatenated sequence and target sequence
#'
#' @description
#' Identifies the longest matching substring between a concatenated recognition sequence and a target sequence.
#'
#' @param concatenated_RecSeq A character string representing the concatenated recognition sequence.
#' @param target_seq A character string representing the target sequence.
#'
#' @return The longest matching substring or `NA` if no match is found.
#'
#' @examples
#' find_intersection_in_target("CCWGG", "CCAGG")
#'
#' @export
find_intersection_in_target <- function(concatenated_RecSeq, target_seq) {
  if (is.na(concatenated_RecSeq) || is.na(target_seq)) {
    return(NA_character_)
  }
  L_rec <- nchar(concatenated_RecSeq)
  L_tar <- nchar(target_seq)
  for (len in L_rec:1) {
    for (start in 1:(L_rec - len + 1)) {
      candidate <- substr(concatenated_RecSeq, start, start + len - 1)
      if (grepl(candidate, target_seq, perl = TRUE)) {
        return(candidate)
      }
    }
  }
  return(NA_character_)
}

#' Check if extended sequence fully matches target sequence
#'
#' @description
#' Compares an extended sequence with a target sequence, considering ambiguous nucleotide matching.
#'
#' @param extended_seq A character string representing the extended sequence.
#' @param target_seq A character string representing the target sequence.
#'
#' @return A logical value (`TRUE` or `FALSE`) indicating whether the sequences match fully.
#'
#' @examples
#' check_full_extended_match("CCWGG", "CCAGG")
#'
#' @export
check_full_extended_match <- function(extended_seq, target_seq) {
  if (nchar(extended_seq) != nchar(target_seq)) {
    return(FALSE)
  }
  for (i in seq_len(nchar(extended_seq))) {
    char_ext <- substr(extended_seq, i, i)
    char_tar <- substr(target_seq, i, i)
    allowed_ext <- ambiguous_dna_values[[char_ext]]
    allowed_tar <- ambiguous_dna_values[[char_tar]]
    if (is.null(allowed_ext)) allowed_ext <- char_ext
    if (is.null(allowed_tar)) allowed_tar <- char_tar
    set_ext <- strsplit(allowed_ext, "")[[1]]
    set_tar <- strsplit(allowed_tar, "")[[1]]
    if (length(intersect(set_ext, set_tar)) == 0) {
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Extend sequence to match target sequence and filter
#'
#' @param df_out Data frame with concatenated sequences.
#' @param target_seq Character string representing the target DNA sequence.
#'
#' @return Filtered data frame with extended sequences.
#'
#' @examples
#' df_extended <- extend_with_target_seq_and_filter(df_out, "CCWGG")
#'
#' @importFrom dplyr mutate filter
#' @importFrom purrr pmap_chr
#' @export
extend_with_target_seq_and_filter <- function(df_out, target_seq) {
  df_out %>%
    dplyr::mutate(
      extended_seq = purrr::pmap_chr(
        list(concatenated_RecSeq, mod_pos, forward_RecSeq),
        function(seq, mod_pos_val, forward_seq) {
          if (is.na(seq) || is.na(mod_pos_val)) {
            return(NA_character_)
          }
          extended <- extend_to_target_seq(seq, target_seq)
          if (!is.na(extended) && check_full_extended_match(extended, target_seq)) {
            return(extended)
          }
          return(NA_character_)
        }
      )
    ) %>%
    dplyr::filter(!is.na(extended_seq))
}

#' Extend concatenated sequence to target sequence
#'
#' @description
#' Extends a concatenated recognition sequence to match the target sequence using ambiguous nucleotide matching.
#'
#' @param concatenated_RecSeq A character string representing the concatenated recognition sequence.
#' @param target_seq A character string representing the target DNA sequence.
#'
#' @return A character string representing the extended sequence or `NA` if no extension is possible.
#'
#' @examples
#' extend_to_target_seq("CCWGG", "CCAGG")
#'
#' @export
extend_to_target_seq <- function(concatenated_RecSeq, target_seq) {
  L_concat <- nchar(concatenated_RecSeq)
  L_target <- nchar(target_seq)
  max_overlap_len <- 0
  for (k in seq_len(min(L_concat, L_target))) {
    if (ambiguous_match(
      substr(concatenated_RecSeq, L_concat - k + 1, L_concat),
      substr(target_seq, 1, k)
    )) {
      max_overlap_len <- k
    }
  }
  if (max_overlap_len == 0) {
    return(NA_character_)
  }
  extended_seq <- paste0(
    concatenated_RecSeq,
    substr(target_seq, max_overlap_len + 1, L_target)
  )
  return(extended_seq)
}

#' Match enzyme sequences with target sequence
#'
#' @description
#' Main function to match enzymes against a target sequence.
#' Cross joins windows with REBASE data, then finds overlaps and extends to target sequence.
#' Uses mod_pos (window's modification position) for overlap calculation.
#'
#' @param mod_type Character string (e.g., "4mC", "5mC", "6mA").
#' @param target_seq Character string (target DNA sequence).
#' @return Data frame with matched enzyme sequences.
#' @importFrom dplyr mutate filter
#' @importFrom purrr pmap_int pmap_chr pmap_lgl
#' @export
match_enzyme_sequences <- function(mod_type, target_seq) {
  df_meth_blocked <- filter_meth_blocked(methylation_sensitivity)
  df_meth_filtered <- filter_meth_by_mod_type(df_meth_blocked, mod_type)
  df_joined <- cross_join_and_filter_na(methylasensitive_window, df_meth_filtered)

  # KEY FILTER: window의 mod_pos와 REBASE의 mod_position이 일치해야 함
  df_joined <- df_joined %>%
    dplyr::filter(mod_pos == mod_position)

  df_out <- df_joined %>%
    dplyr::mutate(
      # mod_pos (window의 modification position) 사용
      overlap_len = purrr::pmap_int(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        longest_overlap
      ),
      concatenated_RecSeq = purrr::pmap_chr(
        list(mod_containing_seq, forward_RecSeq, mod_pos),
        concat_with_overlap
      ),
      # mod_pos (window의 modification position) 사용
      mod_in_overlap = purrr::pmap_lgl(
        list(mod_containing_seq, forward_RecSeq, mod_pos, overlap_len),
        check_mod_position_in_overlap
      )
    ) %>%
    dplyr::filter(mod_in_overlap)

  result_df <<- extend_with_target_seq_and_filter(df_out, target_seq)
  return(result_df)
}
