#' Ambiguous DNA values for nucleotide matching.
#'
#' A list mapping each nucleotide (or ambiguous code) to the set of acceptable bases.
#'
#' @return A list where each element is a string of acceptable bases.
#'
#' @examples
#' ambiguous_dna_values$N  # "ACGT"
ambiguous_dna_values <- list(
  A = "A",
  C = "C",
  G = "G",
  T = "T",
  M = "AC",
  R = "AG",
  W = "AT",
  S = "CG",
  Y = "CT",
  K = "GT",
  V = "ACG",
  H = "ACT",
  D = "AGT",
  B = "CGT",
  N = "ACGT"
)

#########################################################################
#' Reverse complement helper for uppercase A/T/G/C.
#'
#' @param dna A character string representing a DNA sequence.
#'
#' @return The reverse complement of \code{dna}.
#'
#' @examples
#' rev_comp("ATGC")  # returns "GCAT"
#'
#' @export
rev_comp <- function(dna) {
  rev_dna <- paste0(rev(strsplit(dna, "")[[1]]), collapse = "")
  dna_revcomp <- chartr("ATGC", "TACG", rev_dna)
  return(dna_revcomp)
}

#########################################################################
#' Find the longest overlap between two strings with ambiguous matching.
#'
#' @description
#' For two input strings \code{a} and \code{b}, this function finds the length of 
#' the largest suffix of \code{a} that is also a prefix of \code{b}. Ambiguous nucleotides
#' are matched according to the values in \code{ambiguous_dna_values}. If either string is
#' \code{NA} or empty, the function returns 0.
#'
#' @param a A character string.
#' @param b A character string.
#'
#' @return An integer representing the length of the longest overlap.
#'
#' @examples
#' longest_overlap("GAA", "AAGT")  # returns 2 because "AA" overlaps.
#'
#' @importFrom stringr str_detect fixed
#' @export
longest_overlap <- function(a, b) {
  a <- trimws(a)
  b <- trimws(b)
  if (is.na(a) || is.na(b) || a == "" || b == "") {
    return(0)
  }
  max_len <- 0
  limit <- min(nchar(a), nchar(b))
  for (k in seq_len(limit)) {
    overlap_match <- TRUE
    for (i in seq_len(k)) {
      char_a <- substr(a, nchar(a) - k + i, nchar(a) - k + i)
      char_b <- substr(b, i, i)
      char_a <- toupper(char_a)
      char_b <- toupper(char_b)
      
      if (is.na(char_a) || is.na(char_b)) {
        overlap_match <- FALSE
        break
      }
      
      allowed_a <- ambiguous_dna_values[[char_a]]
      if (is.null(allowed_a)) allowed_a <- char_a
      allowed_b <- ambiguous_dna_values[[char_b]]
      if (is.null(allowed_b)) allowed_b <- char_b
      
      if (!(stringr::str_detect(allowed_a, stringr::fixed(char_b)) ||
            stringr::str_detect(allowed_b, stringr::fixed(char_a)))) {
        overlap_match <- FALSE
        break
      }
    }
    if (overlap_match) {
      max_len <- k
    }
  }
  return(max_len)
}

#########################################################################
#' Concatenate two strings by removing duplicated ambiguous overlap.
#'
#' @description
#' Given two strings \code{a} and \code{b}, this function uses \code{longest_overlap()} 
#' to determine the length of the overlapping region between the suffix of \code{a} and 
#' the prefix of \code{b} (with ambiguous nucleotide matching). It then concatenates 
#' \code{a} with the remainder of \code{b} after the overlap.
#'
#' @param a A character string.
#' @param b A character string.
#'
#' @return A concatenated character string.
#'
#' @examples
#' concat_with_overlap("GAA", "AAGT")  # returns "GAAGT"
#'
#' @export
concat_with_overlap <- function(a, b) {
  k <- longest_overlap(a, b)
  if (k == 0) {
    return(paste0(a, b))
  } else {
    return(paste0(a, substr(b, k+1, nchar(b))))
  }
}

#########################################################################
filter_meth_blocked <- function(methylation_sensitivity) {
  methylation_sensitivity %>%
    dplyr::filter(cleavage == "blocked")
}

filter_meth_by_mod_type <- function(methylation_sensitivity, mod_type) {
  methylation_sensitivity %>%
    dplyr::filter(modification == mod_type)
}

cross_join_and_filter_na <- function(df_filtered, df_meth_blocked) {
  tidyr::crossing(df_filtered, df_meth_blocked) %>%
    dplyr::filter(!is.na(mod_containing_seq) & !is.na(forward_RecSeq))
}

check_mod_position_in_overlap <- function(context, forward_RecSeq, mod_position, overlap_len) {
  if (is.na(overlap_len) || overlap_len == 0 || is.na(mod_position)) {
    return(FALSE)
  }
  if (mod_position <= overlap_len) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

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

check_mod_context_in_target <- function(concatenated_RecSeq, target_seq, mod_position) {
  if (is.na(concatenated_RecSeq) || is.na(target_seq) || is.na(mod_position)) {
    return(NA_character_)
  }
  
  L_rec <- nchar(concatenated_RecSeq)
  L_tar <- nchar(target_seq)
  
  if (mod_position > L_rec) {
    return(NA_character_)
  }
  
  # mod_position 기준 최소 5bp 앞뒤 확보해서 context 만들어보자
  context_window_size <- 5
  start_pos <- max(1, mod_position - (context_window_size - 1))
  end_pos <- min(L_rec, mod_position + (context_window_size - 1))
  mod_context <- substr(concatenated_RecSeq, start_pos, end_pos)
  
  # 이제 target_seq에 mod_context가 ambiguous하게 포함되는지 본다.
  L_context <- nchar(mod_context)
  
  for (tar_start in 1:(L_tar - L_context + 1)) {
    target_window <- substr(target_seq, tar_start, tar_start + L_context - 1)
    
    if (ambiguous_match(mod_context, target_window)) {
      # 매칭된다면 확장 서열 만들어서 반환
      left_part <- substr(target_seq, 1, tar_start - 1)
      right_part <- substr(target_seq, tar_start + L_context, L_tar)
      return(paste0(left_part, mod_context, right_part))
    }
  }
  
  return(NA_character_)
}

ambiguous_match <- function(seq1, seq2) {
  if (nchar(seq1) != nchar(seq2)) {
    return(FALSE)
  }
  
  for (i in seq_len(nchar(seq1))) {
    char1 <- substr(seq1, i, i)
    char2 <- substr(seq2, i, i)
    
    allowed1 <- ambiguous_dna_values[[char1]]
    allowed2 <- ambiguous_dna_values[[char2]]
    
    if (is.null(allowed1)) allowed1 <- char1
    if (is.null(allowed2)) allowed2 <- char2
    
    if (char1 == "N") allowed1 <- "ACGT"
    if (char2 == "N") allowed2 <- "ACGT"
    
    set1 <- strsplit(allowed1, "")[[1]]
    set2 <- strsplit(allowed2, "")[[1]]
    
    if (length(intersect(set1, set2)) == 0) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

check_full_extended_match <- function(extended_seq, target_seq) {
  if (nchar(extended_seq) != nchar(target_seq)) {
    return(FALSE)
  }
  
  for (i in seq_len(nchar(extended_seq))) {
    char_ext <- toupper(substr(extended_seq, i, i))
    char_tar <- toupper(substr(target_seq, i, i))
    
    allowed_ext <- ambiguous_dna_values[[char_ext]]
    allowed_tar <- ambiguous_dna_values[[char_tar]]
    
    if (is.null(allowed_ext)) allowed_ext <- char_ext
    if (is.null(allowed_tar)) allowed_tar <- char_tar
    
    if (char_ext == "N") allowed_ext <- "ACGT"
    if (char_tar == "N") allowed_tar <- "ACGT"
    
    set_ext <- strsplit(allowed_ext, "")[[1]]
    set_tar <- strsplit(allowed_tar, "")[[1]]
    
    if (length(intersect(set_ext, set_tar)) == 0) {
      cat(sprintf(
        "Mismatch at position %d: %s(%s) vs %s(%s)\n",
        i, char_ext, allowed_ext, char_tar, allowed_tar
      ))
      cat(sprintf("Extended_seq: %s\nTarget_seq:   %s\n\n", extended_seq, target_seq))
      return(FALSE)
    }
  }
  
  return(TRUE)
}

extend_with_target_seq_and_filter <- function(df_out, target_seq) {
  df_out %>%
    dplyr::mutate(
      extended_seq = purrr::pmap_chr(
        list(concatenated_RecSeq, mod_position, forward_RecSeq),
        function(seq, mod_pos, forward_seq) {
          if (is.na(seq) || is.na(mod_pos)) {
            return(NA_character_)
          }
          
          # mod_position이 concatenated_RecSeq 내에서 어디 있는지 찾는다
          mod_pos_in_concat <- nchar(forward_seq) - (nchar(seq) - mod_pos)
          
          if (mod_pos_in_concat <= 0 || mod_pos_in_concat > nchar(seq)) {
            return(NA_character_)
          }
          
          # concatenated_RecSeq와 target_seq 오버랩을 고려해서 확장한 최종 시퀀스를 만든다.
          extended <- extend_to_target_seq(seq, target_seq)
          
          if (!is.na(extended)) {
            # 확장된 것이 target_seq와 ambiguous하게 완전히 일치하는지 본다.
            if (check_full_extended_match(extended, target_seq)) {
              return(extended)
            }
          }
          
          return(NA_character_)
        }
      )
    ) %>%
    dplyr::filter(!is.na(extended_seq))
}

extend_to_target_seq <- function(concatenated_RecSeq, target_seq) {
  L_concat <- nchar(concatenated_RecSeq)
  L_target <- nchar(target_seq)
  
  # 1. concatenated_RecSeq의 끝과 target_seq의 앞부분이 어디까지 겹칠 수 있는지 찾는다.
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
  
  # 2. 겹친 부분을 기준으로 확장된 서열 만든다.
  extended_seq <- paste0(
    concatenated_RecSeq,
    substr(target_seq, max_overlap_len + 1, L_target)
  )
  
  return(extended_seq)
}


match_enzyme_sequences <- function(mod_type, target_seq) {
  df_meth_blocked <- filter_meth_blocked(methylation_sensitivity)
  df_meth_filtered <- filter_meth_by_mod_type(df_meth_blocked, mod_type)
  
  df_joined <- cross_join_and_filter_na(methylasensitive_window, df_meth_filtered)
  
  df_out <- df_joined %>%
    dplyr::mutate(
      overlap_len = purrr::map2_int(mod_containing_seq, forward_RecSeq, longest_overlap),
      concatenated_RecSeq = purrr::map2_chr(mod_containing_seq, forward_RecSeq, concat_with_overlap),
      mod_in_overlap = purrr::pmap_lgl(
        list(mod_containing_seq, forward_RecSeq, mod_position, overlap_len),
        check_mod_position_in_overlap
      )
    ) %>%
    dplyr::filter(mod_in_overlap)
  
  result_df <<- extend_with_target_seq_and_filter(df_out, target_seq)
  
  return(result_df)
}

