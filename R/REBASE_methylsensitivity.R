#' Fetch and Process REBASE Methylation Sensitivity Data
#'
#' @description
#' This function fetches methylation sensitivity data for 4mC, 5mC, and 6mA modifications from REBASE,
#' processes the data into a tidy format, and optionally saves it as an Excel file.
#'
#' @param save_to_file Logical. If TRUE, saves the output to an Excel file. Default is FALSE.
#' @param output_path Character. Path to save the Excel file if save_to_file is TRUE. Default is "REBASE_MethylSensitivity.xlsx".
#'
#' @return A data frame containing processed methylation sensitivity data with columns for enzyme name,
#'         recognition sequence, modification position, complementary modification position, cleavage information, and modification type.
#'
#' @examples
#' df <- fetch_rebase_methylation()
#' df_with_file <- fetch_rebase_methylation(save_to_file = TRUE, output_path = "methylation_data.xlsx")
#'
#' @import httr
#' @import xml2
#' @import stringr
#' @import dplyr
#' @import tidyr
#' @import xlsx
#' @export

fetch_rebase_methylation <- function(save_to_file = FALSE, output_path = "REBASE_MethylSensitivity.xlsx") {
  parse_rebase_records <- function(url) {
    response <- httr::GET(url)
    raw_html <- httr::content(response, as = "text", encoding = "UTF-8")
    doc <- xml2::read_html(raw_html)
    body_node <- xml2::xml_find_first(doc, "//body")
    html_content <- as.character(body_node)
    header_pattern <- "</font>\\s*<br>\\s*<br>"
    header_loc <- stringr::str_locate(html_content, header_pattern)
    if (!is.na(header_loc[1, "end"])) {
      data_text <- substr(html_content, header_loc[1, "end"] + 1, nchar(html_content))
    } else {
      data_text <- html_content
    }
    lines <- stringr::str_split(data_text, "<br>")[[1]]
    lines <- stringr::str_trim(lines)
    lines <- lines[lines != ""]
    enzyme_pattern <- "^[A-Za-z0-9]+I{1,2}\\b"
    records <- list()
    current_record <- NULL
    for (line in lines) {
      if (stringr::str_detect(line, enzyme_pattern)) {
        if (!is.null(current_record)) {
          records <- append(records, list(current_record))
        }
        current_record <- line
      } else {
        if (is.null(current_record)) {
          current_record <- line
        } else {
          current_record <- paste(current_record, line, sep = "\n")
        }
      }
    }
    if (!is.null(current_record)) {
      records <- append(records, list(current_record))
    }
    df <- base::data.frame(Record = unlist(records), stringsAsFactors = FALSE)
    df <- df %>%
      dplyr::mutate(
        Enzyme = stringr::str_extract(Record, enzyme_pattern),
        Info = stringr::str_trim(stringr::str_replace(Record, enzyme_pattern, ""))
      ) %>%
      dplyr::select(Enzyme, Info) %>%
      tidyr::separate_rows(Info, sep = ";\\s*") %>%
      dplyr::mutate(Info = stringr::str_trim(Info))
    return(df)
  }
  
  process_info <- function(info_string) {
    s_mod <- stringr::str_replace_all(info_string, "<font color=\\\"#00BBBB\\\"[^>]*>([^<]+)</font>", "[[MOD]]\\1[[/MOD]]")
    s_mod <- stringr::str_replace_all(s_mod, "<font color=\\\"orange\\\"[^>]*>([^<]+)</font>", "[[COMP]]\\1[[/COMP]]")
    s_clean <- stringr::str_replace_all(s_mod, "<[^>]+>", "")
    s_clean <- stringr::str_trim(s_clean)
    forward <- stringr::str_extract(s_clean, "^[^(]+")
    forward <- stringr::str_trim(forward)
    forward_seq_no_markers <- stringr::str_replace_all(forward, "\\[\\[MOD\\]\\]|\\[\\[/MOD\\]\\]|\\[\\[COMP\\]\\]|\\[\\[/COMP\\]\\]", "")
    mod_loc <- stringr::str_locate(s_clean, "\\[\\[MOD\\]\\]")
    mod_position <- if (!any(is.na(mod_loc))) nchar(stringr::str_replace_all(stringr::str_sub(s_clean, 1, mod_loc[1] - 1), "\\[\\[MOD\\]\\]|\\[\\[/MOD\\]\\]|\\[\\[COMP\\]\\]|\\[\\[/COMP\\]\\]", "")) + 1 else NA
    comp_loc <- stringr::str_locate(s_clean, "\\[\\[COMP\\]\\]")
    comp_position <- if (!any(is.na(comp_loc))) nchar(stringr::str_replace_all(stringr::str_sub(s_clean, 1, comp_loc[1] - 1), "\\[\\[MOD\\]\\]|\\[\\[/MOD\\]\\]|\\[\\[COMP\\]\\]|\\[\\[/COMP\\]\\]", "")) + 1 else NA
    return(list(forward = forward_seq_no_markers, mod_position = mod_position, comp_position = comp_position))
  }
  
  endpoints <- c("4", "5", "6")
  mod_names <- c("4mC", "5mC", "6mA")
  methylation_sensitivity_list <- list()
  
  for (i in seq_along(endpoints)) {
    mod_endpoint <- endpoints[i]
    mod_name <- mod_names[i]
    df <- parse_rebase_records(paste0("https://rebase.neb.com/cgi-bin/ms1flatlist?", mod_endpoint))
    df_processed <- df %>%
      dplyr::mutate(
        forward_RecSeq = sapply(Info, function(x) process_info(x)$forward),
        mod_position = sapply(Info, function(x) process_info(x)$mod_position),
        comp_mod_position = sapply(Info, function(x) process_info(x)$comp_position),
        cleavage = stringr::str_extract(Info, "\\(([^)]+)\\)$"),
        cleavage = stringr::str_replace_all(cleavage, "[()]", ""),
        Info = stringr::str_replace(Info, "\\s*\\([^)]*\\)$", ""),
        modification = mod_name
      )
    methylation_sensitivity_list[[mod_name]] <- df_processed
  }
  
  methylation_sensitivity <- dplyr::bind_rows(methylation_sensitivity_list)
  assign("methylation_sensitivity", methylation_sensitivity, envir = .GlobalEnv)
  if (save_to_file) {
    xlsx::write.xlsx(methylation_sensitivity, output_path)
  }
  return(methylation_sensitivity)
}
