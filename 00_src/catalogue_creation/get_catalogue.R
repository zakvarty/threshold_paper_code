

#' Write remote Groningen earthquake catalogue to a local json file
#'
#' @param file_path path to directory in which to write now_cat.json
#'
#' @return silently returns downloaded JSON
#' @export
#'
#' @examples write_groningen_json(file_path = "~/Desktop/")
write_groningen_json <- function(file_path = NULL){
  JSON <- readLines("http://cdn.knmi.nl/knmi/map/page/seismologie/all_induced.json",)

  now_string <- Sys.time()
  now_string <- stringr::str_replace_all(now_string ,':','-')
  now_string <- stringr::str_replace(now_string, " ","_")
  file_name  <- paste0(now_string,'_cat.json')

  JSON <- readLines("http://cdn.knmi.nl/knmi/map/page/seismologie/all_induced.json")
  writeLines(text = JSON, con = paste0(file_path,file_name))
  invisible(JSON)
}

#' Read remote Groningen earthquake catalogue into memory
#'
#' @param ... extra arguments passed to jsonlite::read_json
#'
#' @return list containing event data.frame
#' @export
#'
#' @examples cat <- download_groningen_json()
download_groningen_json <- function(...){
  JSON <- jsonlite::read_json("http://cdn.knmi.nl/knmi/map/page/seismologie/all_induced.json", ...)
  invisible(JSON)
}

