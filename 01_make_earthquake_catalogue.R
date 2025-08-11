# 01_make_earthquake_catalogue.R -----------------------

DOWNLOAD_NEW_JSON = FALSE
OUTPUT_PATH = "00_data/eathquake_catalogues/"

library(threshold)
library(stringr)

if (DOWNLOAD_NEW_JSON) { write_groningen_json(file_path = OUTPUT_PATH) }

# Get paths for latest json and csv files
json_path <- threshold:::latest_json_file(OUTPUT_PATH)
csv_path  <- stringr::str_replace(string = json_path, ".json", ".csv")

gf_cat <- threshold::read_JSON_cat(
  JSON_path = json_path,
  field_outline_with_buffer = threshold::gfo_buffer_1000)

if (DOWNLOAD_NEW_JSON) {
  readr::write_csv(gf_cat, file = csv_path)
}

