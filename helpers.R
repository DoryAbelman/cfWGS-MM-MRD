# Collection of helper functions used across scripts

clean_sample_id <- function(x) {
  # remove spaces and convert to uppercase
  toupper(gsub("\\s+", "", x))
}
