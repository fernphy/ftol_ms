# Filter a list of references in YAML format to those cited in an RMD file
library(tidyverse)
library(rmdref)

# Filter YAML references from Zotero library.
# main_library.yaml has been exported from Zotero
# like this: file -> "export library" -> "Better CSL YAML"
filter_refs_yaml(
  rmd_file = "ms/manuscript.Rmd",
  yaml_in = "_targets/user/data_raw/main_library.yaml",
  yaml_out = "ms/references.yaml", silent = TRUE)

# Same for SI
filter_refs_yaml(
  rmd_file = c(
    "ms/si.rmd", # SI
    "_targets/user/data_raw/pteridocat_ppgi_diff_notes.csv" # refs in SI tables
  ),
  yaml_in = "_targets/user/data_raw/main_library.yaml",
  yaml_out = "ms/si_references.yaml", silent = TRUE)

# Write out 'aux' file with citation keys to make collection in Zotero
# so reference data are easier to work with.
# The collection can be made with the betterbibtex for Zotero plugin:
# https://retorque.re/zotero-better-bibtex/citing/aux-scanner/
extract_citations("ms/manuscript.Rmd") %>%
  # aux file just consists of lines formatted like
  # \citation{KEY}
  # where KEY is the citation key
  mutate(latex = glue::glue("\\citation{[key]}", .open = "[", .close = "]")) %>%
  pull(latex) %>%
  write_lines("working/ftol.aux")

extract_citations(
   c(
     "ms/si.rmd", # SI
     "_targets/user/data_raw/pteridocat_ppgi_diff_notes.csv" # refs in SI tables
   )
  ) %>%
  filter(str_detect(key, "figureHmakeatother", negate = TRUE)) %>%
  filter(key != "") %>%
  mutate(latex = glue::glue("\\citation{[key]}", .open = "[", .close = "]")) %>%
  pull(latex) %>%
  write_lines("temp/ftol_si1.aux")
