# Formatting ----

# Specify percentage format
percent <- function(...) {scales::percent(...)}

# Specify number format
number <- function(...) {scales::number(big.mark = ",", ...)}

# Abbreviations
ie <- "*i.e.*"
eg <- "*e.g.*"
ca <- "*ca.*"

# References ----

#' Extract citations (formatted like `[@key]`) from an Rmd file
#'
#' @param rmd_file Character vector; path to Rmd file
#'
#' @return Data frame with one column 'key'
#' 
extract_citations <- function(rmd_file) {
	read_lines(rmd_file) %>%
		stringr::str_split(" |;") %>% 
		unlist %>% 
		magrittr::extract(., stringr::str_detect(., "@")) %>% 
		stringr::str_remove_all("\\[|\\]|\\)|\\(|\\.$|,|\\{|\\}|\\\\") %>% 
		magrittr::extract(., stringr::str_detect(., "^@|^-@")) %>% 
		stringr::str_remove_all("^@|^-@") %>% 
		unique %>% 
		sort %>%
		tibble(key = .)
}

#' Filter a list of references in YAML format to those occurring in an Rmd file
#'
#' @param rmd_file Character vector; Path to Rmd file(s)
#' @param yaml_in String or list; if string, the path to the YAML file to filter.
#' If list, should be result of reading in a YAML file with yaml::read_yaml()
#' @param yaml_out Path to write filtered YAML reference file
#'
#' @return NULL; externally, the filtered YAML will be written to `yaml_out`
#' 
filter_refs_yaml <- function(rmd_file, yaml_in = "ms/main_library.yaml", yaml_out = "ms/references.yaml", silent = FALSE) {
	
	# Parse RMD file and extract citation keys
	citations <- purrr::map_df(rmd_file, extract_citations)
	
	# Read in YAML including all references exported from Zotero
	if (inherits(yaml_in, "character")) {
		ref_yaml <- yaml::read_yaml(yaml_in)
	} else if (inherits(yaml_in, "list")) {
		ref_yaml <- yaml_in
	} else {
		stop("`yaml_in` must be a path to a YAML file (string) or a list read in with yaml::read_yaml()")
	}
	
	# Extract all citation keys from full YAML
	cite_keys_all <- purrr::map_chr(ref_yaml$references, "id") %>%
		tibble::tibble(
			key = .,
			order = 1:length(.)
		)
	
	# Check that all keys in the yaml are present in the input YAML
	missing <- citations %>% dplyr::anti_join(cite_keys_all, by = "key")
	
	if (nrow(missing) > 0 && silent == FALSE)
		warning(glue::glue("The following ref keys are present in the Rmd but missing from the input YAML: {missing$key}"))
	
	cite_keys_filtered <- citations %>% dplyr::inner_join(cite_keys_all, by = "key")
	
	# Filter YAML to only those citation keys in the RMD
	ref_yaml_filtered <- list(references = ref_yaml$references[cite_keys_filtered$order])
	
	# Write out the YAML file
	yaml::write_yaml(ref_yaml_filtered, file = yaml_out)
}

# Figure output ----

#' Generate a path to save a results file
#' 
#' Only works for figures or tables cited in text, and outputs
#' files to the "results" folder
#'
#' @param result_num Number of result, e.g. "Fig. 2"
#' @param extension Extension to use for file
#'
#' @return String.
#' 
result_file <- function (result_num, extension) {
	
	fs::path(
		here::here("results"),
		result_num %>%
			str_remove_all("\\.") %>% 
			str_replace_all(" ", "_")
	) %>%
		fs::path_ext_set(extension)
}

# Captions ----

# Follow Frontiers in Plant Sciences style
# Captions should be preceded by the appropriate label,
# for example "Figure 1." 
# Figure panels are referred to by bold capital letters in brackets: 
# (A), (B), (C), (D), etc.

# - First define the "full" version, which would include a caption
# (except I never use the caption in the function, and instead replace with 'blank')
figure_full <- captioner::captioner(prefix = "Figure ", auto_space = FALSE)
table_full <- captioner::captioner(prefix = "Table ", auto_space = FALSE)
s_figure_full <- captioner::captioner(prefix = "Figure S ", auto_space = FALSE)
s_table_full <- captioner::captioner(prefix = "Table S ", auto_space = FALSE)

# - Make a short function that prints only the object type and number, e.g., "Fig. 1"
figure <- pryr::partial(figure_full, display = "cite", caption = "blank")
table <- pryr::partial(table_full, display = "cite", caption = "blank")
s_figure <- pryr::partial(s_figure_full, display = "cite", caption = "blank")
s_table <- pryr::partial(s_table_full, display = "cite", caption = "blank")

# - Make a short function that prints only the number (e.g., "1")
figure_num <- function(name) {figure(name) %>% str_remove("Figure ")}
table_num <- function(name) {table(name) %>% str_remove("Table ")}
s_figure_num <- function(name) {s_figure(name) %>% str_remove("Figure S ")}
s_table_num <- function(name) {s_table(name) %>% str_remove("Table S ")}

# Pagebreaks ----

#' Pagebreak for PDF only
#' 
#' Meant for use within Rmd when rendering. Requires parameter variable
#' `doc_type` to be defined (either 'doc' for MS Word or 'pdf' for PDF output)
#'
#' Breaks a page depending if the output format is PDF
pagebreak_pdf <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "pdf", return("\\clearpage"), return(""))}
