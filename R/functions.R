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
