
# GenBank ----

#' Fetch metadata from GenBank
#'
#' @param query String to use for querying GenBank
#' @param col_select Character vector; columns of metadata to retain in output
#'
#' @return Tibble
#' 
fetch_metadata <- function(
	query = NULL,
	col_select = c("gi", "caption", "taxid", "title", "slen", "subtype", "subname")) {
	
	assertthat::assert_that(assertthat::is.string(query))
	
	assertthat::assert_that(is.character(col_select))
	
	# Do an initial search without downloading any IDs to see how many hits
	# we get.
	initial_genbank_results <- rentrez::entrez_search(
		db = "nucleotide",
		term = query,
		use_history = FALSE
	)
	
	# If no results, return empty tibble
	if (initial_genbank_results$count < 1) return(tibble(taxid = NA))

	# Download IDs with maximum set to 1 more than the total number of hits.
	genbank_results <- rentrez::entrez_search(
		db = "nucleotide",
		term = query,
		use_history = FALSE,
		retmax = initial_genbank_results$count + 1
	)
	
	# Define internal function to download genbank data into tibble
	entrez_summary_gb <- function(id, col_select) {
		# Download data
		rentrez::entrez_summary(db = "nucleotide", id = id) %>%
			# Extract selected columns from result
			purrr::map_dfr(magrittr::extract, col_select) %>%
			# Make sure taxid column is character
			mutate(taxid = as.character(taxid)) %>%
			assert(not_na, taxid)
	}
	
	# Extract list of IDs from search results
	genbank_ids <- genbank_results$ids
	
	# Fetch metadata for each ID and extract selected columns
	if (length(genbank_ids) == 1) {
		rentrez_results <- rentrez::entrez_summary(db = "nucleotide", id = genbank_ids) %>%
			magrittr::extract(col_select) %>%
			tibble::as_tibble() %>%
			mutate(taxid = as.character(taxid)) %>%
			assert(not_na, taxid)
	} else {
		# Split input vector into chunks
		n <- length(genbank_ids)
		chunk_size <- 200
		r <- rep(1:ceiling(n/chunk_size), each = chunk_size)[1:n]
		genbank_ids_list <- split(genbank_ids, r) %>% magrittr::set_names(NULL)
		# Download results for each chunk
		rentrez_results <- map_df(genbank_ids_list, ~entrez_summary_gb(., col_select = col_select))
	}
	
	return(rentrez_results)
	
}

# Get a tibble of taxids for a single year from genbank
fetch_taxa_by_year <- function(query, year, type) {
	full_query <- glue::glue('{query} AND ("{year}"[Publication Date] : "{year+1}"[Publication Date]) ')
	fetch_metadata(full_query, "taxid") %>%
		count(taxid) %>%
		mutate(type = type, year = year)
}

# Format a dataframe for searching GenBank for
# three types of fern sequences: plastid, nuclear, michondria
# from 1990 to 2021
make_gb_query <- function() {
	list(
		year = 1990:2021,
		query = c(
			"(Polypodiopsida[Organism] AND gene_in_plastid[PROP])",
			"(Polypodiopsida[Organism] gene_in_genomic[PROP])",
			"(Polypodiopsida[Organism] gene_in_mitochondrion[PROP])"
		)
	) %>%
		cross_df() %>%
		mutate(
			type = case_when(
				str_detect(query, "plastid") ~ "plastid",
				str_detect(query, "genomic") ~ "nuclear",
				str_detect(query, "mito") ~ "mitochondrial",
			)
		)
}

#' Load a dataframe of NCBI species names corresponding to taxon IDs
#' 
#' Excludes any taxon names that are not fully identified to species,
#' hybrid formulas, and environmental samples
#'
#' @param taxdump_zip_file Path to zip file with NCBI taxonomy database; must
#' contain a file called "names.dmp".
#' @param taxid_keep Character vector including the NCBI
#' taxids of the names to be extracted.
#'
#' @return Tibble with two columns, "taxid" and "species"
#' 
load_ncbi_names <- function(taxdump_zip_file, taxid_keep) {
	
	# Unzip names.dmp to a temporary directory
	temp_dir <- tempdir(check = TRUE)
	
	utils::unzip(
		taxdump_zip_file, files = "names.dmp",
		overwrite = TRUE, junkpaths = TRUE, exdir = temp_dir)
	
	# Load raw NCBI data
	ncbi_raw <-
		fs::path(temp_dir, "names.dmp") %>%
		readr::read_delim(
			delim = "\t|\t", col_names = FALSE,
			col_types = cols(.default = col_character())
		)
	
	# Delete temporary unzipped file
	fs::file_delete(fs::path(temp_dir, "names.dmp"))
	
	# Prune raw NBCI names to names in metadata
	ncbi_raw %>%
		# Select only needed columns
		transmute(
			taxid = as.character(X1),
			name = X2,
			class = X4) %>%
		# Filter to only taxid in genbank data
		filter(taxid %in% unique(taxid_keep)) %>%
		# Make sure there are no hidden fields in `class`
		verify(all(str_count(class, "\\|") == 1)) %>%
		# Drop field separators in `class`
		mutate(class = str_remove_all(class, "\\\t\\|")) %>%
		# Only keep accepted name
		filter(class == "scientific name") %>%
		# Exclude names from consideration that aren't fully identified to species, 
		# environmental samples, or hybrid formulas.
		# Hybrid names *can* be parsed:
		# - "Equisetum x ferrissii" (x before specific epithet)
		# - "x Cystocarpium roskamianum" (x before nothogenus)
		# Hybrid formulas *can't* be parsed:
		# - "Cystopteris alpina x Cystopteris fragilis" (x before another species)
		mutate(
			exclude = case_when(
				str_detect(name, " sp\\.| aff\\.| cf\\.| × [A-Z]| x [A-Z]|environmental sample") ~ TRUE,
				str_count(name, " ") < 1 ~ TRUE,
				TRUE ~ FALSE
			)
		) %>%
		filter(exclude == FALSE) %>%
		select(-exclude) %>%
		select(taxid, species = name)
}

#' Count the number of species accumulated in GenBank each year by
#' genomic compartment
#'
#' @param gb_taxa Tibble with columns "taxid", "n" (number of accessions with
#' that ID), "type" (plastid, nuclear, or mitochondrial), and "year"
#' @param ncbi_names Tibble with columns "taxid" and "species" 
#' @param year_range Range of years to calculate
#'
#' @return Tibble with total number of species and accessions in genbank by year
#' for each type of genomic compartment
#' 
count_ncbi_species_by_year <- function(gb_taxa, ncbi_names, year_range) {
	
	# Filter GenBank taxaids to only those identified to species,
	# join to species names
	gb_species <-
		gb_taxa %>%
		filter(!is.na(taxid)) %>%
		inner_join(ncbi_names, by = "taxid")
	
	# Helper function to sum total species in each dataset
	# by year
	sum_species <- function(gb_species, year_select) {
		bind_rows(
			gb_species %>% filter(year <= year_select) %>%
				group_by(type) %>%
				summarize(
					n_species = n_distinct(species),
					n_acc = sum(n)
				),
			gb_species %>% filter(year <= year_select) %>%
				summarize(
					n_species = n_distinct(species),
					n_acc = sum(n)
				) %>%
				mutate(type = "total")
		) %>%
			mutate(year = year_select)
	}
	
	# Count total number of accumulated species per year
	map_df(year_range, ~sum_species(gb_species, .))
}

# Etc ----

#' Parse divergene dates in the Testo and Sundue 2016 SI file
#' 
#' Testo WL, Sundue MA (2016) A 4000-species dataset provides new insight into 
#' the evolution of ferns. Molecular Phylogenetics and Evolution 105:200–211. 
#' https://doi.org/10.1016/j.ympev.2016.09.003
#'
#' @param testo_sundue_2016_si_path Path to Testo and Sundue 2016 SI file
#'
#' @return Tibble
#' 
parse_ts_dates <- function(testo_sundue_2016_si_path) {
	readxl::read_excel(
		testo_sundue_2016_si_path,
		# Divergence dates are in the 5th sheet
		sheet = 5,
		skip = 1
	) %>%
		janitor::clean_names() %>%
		rename(
			ts_median = median_age_2,
			ts_low = x95_percent_hpd_low_3,
			ts_high = x95_percent_hpd_high_4,
			rothfels_median = median_age_5,
			rothfels_low = x95_percent_hpd_low_6,
			rothfels_high = x95_percent_hpd_high_7,
			schuettpelz_best = best_age
		)
}
