
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
	# - make "safe" version of rentrez::entrez_summary to catch errors
	safe_entrez_summary <- purrr::safely(rentrez::entrez_summary)
	entrez_summary_gb <- function(id, col_select) {

		# Download data
		res <- safe_entrez_summary(db = "nucleotide", id = id)
		
		# Early exit if error with entrez
		if (!is.null(res$error)) {
			warning("No esummary records found in file, returning empty tibble")
			return(tibble(taxid = NA))
		}
		
		res$result %>%
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

# Data loading ----

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

# Rmarkdown ----

# Specify format for percent and number in Rmd
percent <- function(...) {scales::percent(...)}
number <- function(...) {scales::number(big.mark = ",", ...)}

# Monophyly and ages ----

#' Get parent node from a tip taxon
#'
#' @param tree Phylogenetic tree.
#' @param species Single tip of the tree.
#'
#' @return Number of the parent node of the species
#' 
get_parent <- function(tree, species) {
	node <- which(tree$tip.label == species)
	phangorn::Ancestors(tree, node, type = "parent")
}

#' Add "major clade" to Sanger sampling table
#'
#' "major clade" includes pteridophyte suborder or order, and combines
#' Ophioglossales + Psilotales into one group
#'
#' @param data Tibble including columns "suborder" and "order".
#'
#' @return Tibble withh column "major_clade" added
#'
add_major_clade <- function(data) {
	data %>%
		mutate(
			major_clade = coalesce(suborder, order),
			major_clade = case_when(
				major_clade %in% c("Ophioglossales", "Psilotales") ~ "Ophioglossales + Psilotales", #nolint
				TRUE ~ major_clade
			)
		)
}

#' Make tibble summarizing sampling of Sanger dataset
#'
#' @param plastome_metadata_renamed Plastome data with final (resolved) species
#' names.
#' @param sanger_tree Sanger ML tree.
#' @param ppgi_taxonomy PPGI taxonomy
#'
#' @return Tibble with columns "species",  "genus", "order", "suborder",
#' "family"  "subfamily"  "major_clade" "outgroup"
#'
make_sanger_sampling_tbl <- function(
	plastome_metadata_renamed,
	sanger_tree, ppgi_taxonomy
) {
	
	# check monophyly ----
	# Make tibble of outgroup species
	og_species <-
		plastome_metadata_renamed %>%
		select(species, outgroup) %>%
		filter(outgroup == TRUE)
	
	# Make tibble with one row per species in Sanger sampling
	tibble(species = sanger_tree$tip.label) %>%
		# Add higher-level taxonomy
		mutate(
			genus = str_split(species, "_") %>% map_chr(1)
		) %>%
		left_join(
			select(
				ppgi_taxonomy, order, suborder, family, subfamily, genus), by = "genus"
		) %>%
		# Add major_clade
		add_major_clade() %>%
		# Add outgroup status
		left_join(og_species, by = "species") %>%
		mutate(outgroup = replace_na(outgroup, FALSE)) %>%
		verify(sum(outgroup) == nrow(og_species)) %>%
		# Check for match for tips with tree
		verify(all(species %in% sanger_tree$tip.label)) %>%
		verify(all(sanger_tree$tip.label %in% .$species))
}

#' Get results of monophyly test for various taxa
#'
#' @param solution Result of assessing monophyly with assess_monophy().
#' @param taxlevels Numeric vector: taxonomic levels to extract.
#'
#' @return Tibble
#'
get_result_monophy <- function(solution, taxlevels) {
	MonoPhy::GetResultMonophyly(solution, taxlevels = taxlevels) %>%
		magrittr::extract2(1) %>%
		rownames_to_column("taxon") %>%
		as_tibble() %>%
		janitor::clean_names()
}

#' Get summary of monophyly test
#'
#' @param solution Result of assessing monophyly with assess_monophy().
#' @param taxlevels Numeric vector: taxonomic levels to extract.
#'
#' @return Tibble
get_summary_monophy <- function(solution, taxlevels) {
	mp_sum <- MonoPhy::GetSummaryMonophyly(solution, taxlevels = taxlevels)
	
	mp_sum %>%
		magrittr::extract2(1) %>%
		rownames_to_column("var") %>%
		as_tibble() %>%
		mutate(tax_level = names(mp_sum)) %>%
		janitor::clean_names() %>%
		select(tax_level, var, taxa, tips)
}

#' Assess monophyly
#'
#' Wrapper around MonoPhy::AssessMonophyly()
#'
#' @param taxon_sampling Dataframe of taxa to assess for monophyly. Must
#' include column "species"
#' @param tree Phylogenetic tree.
#' @param og_taxa Character vector; pair of taxa to define the outgroup to
#' root the tree.
#' @param tax_levels Character vector; names of columns in `taxon_sampling`
#' to check for monophyly.
#'
#' @return List; results of MonoPhy::AssessMonophyly()
#'
assess_monophy <- function(
	taxon_sampling, tree,
	og_taxa = NULL,
	tax_levels) {
	tax_levels <- c("species", tax_levels) %>% unique()
	# Root tree
	if (!is.null(og_taxa)) {
		tree <- phytools::reroot(
			tree,
			getMRCA(tree, og_taxa)
		)
	}
	# Check monophyly
	taxon_sampling %>%
		verify("species" %in% colnames(.)) %>%
		select(species, all_of(tax_levels)) %>%
		as.data.frame() %>%
		MonoPhy::AssessMonophyly(tree, .)
}


#' Calculate crown age based on two genera
#'
#' @param tip_pair Character vector of length two: the two genera to be used
#' to define the clade
#' @param sanger_sampling Tibble with one row per species in the tree, including
#' columns "genus", "family", etc.
#' @param plastid_tree_dated Dated phylogeny
#'
#' @return Character vector: age of the clade defined by the two genera
#' 
get_crown_age_from_genus_pair <- function(tip_pair, sanger_sampling, plastid_tree_dated) {
	tree_height <- max(phytools::nodeHeights(plastid_tree_dated))
	sanger_sampling %>%
		filter(genus %in% tip_pair) %>%
		group_by(genus) %>%
		slice(1) %>%
		pull(species) %>%
		ape::getMRCA(plastid_tree_dated, .) %>%
		phytools::nodeheight(plastid_tree_dated, .) %>%
		magrittr::subtract(tree_height, .) %>%
		number(accuracy = 0.1)
}


#' Calculate stem age of families
#' 
#' Will only include monophyltic and monotypic families
#'
#' @param sanger_sampling Tibble with one row per species in the tree, including
#' columns "genus", "family", etc.
#' @param plastid_tree_dated Dated phylogeny 
#' @param ppgi_taxonomy PPGI taxonomic system for ferns and lycophytes
#'
#' @return Tibble with columns "family" and "age" (stem age of the family)
get_stem_family_age <- function(
	sanger_sampling, 
	plastid_tree_dated, ppgi_taxonomy) {
	
	# Get overall tree height
	tree_height <- max(phytools::nodeHeights(plastid_tree_dated))
	
	# Get ages for families
	family_monophy <- 
		select(sanger_sampling, species, family) %>%
		assess_monophy(
			taxon_sampling = .,
			tree = plastid_tree_dated,
			tax_levels = "family"
		) %>%
		get_result_monophy(1) %>%
		rename(family = taxon) %>%
		# Filter to ferns
		left_join(unique(select(ppgi_taxonomy, family, class)), by = "family") %>%
		filter(class == "Polypodiopsida") %>%
		select(-class)
	
	# Get stem nodes of monotypic families
	family_monotypic_stem_node <-
		family_monophy %>%
		filter(monophyly == "Monotypic") %>%
		left_join(select(sanger_sampling, family, species), by = "family") %>%
		assert(is_uniq, family) %>%
		mutate(stem_node = map_dbl(species, ~get_parent(plastid_tree_dated, .))) %>%
		assert(not_na, stem_node) %>%
		select(family, stem_node)
	
	# Get stem nodes of monophyletic families
	family_monophyletic_stem_node <-
		family_monophy %>%
		filter(monophyly == "Yes") %>%
		mutate(
			mrca = parse_number(mrca),
			stem_node = map_dbl(
				mrca,
				~phangorn::Ancestors(plastid_tree_dated, ., type = "parent"))
		) %>%
		select(family, stem_node)
	
	# Combine nodes, get ages
	bind_rows(family_monotypic_stem_node, family_monophyletic_stem_node) %>%
		mutate(
			# Height is distance above root
			height = map_dbl(stem_node, ~phytools::nodeheight(
				tree = plastid_tree_dated, node = .))
		) %>%
		mutate(
			# Age is the total length of the tree minus height
			age = tree_height - height
		) %>%
		assert(is_uniq, family) %>%
		assert(not_na, everything()) %>%
		select(family, age)
}

