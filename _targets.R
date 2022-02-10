library(targets)
library(tarchetypes)
library(tidyverse)
library(assertr)

source("R/functions.R")

tar_plan(
	# Analyze number of fern accessions and species in GenBank by type of DNA 
	# (plastid, mitochondrial, or nuclear)
	#
	# - Make dataframe of query terms by year
	gb_query = make_gb_query(),
	# - Download taxids one year at a time
	tar_target(
		gb_taxa,
		fetch_taxa_by_year(
			query = gb_query$query,
			year = gb_query$year,
			type = gb_query$type),
		pattern = map(gb_query)
		),
	# - Load NCBI taxonomy
	tar_file(taxdump_zip_file, "working/taxdmp_2022-02-01.zip"),
	ncbi_names = load_ncbi_names(
		taxdump_zip_file,
		unique(gb_taxa$taxid)
	),
	# - Count number of species and accessions in GenBank
	# by genomic compartment per year
	gb_species_by_year = count_ncbi_species_by_year(
		gb_taxa, ncbi_names, year_range = 1990:2021),
	# Read in estimated clade ages from Testo and Sundue 2016 SI
	tar_files_input(
		testo_sundue_2016_si_path,
		"_targets/user/data_raw/1-s2.0-S1055790316302287-mmc2.xlsx"
	),
	other_dates = parse_ts_dates(testo_sundue_2016_si_path),
	# Render manuscript ----
	# track files used for rendering MS
	tar_file(ref_files, list.files("ms", "references", full.names = TRUE)),
	tar_render(
		ms,
		"ms/manuscript.Rmd",
		output_dir = "results"
	)
)