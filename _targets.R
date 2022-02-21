library(targets)
library(tarchetypes)
library(tidyverse)
library(assertr)
library(glue)

source("R/functions.R")

# Specify path to FTOL cache
# FIXME: change to ftol_cache when ready
ftol_cache <- here::here("working/_targets")

tar_plan(
  
  # Load data ----
  # - Dates from Testo and Sundue 2016 SI
  # (includes Rothfels 2015, Schuettpelz & Pryer 2009)
  tar_files_input(
    testo_sundue_2016_si_path,
    "_targets/user/data_raw/1-s2.0-S1055790316302287-mmc2.xlsx"
  ),
  other_dates = parse_ts_dates(testo_sundue_2016_si_path),
  # - Notes on differences between Pteridocat and PPGI
  tar_files_input(
    pteridocat_ppgi_diff_notes_path,
    "_targets/user/data_raw/pteridocat_ppgi_diff_notes.csv"
  ),
  pteridocat_ppgi_diff_notes = read_csv(
    pteridocat_ppgi_diff_notes_path
  ),
  # - Fossil ferns data
  tar_files_input(
    fossil_ferns_path,
    fs::path(ftol_cache, "user/data_raw/Fossils_Ferns.csv")
  ),
  fossil_ferns_raw = read_csv(fossil_ferns_path, skip = 1) %>%
    janitor::clean_names(),
  # - Seed plant trees from Barahona et al 2020
  tar_files_input(
    spermato_trees_path,
    "_targets/user/data_raw/angiosperm-time-tree-2.0-v2.0.zip"
  ),
  spermato_trees = load_barahona_trees(spermato_trees_path),
  # - Targets from FTOL workflow
  plastome_metadata_renamed = tar_read(
    plastome_metadata_renamed,
    store = ftol_cache
  ),
  plastome_tree = tar_read(
    plastome_tree,
    store = ftol_cache
  ),
  plastid_tree_dated = tar_read(
    plastid_tree_dated,
    store = ftol_cache
  ),
  ppgi_taxonomy = tar_read(
    ppgi_taxonomy,
    store = ftol_cache
  ),
  pteridocat = tar_read(
    pteridocat,
    store = ftol_cache
  ),
  monophy_by_clade = tar_read(
    monophy_by_clade,
    store = ftol_cache
  ),
  
  # GenBank ----
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
  tar_file(taxdump_zip_file, fs::path(
    ftol_cache, "user/data_raw/taxdmp_2022-02-01.zip")),
  ncbi_names = load_ncbi_names(
    taxdump_zip_file,
    unique(gb_taxa$taxid)
  ),
  # - Count number of species and accessions in GenBank
  # by genomic compartment per year
  gb_species_by_year = count_ncbi_species_by_year(
    gb_taxa, ncbi_names, year_range = 1990:2021),
  
  # Coverage ----
  # Make tibble of accepted species in pteridocat
  accepted_species = get_accepted_species(pteridocat, ppgi_taxonomy),
  # Calculate coverage by rank and overall coverage for FTOL
  coverage_by_rank = calculate_ftol_coverage(
    sanger_sampling, accepted_species, "cov_by_rank"),
  total_coverage = calculate_ftol_coverage(
    sanger_sampling, accepted_species, "total_sampling"),
  # Calculate coverage by rank and overall coverage for backbone tree
  bb_coverage_by_rank = calculate_backbone_coverage(
    plastome_tree, sanger_sampling, accepted_species, "cov_by_rank"),
  bb_total_coverage = calculate_backbone_coverage(
    plastome_tree, sanger_sampling, accepted_species, "total_sampling"),
  
  # Analyze monophyly and ages ----
  # - Make Sanger sampling table
  sanger_sampling = make_sanger_sampling_tbl(
    plastome_metadata_renamed,
    sanger_tree = plastid_tree_dated,
    ppgi_taxonomy = ppgi_taxonomy),
  # - Crown age of various clades
  crown_ages = list(
    ferns = c("Equisetum", "Polypodium"),
    leptos = c("Osmunda", "Polypodium"),
    polypodiales = c("Lindsaea", "Polypodium"),
    eupoly_i = c("Hypodematium", "Polypodium"),
    eupoly_ii = c("Cystopteris", "Asplenium")) %>%
    map(~get_crown_age_from_genus_pair(
      ., sanger_sampling, plastid_tree_dated)),
  # - Stem age of fern families
  family_stem_ages = get_stem_family_age(
    sanger_sampling, plastid_tree_dated, ppgi_taxonomy),
  # - Read in Du 2021 dates
  tar_files_input(
    du_2021_si_path,
    "_targets/user/data_raw/cla12457-sup-0003-tables1-s5.docx"),
  du_dates_all = parse_du_dates(du_2021_si_path),
  # Summarize fern monophyly
  # - filter to ferns, add taxonomic levels
  fern_monophy_by_clade_tax_level = add_tax_levels_to_monophyly(
    monophy_by_clade, ppgi_taxonomy, sanger_sampling
  ),
  # - calculate summary table
  fern_monophy_summ_tbl = summarize_fern_monophyly(
    fern_monophy_by_clade_tax_level, ppgi_taxonomy
  ),
  # - make table of non-monophyletic genera
  fern_nonmono_gen = make_nonmonophy_gen_tab(
    fern_monophy_by_clade_tax_level, ppgi_taxonomy
  ),
  
  # Render manuscript ----
  # track files used for rendering MS
  tar_file(ref_files, list.files("ms", "references", full.names = TRUE)),
  tar_render(
    ms_doc,
    "ms/manuscript.Rmd",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  ),
  tar_render(
    ms_pdf,
    "ms/manuscript.Rmd",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  ),
  tar_render(
    si_doc,
    "ms/si.Rmd",
    knit_root_dir = "ms",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  ),
  tar_render(
    si_pdf,
    "ms/si.Rmd",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  ),
  tar_render(
    si_2_doc,
    "ms/si2.Rmd",
    knit_root_dir = "ms",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  ),
  tar_render(
    si_2_pdf,
    "ms/si2.Rmd",
    knit_root_dir = "ms",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  )
)