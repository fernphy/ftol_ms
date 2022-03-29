library(targets)
library(tarchetypes)
library(tidyverse)
library(assertr)
library(glue)

source("R/functions.R")

# Specify path to FTOL cache
ftol_cache <- here::here("ftol_cache")

# Define Rmd targets outside of main workflow
tmp <- capture.output({
  ms_doc_tar <- tar_render(
    ms_doc,
    "ms/manuscript.Rmd",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  )
  ms_pdf_tar <- tar_render(
    ms_pdf,
    "ms/manuscript.Rmd",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  )
  si_doc_tar <- tar_render(
    si_doc,
    "ms/si.Rmd",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  )
  si_pdf_tar <- tar_render(
    si_pdf,
    "ms/si.Rmd",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  )
  si_doc_2_tar <- tar_render(
    si_doc_2,
    "ms/si2.Rmd",
    output_dir = "results",
    output_format = "bookdown::word_document2",
    params = list(doc_type = "doc")
  )
  si_pdf_2_tar <- tar_render(
    si_pdf_2,
    "ms/si2.Rmd",
    output_dir = "results",
    output_format = "bookdown::pdf_document2",
    params = list(doc_type = "pdf")
  )
  },
  type = "message"
)

tar_plan(
  
  # Load data ----
  # - Dates from Testo and Sundue 2016 SI
  # (includes Rothfels 2015, Schuettpelz & Pryer 2009)
  tar_file_read(
    other_dates,
    "_targets/user/data_raw/1-s2.0-S1055790316302287-mmc2.xlsx",
    parse_ts_dates(file = !!.x)
  ),
  # - Notes on differences between Pteridocat and PPGI
  tar_file_read(
    pteridocat_ppgi_diff_notes,
    "_targets/user/data_raw/pteridocat_ppgi_diff_notes.csv",
    readr::read_csv(file = !!.x)
  ),
  # - Fossil ferns metadata
  tar_file_read(
    fossil_meta,
    "_targets/user/data_raw/fern_fossils_metadata.csv",
    readr::read_csv(file = !!.x)
  ),
  # - Seed plant trees from Barahona et al 2020
  tar_file_read(
    spermato_trees,
    "_targets/user/data_raw/angiosperm-time-tree-2.0-v2.0.zip",
    load_barahona_trees(zip_path = !!.x)
  ),
  # - Targets from FTOL workflow
  tar_file_read(
    plastome_iqtree_log,
    fs::path(ftol_cache, "user/intermediates/iqtree/plastome/plastome_alignment.phy.log"), #nolint
    readr::read_lines(file = !!.x)
  ),
  # original log file from first, single Sanger run
  tar_file_read(
    sanger_iqtree_log_2000, # 2000 iterations were run total
    fs::path(ftol_cache, "user/intermediates/iqtree/sanger/sanger_alignment.phy.2022-02-23.log"), #nolint
    readr::read_lines(file = !!.x),
  ),
  tar_file_read(
    sanger_iqtree_log,
    fs::path(ftol_cache, "objects/sanger_ml_log"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    treepl_cv_results,
    fs::path(ftol_cache, "user/intermediates/treepl/con/treepl_cv_out.txt"),
    readr::read_lines(file = !!.x)
  ),
  tar_file_read(
    plastome_metadata_renamed,
    fs::path(ftol_cache, "objects/plastome_metadata_renamed"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    plastome_tree,
    fs::path(ftol_cache, "objects/plastome_tree"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    sanger_tree_dated,
    fs::path(ftol_cache, "objects/sanger_con_tree_dated"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    ppgi_taxonomy,
    fs::path(ftol_cache, "objects/ppgi_taxonomy"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    pteridocat,
    fs::path(ftol_cache, "objects/pteridocat"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    con_monophy_by_clade,
    fs::path(ftol_cache, "objects/con_monophy_by_clade"),
    readRDS(file = !!.x)
  ),
  tar_file_read(
    ts_sanger_tree_dated,
    fs::path(ftol_cache, "objects/ts_sanger_tree_dated"),
    readRDS(file = !!.x)
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
  tar_file_read(
    ncbi_names,
    fs::path(ftol_cache, "user/data_raw/taxdmp_2022-02-01.zip"),
    load_ncbi_names(taxdump_zip_file = !!.x, taxid_keep = unique(gb_taxa$taxid))
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
    sanger_tree = sanger_tree_dated,
    ppgi_taxonomy = ppgi_taxonomy),
  # - Crown age of various clades
  crown_ages = list(
    ferns = c("Equisetum", "Polypodium"),
    leptos = c("Osmunda", "Polypodium"),
    polypodiales = c("Lindsaea", "Polypodium"),
    eupoly_i = c("Hypodematium", "Polypodium"),
    eupoly_ii = c("Cystopteris", "Asplenium")) %>%
    map(~get_crown_age_from_genus_pair(
      ., sanger_sampling, sanger_tree_dated)),
  # - Stem age of fern families
  family_stem_ages = get_stem_family_age(
    sanger_sampling, sanger_tree_dated, ppgi_taxonomy),
  ts_family_stem_ages = get_stem_family_age(
    sanger_sampling, ts_sanger_tree_dated, ppgi_taxonomy),
  # - Model comparing T&S 2016 to FTOL estimated with T&S 2016 fossils
  ftol_ts_comp_mod = fit_stem_comp_mod(other_dates, ts_family_stem_ages),
  # - Read in Du 2021 dates
  tar_file_read(
    du_dates_all,
    "_targets/user/data_raw/cla12457-sup-0003-tables1-s5.docx",
    parse_du_dates(du_2021_si_path = !!.x)
  ),
  # Summarize fern monophyly
  # - filter to ferns, add taxonomic levels
  fern_con_monophy_by_clade_tax_level = add_tax_levels_to_monophyly(
    con_monophy_by_clade, ppgi_taxonomy, sanger_sampling
  ),
  # - calculate summary table
  fern_monophy_summ_tbl = summarize_fern_monophyly(
    fern_con_monophy_by_clade_tax_level, ppgi_taxonomy
  ),
  # - make table of non-monophyletic genera
  fern_nonmono_gen = make_nonmonophy_gen_tab(
    fern_con_monophy_by_clade_tax_level, ppgi_taxonomy
  ),
  
  # Render manuscript ----
  # track files used for rendering MS
  tar_file(ref_files, list.files("ms", "references", full.names = TRUE)),
  ms_doc_tar,
  # Render MS and each SI in docx (for submission) and pdf (for preprint)
  ms_pdf_tar,
  si_doc_tar,
  si_pdf_tar,
  si_doc_2_tar,
  si_pdf_2_tar
)
