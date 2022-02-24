
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
#' @param file Path to Testo and Sundue 2016 SI file
#'
#' @return Tibble
#' 
parse_ts_dates <- function(file) {
  readxl::read_excel(
    file,
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

#' Parse divergene dates in the Du et al. 2021 SI file
#' 
#' Du X, Lu J, Zhang L, et al (2021) Simultaneous diversification of 
#' Polypodiales and angiosperms in the Mesozoic. Cladistics 37:518–539. 
#' https://doi.org/10/gpb22j
#'
#' @param du_2021_si_path Path to Du et al. 2021 SI file
#'
#' @return Tibble
#' 
parse_du_dates <- function(du_2021_si_path) {
  # Extract dates in original format
  du_dates_raw <-
    docxtractr::read_docx(du_2021_si_path) %>%
    docxtractr::docx_extract_all_tbls(
      guess_header = TRUE, preserve = FALSE, trim = TRUE) %>%
    # Table of interest is the third one in the SI
    magrittr::extract2(3)
  
  # Convert to longish format
  du_dates_raw %>%
    janitor::clean_names() %>%
    pivot_longer(names_to = "author", values_to = "age", -lineage) %>%
    mutate(age = str_remove_all(age, "\\(|,|\\)|ca\\.") %>%
             str_squish()) %>%
    separate(age,
             c("median", "low", "high"),
             sep = " ", fill = "right", convert = TRUE) %>%
    mutate(
      median = na_if(median, "\\") %>%
        parse_number(),
      affinities_group = str_match(
        lineage,
        "stem|crown") %>%
        magrittr::extract(,1),
      affinities = str_remove_all(lineage, "The|stem|crown|of") %>%
        str_squish()
    ) %>%
    assert(not_na, affinities_group) %>%
    select(author, lineage, affinities, affinities_group, median, low, high)
}

#' Load dated angiosperm trees from Barahona et al 2020
#' 
#' Ramírez-Barahona S, Sauquet H, Magallón S (2020) The delayed and 
#' geographically heterogeneous diversification of flowering plant families. 
#' Nature Ecology & Evolution 1–7. https://doi.org/10/gg45rb
#'
#' @param zip_path Path to zip file of supplemental data from 
#' Barahona et al 2020, downloaded from
#' https://zenodo.org/record/4721917/files/spiritu-santi/angiosperm-time-tree-2.0-v2.0.zip
#' 
#' @return List of trees named after analysis scheme
#' 
load_barahona_trees <- function(zip_path) {
  
  # Make temp dir for unzipping
  temp_dir <- fs::path(tempdir(), "angiosperm-time-tree")
  
  if (fs::dir_exists(temp_dir)) fs::dir_delete(temp_dir)
  
  fs::dir_create(temp_dir)
  
  # Files we want are in a zip file within the zip file
  # Unzip Data10_mcctrees.zip first
  
  unzip(
    zipfile = zip_path,
    files = "spiritu-santi-angiosperm-time-tree-2.0-7303f80/data/Supplementary/Data10_mcctrees.zip",
    junkpaths = TRUE,
    exdir = temp_dir
  )
  
  # Tree files of interest
  
  tree_files <- c(
    "CC_complete_MCCv_2.tre",
    "CC_conservative_MCCv_2.tre",
    "RC_complete_MCCv_2.tre",
    "RC_conservative_MCCv_2.tre",
    "UC_complete_MCCv_2.tre",
    "UC_conservative_MCCv_2.tre"
  )
  
  tree_names <- str_remove_all(tree_files, "_MCCv_2\\.tre")
  
  # Next unzip the tree files within Data10_mcctrees.zip
  
  unzip(
    zipfile = fs::path(temp_dir, "Data10_mcctrees.zip"),
    files = tree_files,
    junkpaths = TRUE,
    exdir = temp_dir
  )
  
  # Read in as list of trees, named by analysis scheme
  res <-
    tree_files %>%
    fs::path(temp_dir, .) %>%
    map(~ape::read.nexus(.)) %>%
    set_names(tree_names)
  
  if (fs::dir_exists(temp_dir)) fs::dir_delete(temp_dir)
  
  return(res)
  
}

#' Load fossil fern data
#'
#' @param file Path to fossil fern data (CSV file)
#'
#' @return Tibble
load_fossil_data <- function(file) {
  read_csv(file) %>%
    janitor::clean_names()
}

# Coverage ----

#' Make a tibble of accepted species in pteridocat taxonomic database
#'
#' @param pteridocat Taxonomic database of pteridophytes at species level
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy (genus level
#' and above)
#'
#' @return Tibble of accepted species with higher level taxonomy
#' 
get_accepted_species <- function(pteridocat, ppgi_taxonomy) {
  pteridocat %>% 
    # Only keep accepted names
    filter(str_detect(taxonomicStatus, "accepted")) %>%
    # Only keep species
    # FIXME: need to add rank for newly added names in pteridocat_maker
    # then don't need is.na()
    filter(taxonRank == "species" | is.na(taxonRank)) %>%
    select(scientificName) %>%
    mutate(scientificName = str_squish(scientificName)) %>%
    mutate(
      rgnparser::gn_parse_tidy(scientificName) %>% 
        select(species = canonicalsimple),
      species = str_replace_all(species, " ", "_")
    ) %>%
    # Map on higher level taxonomy
    mutate(
      genus = case_when(
        str_detect(scientificName, "^x |^× ") ~ 
          str_split(scientificName, " ") %>% map_chr(2),
        TRUE ~ str_split(scientificName, " ") %>% map_chr(1)
      )
    ) %>%
    left_join(
      select(ppgi_taxonomy, genus, family, 
             suborder, order, class), by = "genus"
    ) %>%
    # FIXME: update class for nothogenera in PPG
    mutate(class = case_when(
      str_detect(
        scientificName, "Schizoloma|Phlebosia|Cystocarpium") ~ "Polypodiopsida",
      TRUE ~ class
    )) %>%
    # Drop lycophytes
    assert(not_na, class) %>%
    filter(class == "Polypodiopsida") %>%
    select(-class) %>%
    # Add major_clade
    add_major_clade() %>%
    # drop suborder
    select(-suborder) 
}

#' Calculate coverage of species and other taxonomic ranks in FTOL
#'
#' @param sanger_sampling Tibble including sampling statistics in Sanger 
#' dataset (one row per species, with columns for higher-level taxonomy
#' and outgroup status)
#' @param accepted_species Tibble of accepted species in pteridocat
#' @param ret_type Type of output to return:
#' - cov_by_rank: Returns coverage by taxonmic rank
#' - total_sampling: Returns overall sampling coverage
#' 
#' @return Tibble
#'
calculate_ftol_coverage <- function(
  sanger_sampling, 
  accepted_species, ret_type = c("cov_by_rank", "total_sampling")) {
  
  # Tally accepted sampling by rank
  accepted_sampling_by_rank <- tally_accepted_sampling_by_rank(
    accepted_species)
  
  # Tally FTOL sampling by rank
  ftol_sampling_by_rank <-
    sanger_sampling %>%
    # Exclude outgroup
    filter(outgroup == FALSE) %>%
    select(-outgroup) %>%
    # drop suborder
    select(-suborder) %>%
    pivot_longer(names_to = "rank", values_to = "name", -species) %>%
    group_by(rank) %>%
    count(name) %>%
    ungroup
  
  # Calculate coverage by rank
  coverage_by_rank <- left_join(
    ftol_sampling_by_rank,
    rename(accepted_sampling_by_rank, n_accepted = n),
    by = c("rank", "name")
  ) %>%
    mutate(coverage = n/n_accepted)
  
  # Format overall sampling for printing in MS
  total_sampling <-
    # Start with FTOL sampling by genus, family, order
    ftol_sampling_by_rank %>% 
    count(rank) %>%
    # Add n species
    bind_rows(
      sanger_sampling %>%
        filter(outgroup == FALSE) %>%
        select(-outgroup) %>%
        summarize(n = n(), rank = "species")
    ) %>%
    # Join to accepted number of species, genus, family, order
    left_join(
      count(accepted_sampling_by_rank, rank, name = "n_accepted") %>%
        bind_rows(
          summarize(accepted_species, n_accepted = n(), rank = "species")
        ),
      by = "rank"
    ) %>%
    # Calculate coverage
    mutate(
      coverage_p = percent(n / n_accepted, accuracy = 0.1),
      coverage = glue("{number(n, accuracy = 1)}/{number(n_accepted, accuracy = 1)}")) # nolint
  
  switch(
    ret_type,
    "cov_by_rank" = coverage_by_rank,
    "total_sampling" = total_sampling,
    stop("Must choose 'cov_by_rank'  or 'total_sampling' for 'ret_type'")
  )
  
}

#' Tally the number of accepted species in pteridocat by rank
#'
#' @param accepted_species Tibble of accepted species in pteridocat
#'
#' @return Tibble
#' 
tally_accepted_sampling_by_rank <- function(accepted_species) {
  accepted_species %>%
    select(-scientificName) %>%
    pivot_longer(names_to = "rank", values_to = "name", -species) %>%
    group_by(rank) %>%
    count(name) %>%
    ungroup %>%
    # Note that nothogenera have no family or order
    filter(!is.na(name))
}

#' Calculate coverage of species and other taxonomic ranks in plastome
#' (backbone) tree
#'
#' @param plastome_tree Plastome (backbone) tree
#' @param sanger_sampling Tibble including sampling statistics in Sanger 
#' dataset (one row per species, with columns for higher-level taxonomy
#' and outgroup status)
#' @param accepted_species Tibble of accepted species in pteridocat
#' @param ret_type Type of output to return:
#' - cov_by_rank: Returns coverage by taxonmic rank
#' - total_sampling: Returns overall sampling coverage
#' 
#' @return Tibble
#' 
calculate_backbone_coverage <- function(
  plastome_tree, sanger_sampling,
  accepted_species, ret_type = c("cov_by_rank", "total_sampling")) {
  
  # Tally accepted sampling by rank
  accepted_sampling_by_rank <- tally_accepted_sampling_by_rank(accepted_species)
  
  # Tally backbone sampling by rank
  bb_sampling <-
    tibble(species = plastome_tree$tip.label) %>%
    # Add higher taxonomy
    left_join(sanger_sampling, by = "species") %>%
    # Drop OG
    filter(outgroup == FALSE) %>%
    select(-outgroup) %>%# drop suborder
    # drop suborder
    select(-suborder)
  
  # Calculate coverage by rank
  bb_sampling_by_rank <-
    bb_sampling %>%
    pivot_longer(names_to = "rank", values_to = "name", -species) %>%
    group_by(rank) %>%
    count(name) %>%
    ungroup %>%
    left_join(
      rename(accepted_sampling_by_rank, n_accepted = n),
      by = c("rank", "name")
    ) %>%
    mutate(coverage = n/n_accepted)
  
  # Format sampling for printing in MS
  bb_total_sampling <-
    # Start with backbone sampling by genus, family, order
    bb_sampling_by_rank %>% 
    count(rank) %>%
    # Add n species
    bind_rows(
      summarize(bb_sampling, n = n(), rank = "species")
    ) %>%
    # Join to accepted number of species, genus, family, order
    left_join(
      count(accepted_sampling_by_rank, rank, name = "n_accepted") %>%
        bind_rows(
          summarize(accepted_species, n_accepted = n(), rank = "species")
        ),
      by = "rank"
    ) %>%
    # Calculate coverage
    mutate(
      coverage_p = percent(n / n_accepted),
      coverage = glue("{number(n, accuracy = 1)}/{number(n_accepted, accuracy = 1)}"))
  
  switch(
    ret_type,
    "cov_by_rank" = bb_sampling_by_rank,
    "total_sampling" = bb_total_sampling,
    stop("Must choose 'cov_by_rank'  or 'total_sampling' for 'ret_type'")
  )
  
}

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
#' @param sanger_tree_dated Dated phylogeny
#'
#' @return Character vector: age of the clade defined by the two genera
#' 
get_crown_age_from_genus_pair <- function(tip_pair, sanger_sampling, sanger_tree_dated) {
  tree_height <- max(phytools::nodeHeights(sanger_tree_dated))
  sanger_sampling %>%
    filter(genus %in% tip_pair) %>%
    group_by(genus) %>%
    slice(1) %>%
    pull(species) %>%
    ape::getMRCA(sanger_tree_dated, .) %>%
    phytools::nodeheight(sanger_tree_dated, .) %>%
    magrittr::subtract(tree_height, .) %>%
    number(accuracy = 0.1)
}


#' Calculate stem age of families
#' 
#' Will only include monophyltic and monotypic families
#'
#' @param sanger_sampling Tibble with one row per species in the tree, including
#' columns "genus", "family", etc.
#' @param sanger_tree_dated Dated phylogeny 
#' @param ppgi_taxonomy PPGI taxonomic system for ferns and lycophytes
#'
#' @return Tibble with columns "family" and "age" (stem age of the family)
get_stem_family_age <- function(
  sanger_sampling, 
  sanger_tree_dated, ppgi_taxonomy) {
  
  # Get overall tree height
  tree_height <- max(phytools::nodeHeights(sanger_tree_dated))
  
  # Get ages for families
  family_monophy <- 
    select(sanger_sampling, species, family) %>%
    assess_monophy(
      taxon_sampling = .,
      tree = sanger_tree_dated,
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
    mutate(stem_node = map_dbl(species, ~get_parent(sanger_tree_dated, .))) %>%
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
        ~phangorn::Ancestors(sanger_tree_dated, ., type = "parent"))
    ) %>%
    select(family, stem_node)
  
  # Combine nodes, get ages
  bind_rows(family_monotypic_stem_node, family_monophyletic_stem_node) %>%
    mutate(
      # Height is distance above root
      height = map_dbl(stem_node, ~phytools::nodeheight(
        tree = sanger_tree_dated, node = .))
    ) %>%
    mutate(
      # Age is the total length of the tree minus height
      age = tree_height - height
    ) %>%
    assert(is_uniq, family) %>%
    assert(not_na, everything()) %>%
    select(family, age)
}

#' Filter monophyly table to ferns, add taxonomic level
#'
#' @param monophy_by_clade Monophyly status for various taxonomic levels,
#' including ferns and outgroups
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy (genus level
#' and above)
#' @param sanger_sampling Tibble with one row per species in the tree, including
#' columns "genus", "family", etc. 
#'
#' @return Tibble
#' 
add_tax_levels_to_monophyly <- function(
  monophy_by_clade, ppgi_taxonomy,
  sanger_sampling) {
  
  ppgi_tax_levels <-
    ppgi_taxonomy %>%
    select(class:genus) %>%
    pivot_longer(
      names_to = "tax_level",
      values_to = "taxon",
      -class
    ) %>%
    filter(!is.na(taxon)) %>%
    unique()
  
  # Add taxonomic level to monophy by clade table
  monophy_by_clade_tax_level <-
    monophy_by_clade %>%
    assert(not_na, monophyly) %>%
    assert(in_set("Monotypic", "Yes", "No"), monophyly) %>%
    left_join(ppgi_tax_levels, by = "taxon") %>%
    # Drop lycophytes
    filter(class != "Lycopodiopsida") %>%
    # Exclude Equisetum subgenera, outgroups
    filter(!str_detect(taxon, "subgen")) %>%
    assert(not_na, class) %>%
    verify(all(class == "Polypodiopsida")) %>%
    select(-class) %>%
    left_join(
      unique(select(sanger_sampling, taxon = genus, outgroup)),
      by = "taxon") %>%
    mutate(outgroup = replace_na(outgroup, FALSE)) %>%
    filter(outgroup != TRUE) %>%
    assert(not_na, tax_level)
}

#' Summarize (non-)monophyletic status of ferns in FTOL
#'
#' @param monophy_by_clade_tax_level Tibble with monophyly status by clade,
#' including taxonomic level
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy (genus level
#' and above)
#' @param check Logical; check that Polypodioideae is the only taxonomic
#' level above genus that is non-monophyletic
#'
#' @return Tibble
#' 
summarize_fern_monophyly <- function(
  monophy_by_clade_tax_level, ppgi_taxonomy,
  check = TRUE) {
  
  if (isTRUE(check)) {
    monophy_by_clade_tax_level %>%
      filter(tax_level != "genus", monophyly == "No") %>%
      verify(
        taxon == "Polypodioideae",
        success_fun = success_logical,
        error_fun = err_msg(
          "Polypodioideae is not the only non-monophyletic taxon above genus")
      )
  }
  
  # Get order of taxonomic levels
  tax_levels_order <-
    ppgi_taxonomy %>%
    select(order:genus) %>%
    colnames()
  
  # Make summary table
  monophy_by_clade_tax_level %>%
    group_by(tax_level) %>%
    count(monophyly) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(percent = n/total) %>%
    ungroup() %>%
    mutate(
      tax_level = factor(tax_level, levels = tax_levels_order),
      monophyly = factor(monophyly, levels = c("Monotypic", "Yes", "No"))
    ) %>%
    arrange(tax_level, monophyly)
}

#' Make table of non-monophyletic fern genera in FTOL
#'
#' @param fern_monophy_by_clade_tax_level Tibble with monophyletic status
#' of taxa in FTOL including taxonomic level
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy (genus level
#' and above)
#'
#' @return Tibble: nonmonophyletic genera in FTOL
#' 
make_nonmonophy_gen_tab <- function(
  fern_monophy_by_clade_tax_level, ppgi_taxonomy) {
  
  fern_monophy_by_clade_tax_level %>%
    assert(not_na, tax_level, outgroup) %>%
    filter(tax_level == "genus", monophyly == "No", outgroup == FALSE) %>%
    rename(genus = taxon) %>%
    select(-c(monophyly, mrca, tax_level, outgroup)) %>%
    left_join(select(ppgi_taxonomy, genus, family, subfamily), by = "genus") %>%
    select(family, subfamily, everything())
}

# Formatting ----

# Abbreviations
ie <- "*i.e.*"
eg <- "*e.g.*"
ca <- "*ca.*"

### formatters start
# Specify custom percentage format
percent <- function(...) {scales::percent(...)}

# Specify custom number format
number <- function(...) {scales::number(big.mark = ",", ...)}
### formatters end

#' Convert a data.frame of counts to percentages, with automatic formatting
#'
#' Format numbers and percentage to look nice in MS
#'
#' @param ... Arguments passed to janitor::tabyl() 
#'
#' @result Dataframe
tabyl_fmt <- function(...) {
  tabyl(...) %>%
    mutate(
      percent = percent(percent, accuracy = 0.1),
      n = number(n, accuracy = 1)
    )
}

#' Print a single category of fossil data
#'
#' @param fossil_data Fossil data in long form, filtered to a single
#' fossil
#' @param cat_select Category of data to print
#'
#' @return Character vector

print_fossil_cat <- function(fossil_data, cat_select) {
  data <-
    fossil_data %>% filter(category == cat_select)
  
  c(
    md_heading(cat_select, 3),
    md_bullet(data$text),
    md_blank()
  )
}

#' Print data for a single fossil
#'
#' @param fossil_data Fossil data in long form
#' @param fossil_select Selected fossil to print
#'
#' @return Character vector
#'
#' @examples
#' print_fossil(fossils_long, "Acrostichum intertrappeum") %>%
#' 	glue::as_glue()
print_fossil <- function(fossil_data, fossil_select) {
  # Filter data to only selected fossil
  data <-
    fossil_data %>% filter(fossil_taxon == fossil_select)
  
  # Print name of fossil as header in italics
  header <-
    data %>% 
    filter(var == "fossil_taxon_header") %>%
    pull(value) %>%
    str_replace_all("_", " ") %>%
    md_italic() %>%
    md_heading(1)
  
  categories <- data %>%
    filter(!is.na(category)) %>%
    pull(category) %>%
    unique()
  
  body <- map(categories, ~print_fossil_cat(data, .)) %>%
    unlist()
  
  c(header, body, md_rule())
}

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
    stringr::str_remove_all("^[^@]*") %>%
    stringr::str_remove_all('\\[|\\]|\\)|\\(|\\.$|,|\\{|\\}|\\\\|\\"') %>% 
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
#' @param silent Logical; should a warning be issued for missing references?
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
s_figure_full <- captioner::captioner(prefix = "Figure S", auto_space = FALSE)
s_table_full <- captioner::captioner(prefix = "Table S", auto_space = FALSE)

# - Make a short function that prints only the object type and number, e.g., "Fig. 1"
figure <- pryr::partial(figure_full, display = "cite", caption = "blank")
table <- pryr::partial(table_full, display = "cite", caption = "blank")
s_figure <- pryr::partial(s_figure_full, display = "cite", caption = "blank")
s_table <- pryr::partial(s_table_full, display = "cite", caption = "blank")

# - Make a short function that prints only the number (e.g., "1")
figure_num <- function(name) {figure(name) %>% str_remove("Figure ")}
table_num <- function(name) {table(name) %>% str_remove("Table ")}
s_figure_num <- function(name) {s_figure(name) %>% str_remove("Figure ")}
s_table_num <- function(name) {s_table(name) %>% str_remove("Table ")}

# Pagebreaks ----

#' Pagebreak for MS Word (doc) only
#' 
#' Meant for use within Rmd when rendering.
#' Breaks a page depending if the output format is MS Word
#' 
#' @param rmd_params Parameters set in YAML header or rmarkdown::render(). 
#' Must include `doc_type` (either 'doc' for MS Word or 'pdf' for PDF output)
#' 
pagebreak_doc <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "doc", return("\\newpage"), return(""))}

#' Pagebreak for PDF only
#' 
#' Meant for use within Rmd when rendering. Requires parameter variable
#' `doc_type` to be defined (either 'doc' for MS Word or 'pdf' for PDF output)
#'
#' Breaks a page depending if the output format is PDF
pagebreak_pdf <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "pdf", return("\\newpage"), return(""))}

# Phylogeny ----

#' Rescale tree height
#' 
#' from:
#' http://blog.phytools.org/2012/02/quicker-way-to-rescale-total-length-of.html
#'
#' @param tree Input phylogeny (list of class "phylo")
#' @param scale Numeric: total height to scale to
#'
#' @return Phylogeny (list of class "phylo")
#' 
rescale_tree <- function(tree, scale){
  tree$edge.length <-
    tree$edge.length/max(phytools::nodeHeights(tree)[,2])*scale
  return(tree)
}

#' Get node height corresponding to the most recent common
#' ancestor of two or more taxa in a phylogeny
#'
#' @param phy List of class "phylo"
#' @param tips Character vector; tips of the phylogeny
#'
#' @return Height of the MRCA node of the tips
#' 
node_height_from_tips <- function(phy, tips) {
  ape::getMRCA(phy, tips) %>%
    phytools::nodeheight(phy, .)
}

# Etc -----

#' Function to get support value in BS%
#'
#' @param phy an object of class "phylo".
#' @param tips a vector of mode numeric or character specifying the tips
#'
#' @return Bootstrap support as % of that node
get_bs <- function(phy, tips) {
  node_select <- ape::getMRCA(phy, tips)
  as_tibble(phy) %>%
    filter(node == node_select) %>%
    pull(label) %>%
    parse_number() %>%
    magrittr::multiply_by(0.01) %>%
    percent()
}

# Helper function to add plot_group:
# a "rank" that may be at order or suborder level for plotting
add_plot_group <- function(data) {
  data %>%
    mutate(
      plot_group = coalesce(suborder, order),
      plot_group = case_when(
        plot_group %in% c("Ophioglossales", "Psilotales") ~ "Ophioglossales + Psilotales",
        TRUE ~ plot_group
      )
    )
}

#' Get tips of a phylogenetic tree in their plotted order
#'
#' After re-rooting a tree, the order of tips when the tree
#' is plotted no longer match the order of $tip.label. Use
#' this function to get tips in the order they are plotted.
#' @param tree List of class "phylo"
#' @return Character vector
get_tips_in_ape_plot_order <- function(tree) {
  assertthat::assert_that(inherits(tree, "phylo"))
  # First filter out internal nodes
  # from the the second column of the edge matrix
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  # Use this vector to extract the tips in the right order
  tree$tip.label[ordered_tips]
}

#' Calculate frequency of missing bases in a DNA sequence matrix
#'
#' @param seqs DNA sequences in a matrix
#' @param missing_base_codes Codes to count as missing bases
#' @param start Start column of matrix to include in calculation
#' @param end End column of matrix to include in calculation
#' @param freq a logical specifying whether to return the proportions 
#' (the default) or the absolute frequencies (counts).
#'
#' @return Number: the frequency of missing bases in the matrix
#' 
base_freq_missing <- function(
  seqs, missing_base_codes = c("n", "-", "?"), 
  start = 1, end = ncol(seqs),
  freq = FALSE) {
  assertthat::assert_that(is.matrix(seqs))
  assertthat::assert_that(end > start)
  base_freq <- ape::base.freq(seqs[, start:end], all = TRUE, freq = freq)
  sum(base_freq[missing_base_codes])
}

# Generate custom error message for assertr
# https://github.com/ropensci/assertr/issues/70
err_msg <- function(msg) stop(msg, call. = FALSE)

#' Check order of cited figures and tables in an Rmd file
#'
#' @param type Type of citation to check:
#' - "all": everything
#' - "fig": all figures
#' - "s_fig": SI figures
#' - "tab": all tables
#' - "s_tab": SI tables
#' @param rmd_file Path to the Rmd file
#' @return Character vector: the functions citing the tables and references,
#' in the order they appear in the Rmd file.
check_ft_order <- function(type = "all", rmd_file = "ms/manuscript.Rmd") {
  
  # Run inline code when purling
  options(knitr.purl.inline = TRUE)
  
  # Make a temporary file to write out just inline R code from the SI Rmd
  temp_file <- tempfile()
  
  # Generate R script including inline R code from the ms Rmd
  suppressMessages(knitr::purl(rmd_file, output = temp_file))
  
  # Trim this R script down to only functions that define figure and table captions
  refs <- read_lines(temp_file) %>%
    magrittr::extract(
      str_detect(., "^figure\\(|^table\\(|^s_figure\\(|^s_table\\(|^figure_num\\(|^table_num\\(|^s_figure_num\\(|^s_table_num\\(")
    ) %>%
    unique
  
  # Delete the temporary file
  fs::file_delete(temp_file)
  
  switch(
    type,
    "all" = refs,
    "fig" = refs[str_detect(refs, "fig")],
    "s_fig" = refs[str_detect(refs, "s_fig")],
    "tab" = refs[str_detect(refs, "tab")],
    "s_tab" = refs[str_detect(refs, "s_tab")],
    stop("Must choose 'all', 'fig', 's_fig', 'tab', or 's_tab' for type")
  )
  
}

# Count the number of non-missing bases in a DNA sequence
#' @param seq List of class DNAbin of length one.
count_non_missing <- function(seq) {
  bases <- ape::base.freq(seq, all = TRUE, freq = TRUE)
  sum(bases[!names(bases) %in% c("n", "N", "-", "?")])
}

#' Extract the best-scoring model from an IQTREE log
#'
#' @param iqtree_log Raw IQTREE log file (ending in .log) read into R
#'
#' @return Best-scoring model
#' 
extract_iqtree_mod <- function(iqtree_log) {
  iqtree_log[
    # Formatted for IQTREE2
    str_detect(iqtree_log, "Best-fit model: ")] %>%
    str_match("([^ ]+) chosen") %>%
    magrittr::extract(,2)
}

#' Extract number of iterations from an IQTREE log
#'
#' @param iqtree_log Raw IQTREE log file (ending in .log) read into R
#'
#' @return Number of iterations used by IQTREE search
#' 
extract_iqtree_iter <- function(iqtree_log) {
  iqtree_log[
    str_detect(iqtree_log, "TREE SEARCH COMPLETED AFTER")] %>%
    str_match("TREE SEARCH COMPLETED AFTER ([0-9]+) ") %>%
    magrittr::extract(,2) %>%
    parse_number() %>%
    number(accuracy = 1)
}

#' Extract number of bootstrap correlation coefficient from an IQTREE log
#'
#' @param iqtree_log Raw IQTREE log file (ending in .log) read into R
#'
#' @return Number of iterations used by IQTREE search
#' 
extract_iqtree_corr <- function(iqtree_log) {
  iqtree_log[
    str_detect(iqtree_log, "Bootstrap correlation coefficient of split occurrence frequencies: ")] %>%
    str_match("Bootstrap correlation coefficient of split occurrence frequencies: ([^$]+)$") %>%
    magrittr::extract(,2)
}
