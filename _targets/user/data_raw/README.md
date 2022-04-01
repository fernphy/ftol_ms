# README

## Data files used in this analysis
### Data files included in this repo

- 2881_dating_Species_extracted.tree: Phylogenetic tree of angiosperms and outgroup taxa in Li et al (2019). The tree file was extracted from the 2881_dating_Species.tree nexus file downloaded from [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.bq091cg) using Dendroscope.
- 2881_dating_Species_extracted_og.txt: List of outgroup (non-angiosperm) taxa in 2881_dating_Species_extracted.tree.
- fern_fossils_metadata.csv: Metadata for the ferncal [fern fossils dataset v1.0.0](https://github.com/fernphy/ferncal/archive/refs/tags/v1.0.0.zip).
- pteridocat_ppgi_diff_notes.csv: Manually entered notes on differences between pteridocat and PPG I (2016) taxonomic systems for ferns.

### Data files hosted elsewhere

- 1-s2.0-S1055790316302287-mmc2.xlsx: Supplemental Data 1 from [Testo and Sundue (2016)](https://doi.org/10.1016/j.ympev.2016.09.003).
- angiosperm-time-tree-2.0-v2.0.zip: Ramírez-Barahona et al 2020 supplementary data on Zenodo https://zenodo.org/record/4721917 downloaded from https://zenodo.org/record/4721917/files/spiritu-santi/angiosperm-time-tree-2.0-v2.0.zip.
- cla12457-sup-0003-tables1-s5.docx: Supplemental data from Du et al (2021).

## Other data files

- main_library.yaml: A complete [Zotero](https://www.zotero.org/) reference library (or a symlink to it) exported with the [Better BibTeX (BBT)](https://retorque.re/zotero-better-bibtex/) plugin in YAML format. It is not needed for running the code (so is not included in the github repo), but is needed if any new references are cited in the manuscript. After citing new references, `[R/process_refs.R](R/process_refs.R)` should be run to generate the final, filtered reference list.

## References

Du X, Lu J, Zhang L, et al (2021) Simultaneous diversification of Polypodiales and angiosperms in the Mesozoic. Cladistics 37:518–539. https://doi.org/10.1111/cla.12457

Li H-T, Yi T-S, Gao L-M, et al (2019) Origin of angiosperms and the puzzle of the Jurassic gap. Nat Plants 5:461-470. https://doi.org/10/ghdjtm

Ramírez-Barahona S, Sauquet H, Magallón S (2020) The delayed and geographically heterogeneous diversification of flowering plant families. Nature Ecology & Evolution 4:1232-1238. https://doi.org/10/gg45rb

Testo WL, Sundue MA (2016) A 4000-species dataset provides new insight into the evolution of ferns. Molecular Phylogenetics and Evolution 105:200-211. https://doi.org/10.1016/j.ympev.2016.09.003
