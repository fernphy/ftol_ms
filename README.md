# ftol_ms

Code to generate the manuscript "An open and continuously updated fern tree of life (FTOL)" Nitta et al. 2022.

## Data

A summary of all data files is given in the [data README](_targets/user/data_raw/README.md).

Data files not included in this repo that need to be manually added before running code are described below.

### Data hosted elsewhere

The following data files need to be downloaded to `_targets/user/data_raw`:

- 1-s2.0-S1055790316302287-mmc2.xlsx: Supplemental Data 1 from [Testo and Sundue (2016)](https://doi.org/10.1016/j.ympev.2016.09.003).
- angiosperm-time-tree-2.0-v2.0.zip: Supplemental data from [Ramírez-Barahona et al 2020](https://doi.org/10.1038/s41559-020-1241-3) on [Zenodo](https://zenodo.org/record/4721917), available at https://zenodo.org/record/4721917/files/spiritu-santi/angiosperm-time-tree-2.0-v2.0.zip.
- cla12457-sup-0003-tables1-s5.docx: Supplemental data from [Du et al (2021)](https://doi.org/10.1111/cla.12457).

### Data produced by other workflows

- `ftol_cache`: This folder (or a symlink to it) must be present in the project root. It contains the [targets](https://github.com/ropensci/targets) cache for the [FTOL project](https://github.com/fernphy/ftol). In the FTOL project, that folder is called `_targets`. So either the `_targets` folder from the FTOL project (after that analysis has been run to completion) should be copied here as `ftol_cache`, or `ftol_cache` should be a symlink to it.

## Docker

[A docker image](https://hub.docker.com/r/joelnitta/ftol_ms) is provided to run the code.

You can launch an instance of RStudio and access it in a web browser like so:

```
cd ftol_ms # Navigate to the folder containing ftol_ms first
docker run --rm -dt -v ${PWD}:/home/rstudio/ftol_ms -p 8787:8787 -e DISABLE_AUTH=true joelnitta/ftol_ms:latest
```

Navigate to `http://localhost:8787/` and you should have an instance of RStudio running with this project loaded. From there, run `targets::tar_make()` to run the code.

Note that if `ftol_cache` is a (relative) symlink, you will need to mount the parent directory containing its target as well as this repo.

## References

Du X, Lu J, Zhang L, et al (2021) Simultaneous diversification of Polypodiales and angiosperms in the Mesozoic. Cladistics 37:518–539. https://doi.org/10.1111/cla.12457

Ramírez-Barahona S, Sauquet H, Magallón S (2020) The delayed and geographically heterogeneous diversification of flowering plant families. Nature Ecology & Evolution 4:1232-1238. https://doi.org/10/gg45rb

Testo WL, Sundue MA (2016) A 4000-species dataset provides new insight into the evolution of ferns. Molecular Phylogenetics and Evolution 105:200-211. https://doi.org/10.1016/j.ympev.2016.09.003

## Licenses

- Code: [MIT](LICENSE)
- Data (files in `_targets/user/data_raw`): [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)
- Preprint: [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) 
