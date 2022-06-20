# Install latex packages using tinytex

# This would happen automatically anyways when rendering the Rmd to pdf,
# but requires downloading packages and may not work during updates of Tex Live.
# Better to install to the docker image once and keep them there.

tinytex::tlmgr_update()

latex_packages <- c(
  "amsmath",
  "auxhook",
  "bigintcalc",
  "bitset",
  "booktabs",
  "colortbl",
  "etexcmds",
  "etoolbox",
  "euenc",
  "fancyhdr",
  "float",
  "fontspec",
  "geometry",
  "gettitlestring",
  "hycolor",
  "hyperref",
  "iftex",
  "infwarerr",
  "intcalc",
  "kvdefinekeys",
  "kvoptions",
  "kvsetkeys",
  "latex-amsmath-dev",
  "letltxmacro",
  "ltxcmds",
  "mdwtools",
  "multirow",
  "pdfescape",
  "pdftexcmds",
  "refcount",
  "rerunfilecheck",
  "stringenc",
  "tipa",
  "unicode-math",
  "uniquecounter",
  "xcolor",
  "xunicode",
  "zapfding"
)

tinytex::tlmgr_install(latex_packages)
