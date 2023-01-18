options(install.packages.compile.from.source = "always")
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install('ggtree', update=FALSE, ask=FALSE)
cat(file='data/out_rlibs_installed.txt', '')
