rule all:
  input:
    # figures
    'output/fig1.pdf',
    'output/fig2.png',
    # supplemental figs
    'output/sup_fig1.png',
    # deps
    'output/.rlibs_installed.txt'

rule supfig1:
  output:
    'output/sup_fig1.png'
  threads: 1
  script:
    'scripts/metagenomic_heatmap.py'

rule fig2:
  output:
    'output/fig2.png'
  threads: 1
  script:
    'scripts/viralloads_stripplot.py'

rule install_rlibs:
  input:
    'output/fig2.png'
  output:
    'output/.rlibs_installed.txt'
  log:
    'logs/rlibs.log'
  threads: 1
  shell:'''
  Rscript scripts/install_libs.R > {log} 2> {log}
  '''

rule fig1:
  input:
    'output/.rlibs_installed.txt'
  output:
    'output/fig1.pdf'
  threads: 1
  log:
    'logs/treeplots.log'
  shell:'''
  Rscript scripts/treeplots.R
  '''
