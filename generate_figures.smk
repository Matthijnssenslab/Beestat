import os
include: 'rules/fig.smk'
include: 'rules/supp.smk'
include: 'rules/extra.smk'
figures = ['output/fig1.png', 'output/fig2.png', 'output/fig3.png', 'output/fig4.png']
supps = ['output/sup_fig1.png', 'output/sup_fig8.png']
deps = ['data/out_rlibs_installed.txt', 'data/out_roc.tsv', 'data/out_conf.tsv']


rule all:
  input:
    figures, supps, deps
