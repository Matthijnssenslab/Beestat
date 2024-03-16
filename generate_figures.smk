include: 'rules/fig.smk'
include: 'rules/supp.smk'
include: 'rules/extra.smk'

figures = ['output/fig1.pdf', 'output/fig2.pdf', 'output/fig3.pdf', 'output/fig4.pdf']
supps = ['output/sup_fig1.pdf', 'output/sup_fig2.pdf', 'output/sup_fig7.pdf', 'output/sup_fig8.pdf', 'output/sup_fig9.pdf']
deps = ['data/out_rlibs_installed.txt', 'data/out_roc.tsv', 'data/out_conf.tsv', 'data/out_cvsimple_or.tsv']


rule all:
  input:
    figures, supps, deps
