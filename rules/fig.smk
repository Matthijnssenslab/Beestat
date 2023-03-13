rule fig1:
  input:
    'data/out_rlibs_installed.txt'
  output:
    'output/fig1.pdf'
  threads: 1
  log:
    'logs/fig1.log'
  shell:'''
  Rscript scripts/fig1.R > {log} 2> {log}
  '''

rule fig2:
  input:
    'data/out_roc.tsv'
  output:
    'output/fig2.pdf'
  log:
    'logs/fig2.log'
  threads: 1
  script:
    '../scripts/fig2.py'

rule fig3:
  output:
    'output/fig3.pdf'
  log:
    'logs/fig3.log'
  threads: 1
  script:
    '../scripts/fig3.py'

rule fig4:
  output:
    'output/fig4.pdf'
  log:
    'logs/fig4.log'
  threads: 1
  script:
    '../scripts/fig4.py'
