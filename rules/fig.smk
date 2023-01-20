rule fig1:
  input:
    'data/out_rlibs_installed.txt'
  output:
    'output/fig1.png'
  threads: 1
  log:
    'logs/fig1.log'
  shell:'''
  #Rscript scripts/fig1.R > {log} 2> {log}
  Rscript scripts/fig1.R
  '''

rule fig2:
  input:
    'data/out_roc.tsv'
  output:
    'output/fig2.png'
  log:
    'logs/fig2.log'
  threads: 1
  script:
    '../scripts/fig2.py'

rule fig3:
  output:
    'output/fig3.png'
  log:
    'logs/fig3.log'
  threads: 1
  script:
    '../scripts/fig3.py'

