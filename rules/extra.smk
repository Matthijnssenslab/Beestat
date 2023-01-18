# data/installs.
rule install_rlibs:
  output:
    'data/out_rlibs_installed.txt'
  log:
    'logs/rlibs.log'
  threads: 1
  shell:'''
  Rscript scripts/install_libs.R > {log} 2> {log}
  '''

rule logit_data:
  input:
    'data/out_rlibs_installed.txt'
  output:
    'data/out_roc.tsv',
    'data/out_conf.tsv'
  threads: 1
  log:
    'logs/logit_data.log'
  shell:'''
  Rscript scripts/logit.R > {log} 2> {log}
  '''
