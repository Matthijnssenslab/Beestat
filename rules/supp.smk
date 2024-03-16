# supplemental.
rule supfig1:
  output:
    'output/sup_fig1.pdf'
  log:
    'logs/supfig1.log'
  threads: 1
  script:
    '../scripts/supfig1.py'

rule supfig2:
  output:
    'output/sup_fig2.pdf'
  log:
    'logs/supfig2.log'
  threads: 1
  script:
    '../scripts/supfig2.py'


rule supfig7:
  output:
    'output/sup_fig7.pdf'
  log:
    'logs/supfig7.log'
  threads: 1
  script:
    '../scripts/supfig7.py'

rule supfig8:
  input:
    'data/out_cvfull_or.tsv',
    'data/out_cvyearassoc_or.tsv',
    'data/out_cvsimple_or.tsv'
  output:
    'output/sup_fig8.pdf'
  log:
    'logs/supfig8.log'
  script:
    '../scripts/supfig8.py'

rule supfig9:
  output:
    'output/sup_fig9.pdf'
  log:
    'logs/supfig9.log'
  threads: 1
  script:
    '../scripts/supfig9.py'