# supplemental.
rule supfig1:
  output:
    'output/sup_fig1.pdf'
  log:
    'logs/supfig1.log'
  threads: 1
  script:
    '../scripts/supfig1.py'

rule supfig8:
  output:
    'output/sup_fig8.pdf'
  log:
    'logs/supfig8.log'
  threads: 1
  script:
    '../scripts/supfig8.py'

rule supfig9:
  input:
    'data/out_cvfull_or.tsv',
    'data/out_cvyearassoc_or.tsv',
    'data/out_cvsimple_or.tsv'
  output:
    'output/sup_fig9.pdf'
  log:
    'logs/supfig9.log'
  script:
    '../scripts/supfig9.py'

rule supfig10:
  output:
    'output/sup_fig10.pdf'
  log:
    'logs/supfig10.log'
  threads: 1
  script:
    '../scripts/supfig10.py'
