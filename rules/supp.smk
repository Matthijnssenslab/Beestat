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


rule supfig9:
  output:
    'output/sup_fig9.pdf'
  log:
    'logs/supfig9.log'
  threads: 1
  script:
    '../scripts/supfig9.py'

rule supfig10:
  input:
    'data/out_cvfull_or.tsv',
    'data/out_cvyearassoc_or.tsv',
    'data/out_cvsimple_or.tsv'
  output:
    'output/sup_fig10.pdf'
  log:
    'logs/supfig10.log'
  script:
    '../scripts/supfig10.py'

rule supfig11:
  output:
    'output/sup_fig11.pdf'
  log:
    'logs/supfig11.log'
  threads: 1
  script:
    '../scripts/supfig11.py'