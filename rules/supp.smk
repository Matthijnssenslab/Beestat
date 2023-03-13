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
    'output/sup_fig8.pdf'
  log:
    'logs/supfig8.log'
  threads: 1
  script:
    '../scripts/supfig8.py'
