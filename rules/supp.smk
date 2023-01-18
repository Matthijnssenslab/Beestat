# supplemental.
rule supfig1:
  output:
    'output/sup_fig1.png'
  log:
    'logs/supfig1.log'
  threads: 1
  script:
    '../scripts/supfig1.py'

