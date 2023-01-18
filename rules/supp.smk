# supplemental.
rule supfig1:
  output:
    'output/sup_fig1.png'
  log:
    'logs/supfig1.log'
  threads: 1
  script:
    '../scripts/supfig1.py'

rule supfig2:
  output:
    'output/sup_fig2.png'
  log:
    'logs/supfig2.log'
  threads: 1
  script:
    '../scripts/supfig2.py'
