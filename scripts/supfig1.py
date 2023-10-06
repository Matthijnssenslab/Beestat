with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from PIL import Image
    c = pd.read_csv('data/metagenomic_coverage.tsv', sep='\t', index_col=0)
    simpleix = []
    fac = []
    twc = 0
    thc = 0
    for s in list(c.index):
        if '2013' in s:
            thc += 1
            simpleix.append('sample_{}_2013'.format(thc))
            fac.append('2012')
        else:
            twc += 1
            simpleix.append('sample_{}_2012'.format(twc))
            fac.append('2013')

    g = sns.clustermap(
        c.T,
        row_cluster=False,
        cmap='Blues',
        vmin=0,
        col_colors = list(pd.Series(fac).map({'2012': '#D9D0D3', '2013': '#8D8680'})),
        xticklabels=False
    )
    g.savefig('output/sup_fig1.pdf', dpi=300)
    g.savefig('output/sup_fig1.png', dpi=300)
    g.savefig('output/sup_fig1.tiff', dpi=300)
