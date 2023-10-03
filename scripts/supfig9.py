with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    import pandas as pd
    from scipy.stats import spearmanr
    import numpy as np
    from itertools import combinations
    import seaborn as sns

    s2p = pd.read_table('data/sample_to_pool.txt', header=None)
    s2p.columns = ['sample', 'pool']
    s2p.index = s2p['sample']
    del s2p['sample']

    a = pd.read_table('data/metagenomic_coverage.tsv', index_col=0)
    a = a[[
        'apis mellifera filamentous virus',
        'apis orthomyxovirus 1 - PB1',
        'apthili virus',
        'apparli virus',
        'bee-macula like virus',
        'deformed wing virus/varroa-destructor virus 1'
    ]]
    a.columns = ['amfv', 'aov', 'apthili', 'aparli', 'bmlv', 'dwv/vdv']


    b = pd.read_table('data/qpcr_abs.tsv', index_col=0)
    b['dwv/vdv'] = b['dwv'] + b['vdv']
    del b['dwv']
    del b['vdv']

    b = pd.merge(s2p, b,left_index=True, right_index=True)
    b = b.groupby('pool').max()

    a = a.sort_index()
    b = b.sort_index()
    assert list(a.index) == list(b.index)

    correlations = []
    for comb in list(combinations(list(a.columns), 2)):
        c ,p = spearmanr(
            a[comb[0]],
            b[comb[1]]
        )
        c2, p2 = spearmanr(
            a[comb[1]],
            b[comb[0]]
        )
        correlations.append(
            [comb[0], comb[1], c, p]
        )
        correlations.append(
            [comb[1], comb[0], c2, p2]
        )
    for vir in list(a.columns):
        c,p = spearmanr(
            a[vir],
            b[vir]
        )
        correlations.append(
            [vir, vir, c, p]
        )
    cdf = pd.DataFrame(correlations)
    cdf.columns = ['virus1', 'virus2', 'correlation', 'pvalue']
    cordf = cdf[['virus1', 'virus2', 'correlation']].pivot(index='virus1', columns='virus2')
    cordf.columns = ['amfv', 'aov', 'apparli', 'apthili' ,'bmlv', 'dwv/vdv']
    cordf.index = ['amfv', 'aov', 'apparli', 'apthili' ,'bmlv', 'dwv/vdv']
    pvaldf = cdf[['virus1', 'virus2', 'pvalue']].pivot(index='virus1', columns='virus2')
    g = sns.heatmap(cordf, annot=pvaldf, cmap='vlag', center=0)
    g.set_title("spearman correlation")
    g.set_ylabel('metagenomics')
    g.set_xlabel('qpcr')
    g.figure.savefig('output/sup_fig9.pdf', dpi=300)
