with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    import pandas as pd
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress 

    a = pd.read_csv('data/qpcr_ct.tsv', sep='\t', index_col=0)
    a['log10(copies)'] = np.log10(a['Qty'] + 1)

    fig, ax = plt.subplots(ncols=2,nrows=len(a['Virus'].unique()), figsize=(8,15), sharex='col', sharey='col')
    axix = 0
    for virus in sorted(a['Virus'].unique()):
        if axix == 0:
            ax[axix, 0].set_title('Standard curve')
            ax[axix, 1].set_title('Ct distribution for the samples')
        # linear regression params
        slopes = []
        intercepts = []
        rs = []
        effs = []
        tdf = a[(a['Type'] == 'Standard') & (a['Virus'] == virus)]
        for assay in tdf['Assay'].unique():
            linr = linregress(
                x = tdf[tdf['Assay'] == assay]['log10(copies)'],
                y = tdf[tdf['Assay'] == assay]['Ct']
            )
            slopes.append(linr.slope)
            intercepts.append(linr.intercept)
            rs.append(linr.rvalue ** 2)
            effs.append((-1 + 10 ** (-1/linr.slope)))
        sns.regplot(
            data=tdf,
            x='log10(copies)',
            y='Ct',
            ax=ax[axix, 0]
        )
        ax[axix, 0].text(8, 42, f"slope = {round(sum(slopes)/len(slopes), 2)}")
        ax[axix, 0].text(8, 37, f"int = {round(sum(intercepts)/len(intercepts), 2)}")
        ax[axix, 0].text(8, 32, f"R² = {round(sum(rs)/len(rs), 2)}")
        ax[axix, 0].text(8, 27, f"eff = {round(sum(effs)/len(effs), 2)}")
        ax[axix, 0].set_ylabel(virus.replace('AMFV', 'AmFV'))
        ax[axix, 0].set_xlim(0,12)
        ax[axix, 0].set_ylim(0,50)
        tdf2 = a[(a['Type'] != 'Standard') & (a['Virus'] == virus)]
        sns.kdeplot(
            data=tdf2.groupby('Sample').mean(),
            x='Ct',
            ax=ax[axix, 1]
        )
        ax[axix,1].set_ylabel('')
        axix+=1
    plt.savefig('output/sup_fig2.pdf', dpi=300)
    plt.savefig('output/sup_fig2.png', dpi=300)
    plt.savefig('output/sup_fig2.tiff', dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.savefig('output/sup_fig2.eps', dpi=300)
