with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    import pandas as pd
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt

    a = pd.read_csv('data/qpcr_ct.tsv', sep='\t', index_col=0)
    a['log10(copies)'] = np.log10(a['Qty'] + 1)

    fig, ax = plt.subplots(ncols=2,nrows=len(a['Virus'].unique()), figsize=(8,15), sharex='col', sharey='col')
    axix = 0
    for virus in sorted(a['Virus'].unique()):
        if axix == 0:
            ax[axix, 0].set_title('Standard curve')
            ax[axix, 1].set_title('Ct distribution for the samples')
        sns.regplot(
            data=a[(a['Type'] == 'Standard') & (a['Virus'] == virus)],
            x='log10(copies)',
            y='Ct',
            ax=ax[axix, 0]
        )
        ax[axix, 0].set_ylabel(virus)
        sns.kdeplot(
            data=a[(a['Type'] != 'Standard') & (a['Virus'] == virus)],
            x='Ct',
            ax=ax[axix, 1]
        )
        ax[axix,1].set_ylabel('')
        axix+=1
    plt.savefig('output/sup_fig2.pdf', dpi=300)
    plt.savefig('output/sup_fig2.png', dpi=300)
    plt.savefig('output/sup_fig2.tiff', dpi=300)
    plt.savefig('output/sup_fig2.eps', dpi=300)
