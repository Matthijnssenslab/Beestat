with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    import numpy as pd
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from PIL import Image

    fig, ax = plt.subplots(nrows=3, sharex=True, figsize=(8,12))

    # full+1

    cdic = {'NS': '#7f7f7f', 'S':'#899DA4'}
    c = pd.read_table('data/out_cvfull_or.tsv', sep='\t')
    ixupd = []
    for _i in c.index:
        #rename vdv into dwv-b, dwv into dwv-a
        ixupd.append(_i.replace('dwv', 'dwv-a').replace('vdv', 'dwv-b'))
    c.index = ixupd
    c['sig'] = 'NS'
    c.loc[c['pval'] < 0.05, 'sig'] = 'S'
    c['coef'] = c.index
    c = c.sort_values('OR', ascending=False)
    sns.scatterplot(
        data=c,
        y='coef',
        x='OR',
        ax=ax[0],
        hue='sig',
        palette=['#7f7f7f']
    )
    ax[0].legend([],[], frameon=False)
    ax[0].axvline(1, color='red', linestyle='--')
    for i,r in c.iterrows():
        ax[0].plot([r['OR_lower_95'], r['OR_upper_95']], [i,i], c=cdic[r['sig']])
    ax[0].set_ylabel('full model+ 1° int')

    # assoc+1
    cdic = {'NS': '#7f7f7f', 'S':'#C93312'}
    b = pd.read_table('data/out_cvyearassoc_or.tsv', sep='\t')
    ixupd = []
    for _i in b.index:
        #rename vdv into dwv-b, dwv into dwv-a
        ixupd.append(_i.replace('dwv', 'dwv-a').replace('vdv', 'dwv-b'))
    b.index = ixupd
    b['sig'] = 'NS'
    b.loc[b['pval'] < 0.05, 'sig'] = 'S'
    b['coef'] = b.index
    b = b.sort_values('OR', ascending=False)
    sns.scatterplot(
        data=b,
        y='coef',
        x='OR',
        ax=ax[1],
        hue='sig',
        palette=['#C93312', '#7f7f7f']
    )
    ax[1].legend([],[], frameon=False)
    ax[1].axvline(1, color='red', linestyle='--')
    for i,r in b.iterrows():
        ax[1].plot([r['OR_lower_95'], r['OR_upper_95']], [i,i], c=cdic[r['sig']])
    ax[1].set_ylabel('assoc.viruses + 1° model')

    # 'full model'
    cdic = {'NS': '#7f7f7f', 'S':'#DC863B'}
    a = pd.read_table('data/out_cvsimple_or.tsv', sep='\t')
    ixupd = []
    for _i in a.index:
        #rename vdv into dwv-b, dwv into dwv-a
        ixupd.append(_i.replace('dwv', 'dwv-a').replace('vdv', 'dwv-b'))
    a.index = ixupd
    a['sig'] = 'NS'
    a.loc[a['pval'] < 0.05, 'sig'] = 'S'
    a['coef'] = a.index
    a = a.sort_values('OR', ascending=False)
    sns.scatterplot(
        data=a,
        y='coef',
        x='OR',
        ax=ax[2],
        hue='sig',
        palette=['#DC863B', '#7f7f7f']
    )
    ax[2].legend([],[], frameon=False)
    ax[2].axvline(1, color='red', linestyle='--')
    for i,r in a.iterrows():
        ax[2].plot([r['OR_lower_95'], r['OR_upper_95']], [i,i], c=cdic[r['sig']])
    ax[2].set_ylabel('Simple model')
    plt.tight_layout()
    plt.savefig('output/sup_fig8.pdf', dpi=300)
    plt.savefig('output/sup_fig8.png', dpi=300)
    plt.savefig('output/sup_fig8.tiff', dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.savefig('output/sup_fig8.eps', dpi=300)
