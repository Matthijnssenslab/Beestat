with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    from scipy.stats import mannwhitneyu
    from statannotations.Annotator import Annotator
    from collections import Counter

    meta = pd.read_csv('data/metadata.tsv', sep='\t', index_col=0)
    sample_statdic = {}
    sample_coldic = {}
    for i,r in meta.iterrows():
        sample_statdic[i] = r['status_year']
        sample_coldic[i] = r['colony_id']


    a = pd.read_csv('data/qpcr_abs.tsv', sep='\t', index_col=0)
    a.columns = ['DWV', 'AOV', 'apthili', 'BMLV', 'apparli', 'AMFV', 'VDV']
    a['sample'] = a.index
    a['stat_year'] = a['sample'].map(sample_statdic)
    a['colony'] = a['sample'].map(sample_coldic)
    amelt = a.melt(id_vars=['sample', 'stat_year', 'colony'])
    amelt
    amelt.columns = ['sample','stat_year', 'colony', 'virus', 'absval']
    amelt['log_absval'] = np.log10(amelt['absval']+1)

    # Calculate positives.
    percpos = []
    for v in list(amelt.virus.unique()):
        for s in list(amelt.stat_year.unique()):
            percpos.append(
                [
                    '{}'.format(v),
                    len(a[(a[v] > 1000) & (a['stat_year'] == s)])/len(a[a['stat_year'] == s]) * 100,
                    '{}'.format(s)
                ]
            )
    posdf = pd.DataFrame(percpos)
    posdf.columns = ['virus', '% positive samples', 'stat_year']

    # Number of samples with #viruses
    counts_shared = {}
    for s in list(amelt.stat_year.unique()):
        tdf = a[a['stat_year'] == s]
        clis = []
        for i, r in (tdf[['DWV','AOV', 'apthili','BMLV','apparli','AMFV', 'VDV']] > 1000).iterrows():
            clis.append(sum(r))
        counts_shared[s] = dict(Counter(clis))
    cdf = pd.DataFrame(counts_shared)
    cdf = cdf.fillna(0) 
    cdf = (cdf / cdf.sum(axis=0)) * 100
    cdf['count'] = cdf.index
    cdf = cdf.melt(id_vars='count')
    cdf.columns = ['# viruses detected', 'stat_year', '% of samples']

    # Start plot
    fig = plt.figure(constrained_layout=False, figsize=(14,12))
    gs = fig.add_gridspec(4,3)
    ax_l12 = fig.add_subplot(gs[0, 0:2])
    ax_l13 = fig.add_subplot(gs[1, 0:2])
    ax_p12 = fig.add_subplot(gs[0, 2:3])
    ax_p13 = fig.add_subplot(gs[1, 2:3])
    ax_ps12 = fig.add_subplot(gs[2,0])
    ax_ps13 = fig.add_subplot(gs[3,0])
    ax_log = fig.add_subplot(gs[2:4,1:3])

    ##############################################################
    # 2012
    sns.boxplot(
        data=amelt[amelt['stat_year'].isin(['healthy_12', 'wl_12'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        dodge=True,
        ax=ax_l12,
        boxprops=dict(alpha=.3),
        showfliers = False,
        palette=['#CCBA72', '#79402E'],
        hue_order=['healthy_12', 'wl_12']
    )
    sns.stripplot(
        data=amelt[amelt['stat_year'].isin(['healthy_12', 'wl_12'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        dodge=True,
        ax=ax_l12,
        palette=['#CCBA72', '#79402E'],
        hue_order=['healthy_12', 'wl_12']
    )
    ax_l12.legend([],[], frameon=False)
    ax_l12.set_xlabel('')
    ax_l12.set_xticklabels('')
    ax_l12.set_ylabel('log10(absolute viralload + 1)')


    annotator = Annotator(
        ax=ax_l12,
        pairs=[
           (("DWV", "healthy_12"), ("DWV", "wl_12")),
           (("AOV", "healthy_12"), ("AOV", "wl_12")),
           (("apthili", "healthy_12"), ("apthili", "wl_12")),
           (("BMLV", "healthy_12"), ("BMLV", "wl_12")),
           (("apparli", "healthy_12"), ("apparli", "wl_12")),
           (("AMFV", "healthy_12"), ("AMFV", "wl_12")),
           (("VDV", "healthy_12"), ("VDV", "wl_12"))
        ],
        data=amelt[amelt['stat_year'].isin(['healthy_12', 'wl_12'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        hue_order=['healthy_12', 'wl_12']
    )
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', comparisons_correction="BH", verbose=2, correction_format="replace")
    annotator.apply_and_annotate()

    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))

    sns.barplot(
        data=posdf[posdf['stat_year'].isin(['healthy_12', 'wl_12'])],
        x='virus',
        y='% positive samples',
        hue='stat_year',
        ax=ax_p12,
        palette=['#CCBA72', '#79402E']
    )
    ax_p12.legend([],[], frameon=False)
    ax_p12.set_xlabel('')
    ax_p12.set_xticklabels('')

    sns.barplot(
        data=cdf[cdf['stat_year'].isin(['healthy_12', 'wl_12'])],
        x='# viruses detected',
        y='% of samples',
        hue='stat_year',
        ax=ax_ps12,
        palette=['#CCBA72', '#79402E']
    )
    ax_ps12.legend([],[], frameon=False)
    ax_l12.annotate("2012", xy=(0, 0.5), xytext=(-ax_l12.yaxis.labelpad - 5, 0),
                    xycoords=ax_l12.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')
    ax_ps12.set_xticklabels('')
    ax_ps12.annotate("2012", xy=(0, 0.5), xytext=(-ax_ps12.yaxis.labelpad - 5, 0),
                    xycoords=ax_ps12.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')
    ##############################################################
    # 2013
    sns.boxplot(
        data=amelt[amelt['stat_year'].isin(['healthy_13', 'wl_13'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        dodge=True,
        ax=ax_l13,
        boxprops=dict(alpha=.3),
        showfliers = False,
        palette=['#CCBA72', '#79402E'],
        hue_order=['healthy_13', 'wl_13']
    )
    sns.stripplot(
        data=amelt[amelt['stat_year'].isin(['healthy_13', 'wl_13'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        dodge=True,
        ax=ax_l13,
        palette=['#CCBA72', '#79402E'],
        hue_order=['healthy_13', 'wl_13']
    )
    ax_l13.legend([],[], frameon=False)
    ax_l13.set_xlabel('')
    ax_l13.set_ylabel('log10(absolute viralload + 1)')
    annotator = Annotator(
        ax=ax_l13,
        pairs=[
           (("DWV", "healthy_13"), ("DWV", "wl_13")),
           (("AOV", "healthy_13"), ("AOV", "wl_13")),
           (("apthili", "healthy_13"), ("apthili", "wl_13")),
           (("BMLV", "healthy_13"), ("BMLV", "wl_13")),
           (("apparli", "healthy_13"), ("apparli", "wl_13")),
           (("AMFV", "healthy_13"), ("AMFV", "wl_13")),
           (("VDV", "healthy_13"), ("VDV", "wl_13"))
        ],
        data=amelt[amelt['stat_year'].isin(['healthy_13', 'wl_13'])],
        x='virus',
        y='log_absval',
        hue='stat_year',
        hue_order=['healthy_13', 'wl_13']
    )
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', comparisons_correction="BH", verbose=2, correction_format="replace")
    annotator.apply_and_annotate()

    sns.barplot(
        data=posdf[posdf['stat_year'].isin(['healthy_13', 'wl_13'])],
        x='virus',
        y='% positive samples',
        hue='stat_year',
        ax=ax_p13,
        palette=['#CCBA72', '#79402E']
    )
    ax_p13.set_xlabel('')
    ax_p13.legend([],[], frameon=False)

    sns.barplot(
        data=cdf[cdf['stat_year'].isin(['healthy_13', 'wl_13'])],
        x='# viruses detected',
        y='% of samples',
        hue='stat_year',
        ax=ax_ps13,
        palette=['#CCBA72', '#79402E']
    )
    ax_ps13.legend([],[], frameon=False)
    ax_l13.annotate("2013", xy=(0, 0.5), xytext=(-ax_l13.yaxis.labelpad - 5, 0),
                    xycoords=ax_l13.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')
    ax_ps13.annotate("2013", xy=(0, 0.5), xytext=(-ax_ps13.yaxis.labelpad - 5, 0),
                    xycoords=ax_ps13.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')
    # Legend
    health_patch = mpatches.Patch(color='#CCBA72', label='healthy')
    weak_patch = mpatches.Patch(color='#79402E', label='weak')
    ax_l12.legend(title='Colony Status', bbox_to_anchor=(-0.075,0.1), handles = [health_patch, weak_patch])


    # logit
    logroc = pd.read_csv('data/out_roc.tsv', sep='\t')
    sns.lineplot(
        data=logroc,
        x='fullfpr',
        y='fulltpr',
        ax=ax_log,
        color='#899DA4'
    )
    sns.lineplot(
        data=logroc,
        x='assocfpr',
        y='assoctpr',
        ax=ax_log,
        color='#C93312'
    )
    sns.lineplot(
        data=logroc,
        x='simplefpr',
        y='simpletpr',
        ax=ax_log,
        color='#DC863B'
    )
    sns.lineplot(
        x=[0,1],
        y=[0,1],
        color='#7f7f7f',
        alpha=0.2,
        ax=ax_log
    )
    ax_log.set_xlabel('False positive rate')
    ax_log.set_ylabel('True positive rate')

    full_patch = mpatches.Patch(color='#899DA4', label='full model + 1° int.')
    assoc_patch = mpatches.Patch(color='#C93312', label='assoc. viruses + 1° int.')
    simple_patch = mpatches.Patch(color='#DC863B', label='full model')
    ax_log.legend(title='formula', bbox_to_anchor=(0.95,0.35), handles = [full_patch, assoc_patch, simple_patch])


    ax_l12.text(-0.02, 1.05, 'A', transform=ax_l12.transAxes, 
                size=20, weight='bold')
    ax_p12.text(-0.02, 1.05, 'B', transform=ax_p12.transAxes, 
                size=20, weight='bold')
    ax_ps12.text(-0.02, 1.02, 'C', transform=ax_ps12.transAxes, 
                size=20, weight='bold')
    ax_log.text(-0.02, 1.01, 'D', transform=ax_log.transAxes, 
                size=20, weight='bold')
    plt.savefig('output/fig2.pdf', dpi=300)
