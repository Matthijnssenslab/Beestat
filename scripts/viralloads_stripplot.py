import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator

meta = pd.read_csv('data/metadata.tsv', sep='\t', index_col=0)
sample_statdic = {}
sample_coldic = {}
for i,r in meta.iterrows():
    sample_statdic[i] = r['status_year']
    sample_coldic[i] = r['colony_id']

a = pd.read_csv('data/qpcr_abs.tsv', sep='\t', index_col=0)
a['sample'] = a.index
a['stat_year'] = a['sample'].map(sample_statdic)
a['colony'] = a['sample'].map(sample_coldic)
amelt = a.melt(id_vars=['sample', 'stat_year', 'colony'])
amelt
amelt.columns = ['sample','stat_year', 'colony', 'virus', 'absval']
fig, ax = plt.subplots(nrows=2, figsize=(10,8))
amelt['log_absval'] = np.log10(amelt['absval']+1)
sns.boxplot(
    data=amelt[amelt['stat_year'].isin(['healthy_12', 'wl_12'])],
    x='virus',
    y='log_absval',
    hue='stat_year',
    dodge=True,
    ax=ax[0],
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
    ax=ax[0],
    palette=['#CCBA72', '#79402E'],
    hue_order=['healthy_12', 'wl_12']
)
ax[0].legend([],[], frameon=False)
ax[0].set_title('2012')
ax[0].set_xlabel('')
ax[0].set_ylabel('log10(absolute viralload + 1)')

annotator = Annotator(
    ax=ax[0],
    pairs=[
       (("dwv", "healthy_12"), ("dwv", "wl_12")),
       (("orthomyxo", "healthy_12"), ("orthomyxo", "wl_12")),
       (("apthili", "healthy_12"), ("apthili", "wl_12")),
       (("bmlv", "healthy_12"), ("bmlv", "wl_12")),
       (("aparli", "healthy_12"), ("aparli", "wl_12")),
       (("amfv", "healthy_12"), ("amfv", "wl_12")),
       (("vdv", "healthy_12"), ("vdv", "wl_12"))
    ],
    data=amelt[amelt['stat_year'].isin(['healthy_12', 'wl_12'])],
    x='virus',
    y='log_absval',
    hue='stat_year',
    hue_order=['healthy_12', 'wl_12']
)
annotator.configure(test='Mann-Whitney', text_format='simple', loc='inside', comparisons_correction="BH", verbose=2, correction_format="replace")
annotator.apply_and_annotate()

plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))

sns.boxplot(
    data=amelt[amelt['stat_year'].isin(['healthy_13', 'wl_13'])],
    x='virus',
    y='log_absval',
    hue='stat_year',
    dodge=True,
    ax=ax[1],
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
    ax=ax[1],
    palette=['#CCBA72', '#79402E'],
    hue_order=['healthy_13', 'wl_13']
)
ax[1].legend([],[], frameon=False)
ax[1].set_title('2013')
ax[1].set_xlabel('')
ax[1].set_ylabel('log10(absolute viralload + 1)')


annotator = Annotator(
    ax=ax[1],
    pairs=[
       (("dwv", "healthy_13"), ("dwv", "wl_13")),
       (("orthomyxo", "healthy_13"), ("orthomyxo", "wl_13")),
       (("apthili", "healthy_13"), ("apthili", "wl_13")),
       (("bmlv", "healthy_13"), ("bmlv", "wl_13")),
       (("aparli", "healthy_13"), ("aparli", "wl_13")),
       (("amfv", "healthy_13"), ("amfv", "wl_13")),
       (("vdv", "healthy_13"), ("vdv", "wl_13"))
    ],
    data=amelt[amelt['stat_year'].isin(['healthy_13', 'wl_13'])],
    x='virus',
    y='log_absval',
    hue='stat_year',
    hue_order=['healthy_13', 'wl_13']
)
annotator.configure(test='Mann-Whitney', text_format='simple', loc='inside', comparisons_correction="BH", verbose=2, correction_format="replace")
annotator.apply_and_annotate()

health_patch = mpatches.Patch(color='#CCBA72', label='healthy')
weak_patch = mpatches.Patch(color='#79402E', label='weak')
ax[1].legend(title='Colony Status', bbox_to_anchor=(1.15,1.2), handles = [health_patch, weak_patch])
plt.savefig('output/fig2.png', dpi=300, bbox_inches='tight')
