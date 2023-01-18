with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import geopandas as gpd
    import geoplot as gplt
    import seaborn as sns
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from matplotlib_scalebar.scalebar import ScaleBar
    from shapely.geometry.point import Point
    import pysal as ps
    import matplotlib.patches as mpatches
    from pysal.lib.weights import Rook
    from esda import Moran
    import math
    ###################################################### Defs
    # Define a jitter to move the samples around their commune centroid a little bit.
    def jitter(values,j):
        return values + np.random.normal(j,0.01,values.shape)

    # Define dist at 4.25, 51.1 (~= roughly middle of belgium, for a km scale on the map).
    points = gpd.GeoSeries([Point(3.75, 51.1), Point(4.75, 51.1)], crs=4326)  # Geographic WGS 84 - degrees
    points = points.to_crs(32619)
    distance_ms = points[0].distance(points[1])

    ###################################################### getData 

    # Get Flanders
    be = gpd.read_file('data/gadm41_BEL_4.json')
    be['NAME_3'] = be['NAME_3'].str.lower()
    be = be[be['NAME_1'] != 'Bruxelles']
    be = be[be['NAME_1'] != 'Wallonie']

    #Metadata
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
    virdat = pd.merge(a, meta, left_index=True, right_index=True)


    bee = be.explode(index_parts=True)
    vir_bel = pd.merge(virdat, bee, left_on='commune', right_on='NAME_3')
    w_rook = Rook.from_dataframe(vir_bel)
    mi_res = []
    for v in ['dwv', 'aov', 'apthili', 'bmlv', 'aparli', 'amfv', 'vdv']:
        mi_rook = Moran(np.log10(vir_bel[v]+1), w_rook, permutations=10000)
        mi_res.append([v, mi_rook.I, mi_rook.p_sim])
    mirdf = pd.DataFrame(mi_res)
    mirdf.columns = ['virus', 'moran-I', 'p-value']
    mirdf.to_csv('data/out_moranI.tsv', sep='\t')

    ###################################################### Plotter


    ###################################################### 2012
    #fig, ax = plt.subplots(nrows=5,ncols=2, figsize=(10,15))
    fig = plt.figure(constrained_layout=False, figsize=(10,8))
    gs = fig.add_gridspec(21,2)
    ax_12 = fig.add_subplot(gs[0:7, 0])
    ax_13 = fig.add_subplot(gs[7:14, 0])
    ax_mor = fig.add_subplot(gs[14:, 0])
    ax_kde_dwv = fig.add_subplot(gs[0:3,1])
    ax_kde_aov = fig.add_subplot(gs[3:6,1])
    ax_kde_apthili = fig.add_subplot(gs[6:9,1])
    ax_kde_bmlv = fig.add_subplot(gs[9:12,1])
    ax_kde_aparli = fig.add_subplot(gs[12:15,1])
    ax_kde_amfv = fig.add_subplot(gs[15:18,1])
    ax_kde_vdv = fig.add_subplot(gs[18:,1])

    be.plot(
        edgecolor='white',
        color='gainsboro',
        ax=ax_12
    )
    sns.scatterplot(
        y=jitter(virdat[virdat['status_year'].isin(['healthy_12', 'wl_12'])]['lon'], 0.01),
        x=jitter(virdat[virdat['status_year'].isin(['healthy_12', 'wl_12'])]['lat'], 0.01),
        hue=virdat[virdat['status_year'].isin(['healthy_12', 'wl_12'])]['status_year'],
        palette=['#CCBA72', '#79402E'],
        ax=ax_12,
        alpha=0.7
    )
    ax_12.legend([],[], frameon=False)
    ax_12.axis('off')
    ax_12.title.set_text('Samples 2012')
    ax_12.margins(0)

    # Legend
    health_patch = mpatches.Patch(color='#CCBA72', label='healthy')
    weak_patch = mpatches.Patch(color='#79402E', label='weak')
    ax_12.legend(title='Colony Status', bbox_to_anchor=(0,0.1), handles = [health_patch, weak_patch])
    ax_12.text(-0.05, 1.05, 'A', transform=ax_12.transAxes, 
                size=20, weight='bold')
    ###################################################### 2013
    be.plot(
        edgecolor='white',
        color='gainsboro',
        ax=ax_13
    )
    sns.scatterplot(
        y=jitter(virdat[virdat['status_year'].isin(['healthy_13', 'wl_13'])]['lon'], 0.01),
        x=jitter(virdat[virdat['status_year'].isin(['healthy_13', 'wl_13'])]['lat'], 0.01),
        hue=virdat[virdat['status_year'].isin(['healthy_13', 'wl_13'])]['status_year'],
        palette=['#CCBA72', '#79402E'],
        ax=ax_13,
        alpha=0.7
    )
    ax_13.legend([],[], frameon=False)
    ax_13.axis('off')
    ax_13.add_artist(ScaleBar(distance_ms))
    ax_13.title.set_text('Samples 2013')

    # Morans
    sns.barplot(
        data=mirdf,
        x='moran-I',
        y='virus',
        color='#9986A5',
        ax=ax_mor
    )

    ax_mor.set(xlim=(0, 1), xlabel='moran-I')
    ixcount = 0
    for i,r in mirdf.iterrows():
        morI = r['moran-I']
        if morI > 0:
            placeX = morI + 0.005
        elif morI < 0:
            placeX = 0 + 0.005
        if r['p-value'] > 0.05:
            sigstat ='ns'
        else:
            sigstat = 'p-val<1e-{}'.format(abs(math.floor(np.log10(r['p-value']))))
        ax_mor.text(placeX, ixcount, sigstat)
        ixcount +=1 

    ax_mor.text(-0.05, 1.05, 'B', transform=ax_mor.transAxes, 
                size=20, weight='bold')

    ########################################################################### kde
    ################ dwv
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_dwv
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['dwv']+1),
        color='#A42820',
        fill=True,
        alpha=0.5,
        ax=ax_kde_dwv
    )
    ax_kde_dwv.text(-0.05, 1.05, 'C', transform=ax_kde_dwv.transAxes,
                    size=20, weight='bold')
    ax_kde_dwv.add_artist(ScaleBar(distance_ms, location="upper center", border_pad=-2))
    ################ aov
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_aov
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['aov']+1),
        color='#5F5647',
        fill=True,
        alpha=0.5,
        ax=ax_kde_aov
    )

    ################ apthili
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_apthili
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['apthili']+1),
        color='#9B110E',
        fill=True,
        alpha=0.5,
        ax=ax_kde_apthili
    )


    ################ bmlv
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_bmlv
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['bmlv']+1),
        color='#3F5151',
        fill=True,
        alpha=0.5,
        ax=ax_kde_bmlv
    )

    ################ aparli
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_aparli
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['aparli']+1),
        color='#4E2A1E',
        fill=True,
        alpha=0.5,
        ax=ax_kde_aparli
    )

    ################ amfv
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_amfv
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['amfv']+1),
        color='#550307',
        fill=True,
        alpha=0.5,
        ax=ax_kde_amfv
    )

    ################ vdv
    be.plot(
        edgecolor='gainsboro',
        color='white',
        ax=ax_kde_vdv
    )
    sns.kdeplot(
        data=virdat,
        x='lat',
        y='lon',
        weights=np.log10(virdat['vdv']+1),
        color='#0C1707',
        fill=True,
        alpha=0.5,
        ax=ax_kde_vdv
    )

    # set kde axs properly.

    for axname,x in zip(
        ['dwv', 'aov', 'apthili', 'bmlv', 'aparli', 'amfv', 'vdv'],
        [ax_kde_dwv,ax_kde_aov,ax_kde_apthili, ax_kde_bmlv,ax_kde_aparli, ax_kde_amfv, ax_kde_vdv]
    ):
        x.set(xlim=(2, 6.5), ylim=(50.5, 51.5))
        x.legend([],[], frameon=False)
        #x.margins(0)
        # make xaxis invisibel
        x.xaxis.set_visible(False)
        # remove ticks and labels for the left axis
        x.tick_params(left=False, labelleft=False)
        #remove background patch (only needed for non-white background)
        x.patch.set_visible(False)
        x.set(ylabel=axname)
    plt.savefig('output/fig3.png', dpi=300)
