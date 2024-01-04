with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import yaml
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import pandas as pd
    import seaborn as sns
    from PIL import Image

    rename_dict = {
        'Apparli_virus': 'apparli virus',
        'DWV-B': 'DWV-B',
        'DWV-A':'DWV-A',
        'Apthili_virus': 'apthili virus',
        'BMLV':'BMLV',
        'Apis_Orthomyxovirus_1_PB1' : 'PB1',
        'Apis_Orthomyxovirus_1_GP' : 'GP',
        'Apis_Orthomyxovirus_1_PB2': 'PB2',
        'Apis_Orthomyxovirus_1_M': 'M',
        'Apis_Orthomyxovirus_1_PA': 'PA',
        'Apis_Orthomyxovirus_1_NP': 'NP'
    }

    snpslis = []
    with open('data/snps.vcf') as f:
        for line in f:
            if not line.startswith('#'):
                vir = line.strip().split()[0]
                af = float(line.strip().split()[7].split(';')[0].replace('AF=', ''))
                pos = int(line.strip().split()[1])
                if vir != 'AMFV':
                    snpslis.append(
                        [rename_dict[vir], pos, af]
                    )
    snpsdf = pd.DataFrame(snpslis)
    snpsdf.columns = ['vir', 'pos','af']
    with open('data/constellations.yaml') as f:
        a = yaml.safe_load(f)
    fig, ax = plt.subplots(
        nrows=6,
        ncols=2,
        figsize=(16,10),
        sharex=True,
        sharey=False
    )

    axix = 0
    for vir in a:
        if vir != 'apis orthomyxovirus 1':
            # SNPs & amplicon in right plot
            ampl = patches.Rectangle(
               (a[vir]['amplicon'][0][0], - 1.5),
               a[vir]['amplicon'][1][1] - a[vir]['amplicon'][0][0],
               1,
               color='#C93312',
               alpha=0.9
            )
            ax[axix,1].add_patch(ampl)
            sns.histplot(
                data=snpsdf[snpsdf['vir'] == vir],
                x='pos',
                ax=ax[axix,1],
                binwidth=50
            )
            ax[axix,1].set_ylabel('SNP count')
            ax[axix,1].set_ylim(-5,30)
            # constellation in left plot
            ax[axix,0].axhline(linewidth=0.25, y=0, color='#7f7f7f', linestyle='--')
            ax[axix,0].axhline(linewidth=0.25, y=1, color='#7f7f7f', linestyle='--')
            ax[axix,0].axhline(linewidth=0.25, y=2, color='#7f7f7f', linestyle='--')
            ax[axix,0].axhline(linewidth=0.25, y=-1, color='#7f7f7f', linestyle='--')
            ax[axix,0].axhline(linewidth=0.25, y=-2, color='#7f7f7f', linestyle='--')
            ax[axix,0].set_xlim(-100,13500)
            genrec = patches.Rectangle(
                (1,-0.125),
                a[vir]['l'],
                0.25,
                color='darkgray'
            )
            genrec2 = patches.Rectangle(
                (1,-0.125),
                a[vir]['l'],
                0.25,
                color='darkgray'
            )
            ax[axix,0].add_patch(genrec)
            ax[axix,1].add_patch(genrec2)
            for k in a[vir]:
                if k != 'l' and k != 'amplicon':
                    i1 = a[vir][k][0]
                    i2 = a[vir][k][1]
                    if i1 > i2:
                        start = i2
                        stop = i1
                    else:
                        start = i1
                        stop = i2
                    if a[vir][k][2] == '+':
                        reccol = '#C93312'
                        yanc = a[vir][k][3]-1.5
                    else:
                        reccol = '#DC863B'
                        yanc = -a[vir][k][3]+.5
                    rec = patches.Rectangle(
                        (start, yanc),
                        stop-start,
                        1,
                        color=reccol,
                        alpha = 0.9
                    )
                    ax[axix,0].add_patch(rec)
                    ax[axix,0].set_ylabel(vir)
                    ax[axix,0].set_ylim(-4,4)
                    rx, ry = rec.get_xy()
                    cx = rx
                    cy = ry + rec.get_height()+1
                    if a[vir][k][2] == '-':
                        ax[axix,0].annotate("{}\n{}".format(k, a[vir][k][4]), (cx, cy-1), color='black', weight='bold', fontsize=10, ha='left', va='top')
                    else:
                        ax[axix,0].annotate("{}\n{}".format(k, a[vir][k][4]), (cx, cy), color='black', weight='bold', fontsize=10, ha='left', va='top')
                    if vir == 'apthili virus':
                        if k == 'Capsid':
                            ax[axix,0].annotate(a[vir][k][5], (cx, 4), fontsize=8)
                        else:
                            ax[axix,0].annotate(a[vir][k][5], (cx, 3), fontsize=8)
                    if vir == 'BMLV' and a[vir][k][5] != 'NA':
                        if k == 'Cprotein':
                            ax[axix,0].annotate(a[vir][k][5], (cx, 4), fontsize=8)
                        else:
                            ax[axix,0].annotate(a[vir][k][5], (cx, 3), fontsize=8)
                    if vir in ['apparli virus', 'DWV-B']:
                        ax[axix,0].annotate(a[vir][k][5], (cx, 2), fontsize=8)
                    if vir == 'DWV-A':
                        ax[axix,0].annotate(a[vir][k][5], (cx, 3), fontsize=8)
            ax[axix,0].text(13500, 0, 'frame 1', fontsize=8, va='center', ha='center')
            ax[axix,0].text(13500, 1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[axix,0].text(13500, 2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[axix,0].text(13500, -1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[axix,0].text(13500, -2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[axix,0].set_frame_on(False)
            ax[axix,0].set_yticklabels('')
            ax[axix,0].tick_params(left = False, bottom=False)
            axix += 1
        if vir == 'apis orthomyxovirus 1':
            ax[5,0].axhline(linewidth=0.25, y=0, color='#7f7f7f', linestyle='--')
            ax[5,0].axhline(linewidth=0.25, y=1, color='#7f7f7f', linestyle='--')
            ax[5,0].axhline(linewidth=0.25, y=2, color='#7f7f7f', linestyle='--')
            ax[5,0].axhline(linewidth=0.25, y=-1, color='#7f7f7f', linestyle='--')
            ax[5,0].axhline(linewidth=0.25, y=-2, color='#7f7f7f', linestyle='--')
            ax[5,0].set_xlim(-100,13500)
            # Iterate through segments per length.
            for seg, segstart in zip(['PB2', 'PB1', 'PA', 'GP', 'NP', 'M'], [1, 2501, 5001, 7501, 9501, 11501]):
                if seg == 'PB1':
                    ampl = patches.Rectangle(
                       (segstart+a[vir][seg]['amplicon'][0][0], -1.5),
                       a[vir][seg]['amplicon'][1][1] - a[vir][seg]['amplicon'][0][0],
                       1,
                       color='blue',
                       alpha=0.9
                    )
                    ax[5,1].add_patch(ampl)
                ax[5,1].set_ylabel("SNP count")
                ax[5,1].set_ylim(-5,30)
                tdf = snpsdf.copy()
                tdf['pos'] = tdf['pos'] + segstart
                sns.histplot(
                    data=tdf[tdf['vir'] == seg],
                    x='pos',
                    ax=ax[5,1],
                    binwidth=50
                )
                #ax[axix,1].set_ylabel('SNP count')
                genrec = patches.Rectangle(
                    (segstart,-0.125),
                    a[vir][seg]['l'],
                    0.25,
                    color='darkgray'
                )
                genrec2 = patches.Rectangle(
                    (segstart,-0.125),
                    a[vir][seg]['l'],
                    0.25,
                    color='darkgray'
                )
                ax[5,0].add_patch(genrec)
                ax[5,1].add_patch(genrec2)
                for k in a[vir][seg]:
                    if k != 'l' and k != 'amplicon':
                        i1 = a[vir][seg][k][0]
                        i2 = a[vir][seg][k][1]
                        if i1 > i2:
                            start = i2+segstart
                            stop = i1+segstart
                        else:
                            start = i1+segstart
                            stop = i2+segstart
                        if a[vir][seg][k][2] == '+':
                            reccol = '#C93312'
                            yanc = a[vir][seg][k][3]-1.5
                        else:
                            reccol = '#DC863B'
                            yanc = -a[vir][seg][k][3]+.5
                        rec = patches.Rectangle(
                            (start, yanc),
                            stop-start,
                            1,
                            color=reccol,
                            alpha = 0.9
                        )
                        ax[5,0].add_patch(rec)
                        ax[5,0].set_ylabel('Apis orthomyxovirus 1')
                        ax[5,0].set_ylim(-4,4)
                        rx, ry = rec.get_xy()
                        cx = rx
                        cy = ry + rec.get_height()+1
                        if a[vir][seg][k][2] == '-':
                            ax[5,0].annotate("{}\n{}".format(k , a[vir][seg][k][4]), (cx, cy-1), color='black', weight='bold', fontsize=10, ha='left', va='top')
                        else:
                            ax[5,0].annotate("{}\n{}".format(k , a[vir][seg][k][4]), (cx, cy), color='black', weight='bold', fontsize=10, ha='left', va='top')
                        if seg in ['PB2', 'PB1', 'PA', 'GP'] and a[vir][seg][k][5] != 'NA':
                            ax[5,0].annotate(a[vir][seg][k][5], (cx, 4), fontsize=8)
                        elif seg in ['NP', 'M']:
                            ax[5,0].annotate(a[vir][seg][k][5], (cx, 3), fontsize=8)
            ax[5,0].text(13500, 0, 'frame 1', fontsize=8, va='center', ha='center')
            ax[5,0].text(13500, 1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[5,0].text(13500, 2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[5,0].text(13500, -1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[5,0].text(13500, -2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[5,0].set_frame_on(False)
            ax[5,0].set_yticklabels('')
            ax[5,0].tick_params(left = False)
    plt.savefig('output/sup_fig9.pdf', dpi=300)
    plt.savefig('output/sup_fig9.png', dpi=300)
    plt.savefig('output/sup_fig9.tiff', dpi=300)
    plt.savefig('output/sup_fig9.eps', dpi=300)
