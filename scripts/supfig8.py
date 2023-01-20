with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import yaml
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches


    with open('data/constellations.yaml') as f:
        a = yaml.safe_load(f)
    fig, ax = plt.subplots(6, figsize=(8,10), sharex=True, sharey=True)
    #patches.Rectangle()
    vir = 'Apthili virus' 
    axix = 0
    for vir in a:
        if vir != 'Apis Orthomyxovirus 1':
            ax[axix].axhline(linewidth=0.25, y=0, color='#7f7f7f', linestyle='--')
            ax[axix].axhline(linewidth=0.25, y=1, color='#7f7f7f', linestyle='--')
            ax[axix].axhline(linewidth=0.25, y=2, color='#7f7f7f', linestyle='--')
            ax[axix].axhline(linewidth=0.25, y=-1, color='#7f7f7f', linestyle='--')
            ax[axix].axhline(linewidth=0.25, y=-2, color='#7f7f7f', linestyle='--')
            ax[axix].set_xlim(-100,13500)
            genrec = patches.Rectangle(
                (1,-0.125),
                a[vir]['l'],
                0.25,
                color='darkgray'
            )
            ax[axix].add_patch(genrec)
            for k in a[vir]:
                if k != 'l':
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
                    ax[axix].add_patch(rec)
                    ax[axix].set_ylabel(vir)
                    ax[axix].set_ylim(-4,4)
                    rx, ry = rec.get_xy()
                    cx = rx
                    cy = ry + rec.get_height()+1
                    if a[vir][k][2] == '-':
                        ax[axix].annotate("{}\n{}".format(k , a[vir][k][4]), (cx, cy-1), color='black', weight='bold', fontsize=10, ha='left', va='top')
                    else:
                        ax[axix].annotate("{}\n{}".format(k , a[vir][k][4]), (cx, cy), color='black', weight='bold', fontsize=10, ha='left', va='top')
            ax[axix].text(13500, 0, 'frame 1', fontsize=8, va='center', ha='center')
            ax[axix].text(13500, 1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[axix].text(13500, 2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[axix].text(13500, -1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[axix].text(13500, -2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[axix].set_frame_on(False)
            ax[axix].set_yticklabels('')
            ax[axix].tick_params(left = False, bottom=False)
            axix += 1
        if vir == 'Apis Orthomyxovirus 1':
            ax[5].axhline(linewidth=0.25, y=0, color='#7f7f7f', linestyle='--')
            ax[5].axhline(linewidth=0.25, y=1, color='#7f7f7f', linestyle='--')
            ax[5].axhline(linewidth=0.25, y=2, color='#7f7f7f', linestyle='--')
            ax[5].axhline(linewidth=0.25, y=-1, color='#7f7f7f', linestyle='--')
            ax[5].axhline(linewidth=0.25, y=-2, color='#7f7f7f', linestyle='--')
            ax[5].set_xlim(-100,13500)
            # Iterate through segments per length.
            for seg, segstart in zip(['PB2', 'PB1', 'PA', 'GP', 'NP', 'M'], [1, 2501, 5001, 7501, 9501, 11501]):
                genrec = patches.Rectangle(
                    (segstart,-0.125),
                    a[vir][seg]['l'],
                    0.25,
                    color='darkgray'
                )
                ax[5].add_patch(genrec)
                for k in a[vir][seg]:
                    if k != 'l':
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
                        ax[5].add_patch(rec)
                        ax[5].set_ylabel(vir)
                        ax[5].set_ylim(-4,4)
                        rx, ry = rec.get_xy()
                        cx = rx
                        cy = ry + rec.get_height()+1
                        if a[vir][seg][k][2] == '-':
                            ax[5].annotate("{}\n{}".format(k , a[vir][seg][k][4]), (cx, cy-1), color='black', weight='bold', fontsize=10, ha='left', va='top')
                        else:
                            ax[5].annotate("{}\n{}".format(k , a[vir][seg][k][4]), (cx, cy), color='black', weight='bold', fontsize=10, ha='left', va='top')
            ax[5].text(13500, 0, 'frame 1', fontsize=8, va='center', ha='center')
            ax[5].text(13500, 1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[5].text(13500, 2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[5].text(13500, -1, 'frame 2', fontsize=8, va='center', ha='center')
            ax[5].text(13500, -2, 'frame 3', fontsize=8, va='center', ha='center')
            ax[5].set_frame_on(False)
            ax[5].set_yticklabels('')
            ax[5].tick_params(left = False)
    plt.savefig('output/sup_fig8.png', dpi=300)
