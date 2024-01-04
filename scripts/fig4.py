with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from PIL import Image
    #mlpa
    mlpa = mpimg.imread('data/figures/mlpa.tif')
    tem1 = mpimg.imread('data/figures/tem_1.jpg')
    tem2 = mpimg.imread('data/figures/tem_2.jpg')
    fig = plt.figure(constrained_layout=True, figsize=(10,8))
    gs = fig.add_gridspec(4,4)
    ax_mlpadata = fig.add_subplot(gs[0, 0:4])
    ax_tem1 = fig.add_subplot(gs[1:3, 0:2])
    ax_tem2 = fig.add_subplot(gs[1:3, 2:4])
    ax_mlpadata.axis('off')
    ax_mlpadata.set_adjustable('box')
    ax_tem1.axis('off')
    ax_tem2.axis('off')
    ax_mlpadata.imshow(mlpa)
    # ax_tem1.imshow(tem1, aspect='auto')
    # ax_tem2.imshow(tem2, aspect='auto')
    ax_tem1.imshow(tem1)
    ax_tem2.imshow(tem2)
    ax_mlpadata.text(-0.025, 1.05, 'A', transform=ax_mlpadata.transAxes,
                        size=20, weight='bold')
    ax_mlpadata.text(0.15, 1.05, 'B', transform=ax_mlpadata.transAxes,
                        size=20, weight='bold')
    ax_mlpadata.text(-0.025, -0.8, 'C', transform=ax_mlpadata.transAxes,
                        size=20, weight='bold')
    ax_mlpadata.text(0.5, -0.8, 'D', transform=ax_mlpadata.transAxes,
                        size=20, weight='bold')
    # - strand.
    ax_mlpadata.text(
        0.025, -0.1,
        'L',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.06, -0.1,
        '+',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.11, -0.1,
        '-',
        transform=ax_mlpadata.transAxes,
        size=15
    )

    # MLPA
    ax_mlpadata.text(
        0.21, -0.1,
        'L',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.27, -0.1,
        'e',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.335, -0.1,
        'l',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.39, -0.1,
        'l',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.435, -0.1,
        'sl',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.485, -0.1,
        'wp',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.54, -0.1,
        'wp',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.6, -0.1,
        'rp',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.66, -0.1,
        'eb',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.725, -0.1,
        'v',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.77, -0.1,
        'pp',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.83, -0.1,
        'ab',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.9, -0.1,
        '-',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    ax_mlpadata.text(
        0.96, -0.1,
        'L',
        transform=ax_mlpadata.transAxes,
        size=15
    )
    # Arrows.
    ax_tem1.annotate("", xytext=(2350, 50), xy=(2200, 300),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))
    ax_tem1.annotate("", xytext=(2600, 100), xy=(2450, 350),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))
    ax_tem1.annotate("", xytext=(2850, 280), xy=(2700, 530),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))

    ax_tem2.annotate("", xytext=(1500, 650), xy=(1350, 900),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))
    ax_tem2.annotate("", xytext=(1950, 750), xy=(1800, 1000),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))
    ax_tem2.annotate("", xytext=(2400, 1050), xy=(2250, 1300),
                arrowprops=dict(linewidth=3,arrowstyle="->", color='#79402E'))
    plt.savefig('output/fig4.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('output/fig4.png', dpi=300, bbox_inches='tight')
    plt.savefig('output/fig4.tiff', dpi=300, bbox_inches='tight')
    plt.savefig('output/fig4.eps', dpi=300, bbox_inches='tight')
