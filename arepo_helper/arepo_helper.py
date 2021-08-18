def main():
    """Main script."""

    from names import n
    from run import ArepoRun
    from plot_manager import quick_pcolor, quick_radial

    ar = ArepoRun.from_directory("/home/pierre/Desktop/output_higher_basefac")
    bs = ar.snapshots[0].get_from_h5(n.BOXSIZE)
    ibs = 1e10
    a = [bs / 2 - ibs / 2, bs / 2, bs / 2]
    b = [bs / 2 + ibs / 2, bs / 2, bs / 2]

    for i, s in enumerate(ar.snapshots):
        plot = quick_pcolor(s, n.DENSITY, explicit_options={"cbar_lims": [1e-1, 1e7], "inner_boxsize": 1e10})
        plot.save(f"{i}.png")

        rad = quick_radial(s, n.DENSITY, a=a, b=b, explicit_options={"inner_boxsize": ibs})
        rad.save(f"{i}_rad.png")


if __name__ == "__main__":
    main()
