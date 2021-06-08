import matplotlib.pyplot as plt
import numpy as np
import h5py

import pcolor_pierre
from names import n
import snapshot

s = snapshot.ArepoSnapshot("./data/healpixtext_snaps/snapshot_000.hdf5")

coords = s.get_from_h5(n.COORDINATES)
points = s.get_from_h5(n.DENSITY)

data = pcolor_pierre.make_pcolor(coords.astype('float64'),
                                points.astype('float64'),
                                1024, 1024,
                                1e10, 1e10, 0,
                                5e9, 5e9, 5e9,
                                axis0=0, axis1=1,
                                nz=1,
                                include_neighbours_in_output=True,
                                numthreads=1)

if __name__ == '__main__':
    pass
