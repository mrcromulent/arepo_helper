from utilities import common_snapbases, plot_quantities
from h5_file import ArepoSnapshot
from names import ArepoHeader
from datetime import datetime
from tqdm import tqdm
import pickle as pk
import glob
import os
import re


class ArepoRun(object):

    snapshots       = []
    num_snapshots   = 0
    snapbase        = None
    target_dir      = None
    cachedir        = "areporun_cache"

    run_header      = dict()
    maxima          = dict()
    minima          = dict()

    def __init__(self, snapshot_list, save_to_cache=False):

        self.snapshots      = snapshot_list
        self.num_snapshots  = len(snapshot_list)
        
        if self.num_snapshots > 0:
            self.snapbase       = self.find_snapbase(snapshot_list[0])
            self.target_dir     = self.find_target_dir(snapshot_list[0])

            if save_to_cache:
                self.map_header()
                self.map_maxima_minima()
                self.save_to_cache()

    def map_header(self):
        snap = self.snapshots[0]
        d = [ArepoHeader.BOXSIZE]
        for key in d:
            self.run_header[key] = snap.get_from_h5(key)

    def map_maxima_minima(self):

        for quantity in tqdm(plot_quantities):
            self.maxima[quantity] = []
            self.minima[quantity] = []

            for snap in self.snapshots:
                self.maxima[quantity].append(snap.max(quantity))
                self.minima[quantity].append(snap.min(quantity))

    def is_empty(self):
        return self.num_snapshots < 1

    def save_to_cache(self):

        dt      = datetime.now().strftime("%Y%m%d-%H%M%S")
        direc   = f"{self.target_dir}/{self.cachedir}/"
        path    = f"{direc}{dt}.obj"

        if not os.path.exists(direc):
            os.mkdir(direc)

        with open(path, "wb") as handle:
            pk.dump((self.snapshots, self.run_header, self.maxima, self.minima), handle)

        print(f"Cache saved at: {path}")

    @classmethod
    def from_directory(cls, target_dir=".", snapbase=None, from_cache=True):

        target_dir      = os.path.abspath(target_dir)
        snapshot_paths  = []

        if snapbase is not None:
            snapshot_paths = cls.check_directory_for_snapshots(target_dir, snapbase)
        else:
            for sb in common_snapbases:
                snapshot_paths = cls.check_directory_for_snapshots(target_dir, sb)
                if len(snapshot_paths) > 0:
                    snapbase = sb
                    break

        if len(snapshot_paths) > 0:

            if from_cache: 
                from_cache = cls.open_run_from_cache(target_dir)
                
            if not from_cache:
                cls.snapshots = []
                for path in snapshot_paths:
                    cls.snapshots.append(ArepoSnapshot(path))

            save_to_cache = not from_cache
            return cls(cls.snapshots, save_to_cache=save_to_cache)
        else:
            raise FileNotFoundError(f"No run found in {target_dir} with snapbase {snapbase}.")
    
    @classmethod
    def open_run_from_cache(cls, target_dir):

        list_of_files   = glob.glob(f"{target_dir}/{cls.cachedir}/*.obj")

        if len(list_of_files) > 0:
            path = max(list_of_files, key=os.path.getctime)
            with open(path, "rb") as handle:
                cls.snapshots, cls.run_header, cls.maxima, cls.minima = pk.load(handle)

            print(f"Cache found at: {path}")
            return True

        else:
            print("Cache requested and none found.")
            return False

    @staticmethod
    def check_directory_for_snapshots(target_dir, snapbase):

        file_pattern    = f"{target_dir}/{snapbase}_*.hdf5"
        snapshot_paths  = [file for file in sorted(glob.glob(file_pattern))]

        return snapshot_paths

    @staticmethod
    def find_snapbase(snapshot):
        snap_filename = snapshot.filename
        return re.sub(r"(?=_).+$", "", os.path.basename(snap_filename))

    @staticmethod
    def find_target_dir(snapshot):
        return os.path.dirname(snapshot.filename)

    def __str__(self):
        if not self.is_empty():
            return f"{[str(sh) for sh in self.snapshots]}"

        else:
            return "Empty run."
