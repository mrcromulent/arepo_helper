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
    """Class to hold a series of ArepoSnapshots, forming a run."""

    snapshots       = []
    num_snapshots   = 0
    snapbase        = None
    target_dir      = None
    cachedir        = "areporun_cache"

    run_header      = dict()
    maxima          = dict()
    minima          = dict()

    def __init__(self, snapshot_list, save_to_cache=False):
        """ArepoRun constructor.

        :param snapshot_list: List of snapshots
        :type snapshot_list: list
        :param save_to_cache: Switch that determines whether a cache should be saved to disk
        :type save_to_cache: bool

        :return: ArepoRun object
        :rtype: ArepoRun
        """

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
        """Searches all snapshots for the boxsize and records it in self.run_header."""
        snap = self.snapshots[0]
        d = [ArepoHeader.BOXSIZE]
        for key in d:
            self.run_header[key] = snap.get_from_h5(key)

    def map_maxima_minima(self):
        """Maps the maxima and minima of every quantity in utilities.plot_quantities."""

        for quantity in tqdm(plot_quantities):
            self.maxima[quantity] = []
            self.minima[quantity] = []

            for snap in self.snapshots:
                self.maxima[quantity].append(snap.max(quantity))
                self.minima[quantity].append(snap.min(quantity))

    def is_empty(self):
        """Returns if the run is empty."""
        return self.num_snapshots < 1

    def save_to_cache(self):
        """Saves the information contained in the snapshots to a picke file."""

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
        """Alternative (and more popular) constructor for ArepoRun. Gets information from a directory which contains
        snaps

        :param target_dir: Directory containing snapshots
        :type target_dir: str
        :param snapbase: Base of file name (e.g. snap for snap_000.hdf5)
        :type snapbase: str
        :param from_cache: switch that determines whether the cache should be used
        :type from_cache: bool

        :return: ArepoRun object
        :rtype: ArepoRun
        """

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
        """Sets the internal data based on unpicked data from the cache.

        :param target_dir: Directory containing cache
        :type target_dir: str

        :return: Success
        :rtype: bool
        """

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
        """Determines whether directory contains snapshots.

        :param target_dir: Directory
        :type target_dir: str
        :param snapbase: Base of snapshot file (e.g. snap for snap_000.hdf5)
        :type snapbase: str

        :return: List of paths of snapshot files in directory
        :rtype: list
        """

        file_pattern    = f"{target_dir}/{snapbase}_*.hdf5"
        snapshot_paths  = [file for file in sorted(glob.glob(file_pattern))]

        return snapshot_paths

    @staticmethod
    def find_snapbase(snapshot):
        """Gets snapbase given a snapshot object

        :param snapshot: ArepoSnapshot
        :type snapshot: ArepoSnapshot

        :return: Snap base
        :rtype: str
        """
        snap_filename = snapshot.filename
        return re.sub(r"(?=_).+$", "", os.path.basename(snap_filename))

    @staticmethod
    def find_target_dir(snapshot):
        """Gets directory given a snapshot object

        :param snapshot: ArepoSnapshot
        :type snapshot: ArepoSnapshot

        :return: Directory
        :rtype: str
        """
        return os.path.dirname(snapshot.filename)

    def __str__(self):
        if not self.is_empty():
            return f"{[str(sh) for sh in self.snapshots]}"

        else:
            return "Empty run."
