from h5_file import ArepoH5File


class ArepoSnapshot(ArepoH5File):

    def __init__(self, filename):
        # self.data = dict()
        super(ArepoSnapshot, self).__init__(filename)

    def __str__(self):
        return f"ArepoSnapshot with n particles: {self.num_particles()}"
