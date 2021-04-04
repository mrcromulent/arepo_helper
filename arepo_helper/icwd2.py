from createICs import create_particles_healpix
from species import ArepoSpeciesList
from pyhelm_eos import loadhelm_eos
from ics import ArepoICs
from scipy import interpolate
from wd_utils import WDUtils
from ic import create_wd
from const import msol
import numpy as np
import utilities
from collections import OrderedDict
from names import n
import scipy.optimize as opt

spec_file = "./data/eostable/species05.txt"
helm_file = "./data/snapshots/helm_table.dat"
wdmass = 0.55 * msol
pmass = 1e-6 * msol
boxsize = 1e10
c12prop = 0.5
o16prop = 1 - c12prop

def wdCOgetMassFromRhoCExact( rhoc, eos, mass=0., temp=5e5 ):
    import ic
    wd = ic.create_wd( eos, rhoc, temp=temp, xC12=0.5, xO16=0.5 )
    print( rhoc, wd['dm'].sum()/msol )
    return wd['dm'].sum()/msol - mass

def wdCOgetRhoCFromMassExact( mass, eos, temp=5e5, xtol=1e2 ):
    rhoguess = 1e4
    massguess = wdCOgetMassFromRhoCExact( rhoguess, eos, temp ) * msol

    if massguess > mass:
      print( "Already a central density of 1e4 g/ccm produces a WD more massive than %g msun (%g)." % (mass/msol,massguess/msol) )

    while massguess <= mass:
      rhoguess *= 10.
      massguess = wdCOgetMassFromRhoCExact( rhoguess, eos, temp=temp ) * msol

    print( "Guess for central density: %g g/ccm." % rhoguess )
    rhoc = opt.bisect( wdCOgetMassFromRhoCExact, 0.1 * rhoguess, rhoguess, args=(eos,mass/msol,temp,), xtol=xtol )
    print( "Actual central density for %g msun WD: %g g/ccm." % (mass/msol,rhoc) )
    return rhoc

eos = loadhelm_eos(helm_file, spec_file, True)
# rhoc = WDUtils.get_rho_c_from_mass(wdmass / msol, eos)
rhoc = wdCOgetRhoCFromMassExact(wdmass, eos)

print(rhoc)

wd = create_wd(eos, rhoc, xC12=c12prop, xO16=o16prop)
frho = interpolate.interp1d(wd['r'], wd['rho'], kind='cubic')
fpres = interpolate.interp1d(wd['r'], wd['p'], kind='cubic')

sp = ArepoSpeciesList(spec_file)
wd['v'] = np.zeros(wd['ncells'])
wd['xnuc'] = np.zeros(sp.num_species)
wd['xnuc'][sp.index_of('C12')] = c12prop
wd['xnuc'][sp.index_of('O16')] = o16prop
wd['count'] = wd['ncells']

# Map the 1d distribution into a set of 3d points using a healpix mapping
data = create_particles_healpix(wd, eos,
                                minenergy=1e14,
                                boxfactor=10.,
                                boxres=32,
                                npart=0,
                                randomizeshells=True,
                                randomizeradii=False,
                                makebox=True,
                                nspecies=sp.num_species,
                                boxsize=boxsize,
                                pmass=pmass)

# Find the radius of each cell using the distance from the centre of the box
rad = np.sqrt(((data['pos'] - 0.5 * boxsize) ** 2.).sum(axis=1))
i, = np.where(rad < wd['r'].max())
rho = frho(rad[i])
pres = fpres(rad[i])
xnuc = np.zeros(sp.num_species)
xnuc[sp.index_of('C12')] = c12prop
xnuc[sp.index_of('O16')] = o16prop

# For the cells included in the WD, find the internal energy and temperature
for index in range(np.size(i)):
    idx = i[index]
    temp, data['u'][idx] = eos.pgiven(rho[index], xnuc, pres[index])

# Set mass equal to density in the WD, zero everywhere else
data['mass'][:] = 0.
data['mass'][i] = rho


def gadget_write_ics(filename, data, format='hdf5', transpose=True, time=0., skipxnuc=False, double=False,
                     longids=False, num_files=1, boxsize=0., verbose=False, masses=None):
    if format == 'gadget1':
        # gadget_write_ics_format1(filename, data, transpose, time, skipxnuc, double, longids, num_files, boxsize,
        #                          verbose, masses)
        pass
    elif format == 'gadget2':
        # gadget_write_ics_format2(filename, data, transpose, time, skipxnuc, double, longids, num_files, boxsize,
        #                          verbose, masses)
        pass
    elif format == 'hdf5':
        gadget_write_ics_format3(filename, data, time, double, longids, num_files, boxsize, masses)
    else:
        raise ValueError('Choose a valid file format')


def gadget_write_ics_format3(filename, data, time, double, longids, num_files, boxsize, masses, skipxnuc=False):
    import h5py

    filename += '.hdf5'
    print("Writing gadget file: ", filename)
    f = h5py.File(filename, 'w')

    npart = np.zeros(6, dtype=np.int32)
    nparthighword = np.zeros(6, dtype=np.int32)
    offset = np.zeros(6, dtype=np.int32)
    npartmass = np.zeros(6, dtype=np.int32)
    massoffset = np.zeros(6, dtype=np.int32)

    if 'type' in data:
        for ptype in range(6):
            npart[ptype] = np.size(np.where(data['type'] == ptype))
    else:
        npart[0] = data['count']

    offset[1:] = np.cumsum(npart[:-1])

    if not masses is None:
        massarr = masses
    else:
        massarr = np.zeros(6, dtype=np.float64)

    npartmass[:] = npart[:]
    j, = np.where(massarr > 0.0)
    npartmass[j] = 0
    massoffset[1:] = np.cumsum(npartmass[:-1])

    header = f.create_group("/Header")

    header.attrs['NumPart_ThisFile'] = npart
    header.attrs['NumPart_Total'] = npart
    header.attrs['NumPart_Total_HighWord'] = nparthighword
    header.attrs['MassTable'] = massarr
    header.attrs['Time'] = time
    header.attrs['NumFilesPerSnapshot'] = np.array(num_files, dtype=np.int32)
    header.attrs['BoxSize'] = boxsize
    header.attrs['Flag_DoublePrecision'] = np.array(double, dtype=np.int32)
    header.attrs['Flag_IC_info'] = np.array(0, dtype=np.int32)
    header.attrs['Flag_Entropy_ICs'] = np.array(0, dtype=np.int32)
    header.attrs['Redshift'] = np.array(0, dtype=np.float64)
    header.attrs['Omega0'] = np.array(0, dtype=np.float64)
    header.attrs['OmegaLambda'] = np.array(0, dtype=np.float64)
    header.attrs['HubbleParam'] = np.array(0, dtype=np.float64)
    header.attrs['Flag_Sfr'] = np.array(0, dtype=np.int32)
    header.attrs['Flag_Cooling'] = np.array(0, dtype=np.int32)
    header.attrs['Flag_StellarAge'] = np.array(0, dtype=np.int32)
    header.attrs['Flag_Metals'] = np.array(0, dtype=np.int32)
    header.attrs['Flag_Feedback'] = np.array(0, dtype=np.int32)

    if double:
        dtype = np.float64
    else:
        dtype = np.float32

    if longids:
        dtypeids = np.uint64
    else:
        dtypeids = np.uint32

    fields = []
    if 'bfld' in data:
        fields += ['bfld']
    if not skipxnuc and ('xnuc' in data):
        fields += ['xnuc']
    if 'temp' in data:
        fields += ['temp']
    if 'rho' in data:
        fields += ["rho"]
    if 'pass' in data:
        fields += ["pass"]
    if 'erad' in data:
        fields += ["erad"]

    fields_to_names = {
        'bfld': 'MagneticField',
        'xnuc': 'NuclearComposition',
        'temp': 'Temperature',
        'rho': 'Density',
        'pass': 'PassiveScalars',
        'erad': 'Erad'
    }

    for ptype in range(6):
        if npart[ptype] > 0:
            group = f.create_group('/PartType%d' % ptype)

            pos = group.create_dataset("Coordinates", (npart[ptype], 3), dtype)
            pos[:, :] = data['pos'][offset[ptype]:offset[ptype] + npart[ptype], :]

            vel = group.create_dataset("Velocities", (npart[ptype], 3), dtype)
            vel[:, :] = data['vel'][offset[ptype]:offset[ptype] + npart[ptype], :]

            id = group.create_dataset("ParticleIDs", (npart[ptype],), dtypeids)
            if 'id' in data:
                id[:] = data['id'][offset[ptype]:offset[ptype] + npart[ptype]]
            else:
                id[:] = np.arange(offset[ptype] + 1, offset[ptype] + npart[ptype] + 1, dtype=dtypeids)

            if massarr[ptype] == 0:
                mass = group.create_dataset("Masses", (npart[ptype],), dtype)
                mass[:] = data['mass'][massoffset[ptype]:massoffset[ptype] + npart[ptype]]

            if ptype == 0:
                u = group.create_dataset("InternalEnergy", (npart[ptype],), dtype)
                u[:] = data['u'][offset[ptype]:offset[ptype] + npart[ptype]]
                # loop over the other fields
                for field in fields:
                    hdf5key = fields_to_names[field]
                    if field in ['bfld', 'xnuc', 'pass']:
                        val = group.create_dataset(hdf5key, (npart[ptype], data[field].shape[1]), dtype)
                        val[:, :] = data[field][offset[ptype]:offset[ptype] + npart[ptype], :]
                    else:
                        val = group.create_dataset(hdf5key, (npart[ptype],), dtype)
                        val[:] = data[field][offset[ptype]:offset[ptype] + npart[ptype]]
    f.close()

    print("Done.")
    return


gadget_write_ics("bin.dat.ic", data, double=True)


# Write initial condition file
# a = ArepoICs('bin.dat.ic.hdf5')
#
# nparticles = data['mass'].shape[0]
# data["id"] = np.array(range(1, nparticles + 1))
#
#
# particle_dict = OrderedDict({"PartType0": utilities.convert_from_ruediger_dict(data)})
# particle_list = np.array([nparticles, 0, 0, 0, 0, 0], dtype=np.int32)
# header = {n.NUMPARTTOTAL: particle_list,
#           n.NUMPARTTHISFILE: particle_list,
#           n.NUMPARTTOTALHIGHWORD: np.zeros(6, dtype=np.int32),
#           n.FLAGDOUBLEPRECISION: np.array(True, dtype=np.int32),
#           n.BOXSIZE: boxsize}
# a.write_ics(particle_dict, header)

if __name__ == '__main__':
    pass
