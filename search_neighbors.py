import numpy as np
import sys
import os
from astroquery.simbad import Simbad
from astropy import coordinates
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.utils.console import color_print
import astropy.units as u
import urllib2

if len(sys.argv) == 1:
    print 'Usage:'
    print 'python search_neighbors.py <camp> <epic id> <arcmin>'
    print 'El catalogo .csv debe estar en la carpeta, de lo contrario se descargara'

camp = int(sys.argv[1])
epic = int(sys.argv[2])
dist = int(sys.argv[3]) * u.arcminute

#catalog = 'K2Campaign%dtargets.csv' % camp
color_print('\nBuscando catalogo EPIC', 'yellow')
catalog = 'K2Campaign%dtargets.csv' % camp
if not os.path.isfile(catalog):
    print 'Descargando...'
    catalog = urllib2.urlopen('http://keplerscience.arc.nasa.gov/K2/docs/Campaigns/C%d/K2Campaign%dtargets.csv' % (camp,camp))
else:
    print 'Encontrada copia local'

data = np.genfromtxt(catalog, delimiter=',', usecols=range(4))
nans = np.isfinite(data[:,1])
data = data[nans]
epid, ra, de, kp = data.T

idx = epid == epic
if idx.sum() == 0:
    color_print('\nEPIC ID no encontrado!', 'lightred')
    sys.exit()
ora = ra[idx][0] * u.degree
ode = de[idx][0] * u.degree

coords = coordinates.SkyCoord(ora, ode, unit=u.degree)
color_print('\nObjeto', 'yellow')
color_print('EPIC %d ' % int(epid[idx][0]), 'red', coords.to_string(style='hmsdms', sep=' '))
ocoord = coordinates.SkyCoord(ra[~idx], de[~idx], unit=u.degree)
osep   = ocoord.separation(coords).arcminute

cerca  = osep < dist.value
osep   = osep[cerca] * u.arcminute
epicc  = epid[~idx][cerca].astype(int)
rac    = Angle(ra[~idx][cerca], unit=u.deg).to_string(unit=u.hourangle, sep=' ', pad=True, precision=3)
dec    = Angle(de[~idx][cerca], unit=u.deg).to_string(unit=u.deg, sep=' ', pad=True, alwayssign=True, precision=2)
kpc    = kp[~idx][cerca] * u.mag

k2 = Table([epicc, rac, dec, kpc, osep], names=('EPIC ID', 'RA', 'DEC', 'Kp', 'Separation'))
k2['RA'].unit = '"h:m:s"'
k2['DEC'].unit = '"d:m:s"'
k2.sort('Separation')
color_print('\nTargets cercanos en el EPIC catalogue', 'yellow')
k2.pprint()

sim = Simbad()
sim.remove_votable_fields('coordinates')
sim.add_votable_fields('ra', 'dec', 'otype', 'flux(V)', 'flux(K)', 'sp')


table = sim.query_region(coords, radius=dist)
allcoords = coordinates.SkyCoord(table['RA'], table['DEC'], unit=(u.hourangle, u.deg))

table['Separation'] = allcoords.separation(coords).arcminute
table['Separation'].unit = u.arcminute
color_print('\nTargets cercanos en SIMBAD', 'yellow')
table.pprint(max_lines=-1, max_width=-1)
