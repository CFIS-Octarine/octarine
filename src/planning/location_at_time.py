__author__ = 'Michele Bannister   git:@mtbannister'

import os
from mp_ephem import BKOrbit, EphemerisReader

date = '2016-10-24T16:00:00'
ossinpath = '/Users/bannisterm/Dropbox/OSSOS/measure3/ossin/tmp/' #'vos:OSSOS/dbaseclone/ast/'  #
outfile = '/Users/bannisterm/Desktop/{}.txt'.format(date)

with open(outfile, 'w') as ofile:
    ofile.write('Target	RA (hrs)	DEC		m_r	delta RA (")	delta DEC (")	Time predicted\n')

for kbo_filename in os.listdir(ossinpath):
    mpc_observations = EphemerisReader.read(ossinpath + kbo_filename).mpc_observations
    orbit = BKOrbit(mpc_observations)
    orbit.predict(date=date)

    with open(outfile, 'a') as ofile:
        ofile.write("{:>10s} {:>10s} {:6.2f} {:6.2f} {:>10s}\n".format(
            orbit.name, orbit.coordinate.to_string('hmsdms', sep=':'),
            # orbit.coordinate.dec.to_string('hmsdms'),  # need to add mag back in
            orbit.dra, orbit.ddec, date))
