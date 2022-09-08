import os
import numpy as np
from .ramses import Ramses

class Dyn(Ramses):

    def read_sink_as_msun_au(self, fname):
        """ Read sink_x.csv as the following unit system (and move to com):
        [M] = 1 Msun
        [L] = 1 AU
        [T] = sqrt(AU**3*(1/(G*Msun))) = 5022642.892 s = 0.1591579490 yr
        [V] = [L]/[T] = 29.78469183 km / s

        time, id, mass, x, y, z, vx, vy, vz
        """

        assert os.path.isfile(fname)
        data = np.loadtxt(fname, delimiter=',')
        particles = data[:, 1:8]

        particles[:, 0] *= self.unit_m / cgs.Msun
        particles[:, 1:4] *= self.unit_l_code / cgs.AU
        particles[:, 4:7] *= self.unit_v / 1e5 # km/s
        
        # shift to the center-of-mass frame with the com in the origin
        Mtot = particles[:, 0].sum()
        x_com = np.matmul(particles[:, 0], particles[:, 1:4]) / Mtot
        v_com = np.matmul(particles[:, 0], particles[:, 4:7]) / Mtot
        for dim in range(3):
            particles[:, 1+dim] -= x_com[dim]
            particles[:, 4+dim] -= v_com[dim]
            
        return particles
        
