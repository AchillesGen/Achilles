""" Overload the vegas integrator for more flexibility. """

import numpy as np

import vegas


class Integrator(vegas.Integrator):
    """ Overload of the vegas integrator. """

    def get_events(self, nevents):
        """ Get a fixed number of events. """

        nhcube = self.nstrat ** self.dim
        dv_y = 1. / nhcube
        neval_hcube = np.empty(nhcube, dtype=np.int)

        evals = nevents // nhcube
        remainder = nevents % nhcube

        neval_hcube[:] = evals
        neval_hcube[:remainder] += 1

        assert np.sum(neval_hcube) == nevents

        y0 = np.empty((self.dim))
        y = np.empty((nevents, self.dim))
        x = np.empty((nevents, self.dim))
        jac = np.empty((nevents))

        yran = self.ran_array_generator((nevents, self.dim))
        i_start = 0
        for ihcube in range(nhcube):
            tmp_hcube = ihcube
            for d in range(self.dim):
                y0[d] = tmp_hcube % self.nstrat
                tmp_hcube = (tmp_hcube - y0[d]) / self.nstrat
            for d in range(self.dim):
                for i in range(i_start, i_start + neval_hcube[ihcube]):
                    y[i, d] = (y0[d] + yran[i, d]) / self.nstrat
            i_start += neval_hcube[ihcube]
        self.map.map(y, x, jac, nevents)

        i_start = 0
        for ihcube in range(nhcube):
            for i in range(i_start, i_start + neval_hcube[ihcube]):
                jac[i] *= dv_y / neval_hcube[ihcube]
            i_start += neval_hcube[ihcube]
        return x, jac
