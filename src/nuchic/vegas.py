""" Overload the vegas integrator for more flexibility. """

import numpy as np

import vegas
import gvar


class Integrator(vegas.Integrator):
    """ Overload of the vegas integrator. """

    def __init__(self, inmap, **kwargs):
        super().__init__(inmap, **kwargs)
        self.maximum = 0

    def _map(self, nevents, nhcube, neval_hcube):

        y_0 = np.empty((self.dim))
        y = np.empty((nevents, self.dim))
        x = np.empty((nevents, self.dim))
        jac = np.empty((nevents))

        yran = self.ran_array_generator((nevents, self.dim))
        i_start = 0
        for ihcube in range(nhcube):
            tmp_hcube = ihcube
            for dim in range(self.dim):
                y_0[dim] = tmp_hcube % self.nstrat
                tmp_hcube = (tmp_hcube - y_0[dim]) / self.nstrat
            for dim in range(self.dim):
                for i in range(i_start, i_start + neval_hcube[ihcube]):
                    y[i, dim] = (y_0[dim] + yran[i, dim]) / self.nstrat
            i_start += neval_hcube[ihcube]
        self.map.map(y, x, jac, nevents)

        return x, jac

    def generate(self, nevents):
        """ Generate a fixed number of points with a weight. """
        nhcube = self.nstrat ** self.dim
        dv_y = 1. / nhcube
        neval_hcube = np.empty(nhcube, dtype=np.int)

        evals = nevents // nhcube
        remainder = nevents % nhcube

        neval_hcube[:] = evals
        if remainder != 0:
            neval_hcube[:remainder] += 1

        assert np.sum(neval_hcube) == nevents

        x, jac = self._map(nevents, nhcube, neval_hcube)

        hcube_array = np.empty(nevents, np.int)

        i_start = 0
        for ihcube in range(nhcube):
            for i in range(i_start, i_start + neval_hcube[ihcube]):
                jac[i] *= dv_y / neval_hcube[ihcube]
                hcube_array[i] = ihcube
            i_start += neval_hcube[ihcube]
        return x, jac, hcube_array

    def get_weighted_events(self, func, nevents):
        """ Get a fixed number of weighted events. """
        if self.sync_ran and self.mpi:
            self.synchronize_random()

        self.set(beta=-1, neval=nevents)

        point, weight, hcube = self.generate(nevents)

        nwf = np.array(np.zeros_like(weight), np.float)
        sum_wf = np.array(np.zeros_like(weight), np.float)
        sum_wf2 = np.array(np.zeros_like(weight), np.float)
        for i, _ in enumerate(point):
            nwf[hcube[i]] += 1
            weight[i] = func(point[i]) * weight[i]
            sum_wf[hcube[i]] += weight[i]
            sum_wf2[hcube[i]] += weight[i]**2

        integral = 0
        variance = 0
        for i, _ in enumerate(sum_wf):
            integral += sum_wf[i]
            variance += (sum_wf2[i] * nwf[i] - sum_wf[i] ** 2) / (nwf[i] - 1.)

        return point, weight, gvar.gvar(integral, variance ** 0.5)

    def get_unweighted_events(self, func, nevents):
        """ Generate a fixed number of unweighted events. """
        if self.sync_ran and self.mpi:
            self.synchronize_random()

        self.set(beta=-1, neval=nevents)

        # Estimate acceptance
        accuracy = nevents ** -0.5
        eta = 1
        eta_0 = 0
        weights = np.array([])
        while abs(eta - eta_0) / eta > 0.001:
            nsamples = int(np.ceil(accuracy ** -2 / eta) - len(weights))
            if nsamples < 0:
                break
            _, weight, _ = self.get_weighted_events(func, nsamples)
            weights = np.sort(np.append(weights, weight))
            cum_weights = np.cumsum(weights)
            cum_weights /= cum_weights[-1]
            index = np.searchsorted(cum_weights, 1. - accuracy)
            print(cum_weights, index, 1. - accuracy)
            self.maximum = weights[index]
            avg_val = np.mean(weights[:index])
            eta_0 = eta
            eta = avg_val / self.maximum
            print(eta_0, eta)
            print(cum_weights[index])
            print(self.maximum)
