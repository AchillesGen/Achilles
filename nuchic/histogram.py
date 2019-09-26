""" Implement the histogram class """

import os
import numpy as np
import matplotlib.pyplot as plt


class Histogram:
    """ Histogram class to generate plots for the nuChic code. """
    def __init__(self, plot_range=None, bins=100, scale='linear', copy=None):
        self.scale = scale
        self.range = None
        self.bins = None
        self.weights = None
        self.weights2 = None
        self.flow = np.array([0, 0])
        if isinstance(copy, Histogram):
            self._copy_construct(copy)
        else:
            self._construct(plot_range, bins, scale)

    def _copy_construct(self, other):
        self.scale = other.scale
        self.range = other.range
        self.bins = other.bins
        self.weights = other.weights
        self.weights2 = other.weights2
        self.flow = other.flow

    def _construct(self, plot_range, bins, scale):
        if plot_range is not None:
            if not isinstance(plot_range, list):
                raise TypeError('Must pass a list with the min and max value '
                                'of the histogram.')
            if len(plot_range) != 2:
                raise ValueError('Range for the histogram must be equal to '
                                 'two. Got a length of {}'.format(
                                     len(plot_range)))
            if not isinstance(bins, int):
                raise TypeError('Must pass integer number of bins if given '
                                'min and max values.')
            self.range = plot_range
            if self.range[0] > self.range[1]:
                raise ValueError('Require range[0] > range[1], '
                                 'got {} and {}'
                                 .format(self.range[0], self.range[1]))
            if scale == 'linear':
                self.bins = np.linspace(self.min, self.max, bins)
            elif scale == 'log':
                if self.min == 0:
                    raise ValueError('Logarithmic requires min_val > 0, '
                                     'got a value of {}'.format(self.min))
                self.bins = np.logspace(np.log10(self.min),
                                        np.log10(self.max),
                                        bins)
            else:
                raise ValueError('Only `linear` and `log` binning exist, '
                                 'got {} binning.'.format(scale))
        elif isinstance(bins, (list, np.ndarray)):
            self.bins = np.array(bins)
            self.range = [bins[0], bins[-1]]
            if np.any(bins[:-1] >= bins[1:]):
                raise ValueError('Bins must be in increasing order, '
                                 'and have non-zero width.')
        else:
            raise ValueError('Either min_val and max_val need to be defined, '
                             'or a bin array needs to be supplied.')

        self.weights = np.array(np.zeros_like(self.bin_centers))
        self.weights2 = np.array(np.zeros_like(self.weights))
        self.flow = np.array([0, 0])

    @property
    def bin_centers(self):
        """ Get the bin centers. """
        return (self.bins[1:] + self.bins[:-1])/2.

    @property
    def min(self):
        """ Get minimum value of histogram. """
        return self.range[0]

    @property
    def max(self):
        """ Get maximum value of histogram. """
        return self.range[1]

    def fill(self, value, weight=1.):
        """ Add a new event fill to the histogram. """
        loc = np.searchsorted(self.bins, value)
        if loc == len(self.bins):
            self.flow[1] += weight
        elif loc == 0:
            self.flow[0] += weight
        else:
            self.weights[loc-1] += weight
            self.weights2[loc-1] += weight**2

    def __str__(self):
        """ Convert histogram to a string. """
        line_format = '{:^15}{:^15}{:^15}{:^15}\n'
        string = line_format.format('left', 'right', 'wgt', 'wgt2')
        for i, _ in enumerate(self.bins[1:]):
            string += line_format.format(self.bins[i],
                                         self.bins[i+1],
                                         self.weights[i],
                                         self.weights2[i])
        return string

    @property
    def underflow(self):
        """ Get the amount of underflow. """
        return self.flow[0]

    @property
    def overflow(self):
        """ Get the amount of overflow. """
        return self.flow[1]

    def integral(self):
        """ Return the integral of the histogram. """
        return np.sum(self.weights*np.diff(self.bins))

    def normalize(self):
        """ Return a histogram normalized to unit area. """
        other = Histogram(copy=self)
        other.weights = other.weights/other.integral()
        return other

    def save(self, name, path=os.getcwd(), mode='text', **kwargs):
        """ Save the histogram either as a text file, or as a figure. """
        if mode == 'text':
            with open(os.path.join(path, name)) as outfile:
                outfile.write(str(self))
        elif mode == 'plot':
            self.weights = np.insert(self.weights, 0, self.weights[0])
            plt.plot(self.bins, self.weights, ds='steps', **kwargs)
            plt.savefig(os.path.join(path, name))
            plt.close()
            self.weights = self.weights[1:]
        else:
            raise ValueError('Expected either mode as `text` or `plot`, '
                             'got mode {}'.format(mode))

    def show(self, **kwargs):
        """ Open window displaying histogram. """
        self.weights = np.insert(self.weights, 0, self.weights[0])
        plt.plot(self.bins, self.weights, ds='steps', **kwargs)
        plt.show()
        self.weights = self.weights[1:]
