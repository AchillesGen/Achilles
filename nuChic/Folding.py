from scipy import interpolate
from absl import flags
from absl import logging
import numpy as np
import os

FLAGS = flags.FLAGS

DIR, FILE = os.path.split(__file__)

class Folding:
    def __init__(self):
        logging.info('Initializing folding function')
        with open(os.path.join(DIR,'pke','folding.in')) as f:
            TA, hwfold, nwfold = f.readline().split()
            TA = float(TA)
            hwfold = float(hwfold)
            nwfold = int(nwfold)
            self.wfold = np.empty(nwfold)
            fold = np.empty(nwfold)
            for i in range(nwfold):
                self.wfold[i], fold[i] = f.readline().split()
                logging.debug('wfold[%d] = %e, fold[%d] = %e' 
                    % (i, self.wfold[i], i, fold[i]))

        logging.info('TA = %e' % TA)
        self.TA = 1.0 - 2.0*np.sum(fold)*hwfold
        logging.info('Norm folding = %e' % TA)

        self.fold = interpolate.interp1d(self.wfold,fold) 

        logging.info('Initializing optical potential')
        with open(os.path.join(DIR,'pke','realOP_12C_EDAI.dat')) as f:
            npot = int(f.readline())
            self.kin = np.empty(npot)
            pot = np.empty(npot)
            for i in range(npot):
                self.kin[i], pot[i] = f.readline().split()
                logging.debug('kin[%d] = %e, potential[%d] = %e' 
                    % (i, self.kin[i], i, pot[i]))
            
        self.kinematic = interpolate.interp1d(self.kin, pot)

    def __call__(self,w,wp):
        if abs(w-wp) <= self.wmax and abs(w-wp) >= self.wmin:
            return self.fold(abs(w-wp))
        else:
            return 0

    @property
    def wmin(self):
        return self.wfold[0]

    @property
    def wmax(self):
        return self.wfold[-1]

    @property
    def kin_min(self):
        return self.kin[0]

    @property
    def kin_max(self):
        return self.kin[-1]
