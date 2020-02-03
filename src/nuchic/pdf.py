""" Implement PDF interface class. """

# from absl import logging

try:
    import lhapdf
except ImportError:
    # logging.warn('Could not find lhapdf. Using dummy class.')
    lhapdf = None


class PDF:
    """ Class for accessing information about the PDF. """
    def __init__(self, name, iset):
        self.name = name
        self.iset = iset

        if lhapdf is None:
            self.pdf = None
        else:
            self.pdf = lhapdf.mkPDF(self.name, self.iset)

    def fxq(self, pid, mom_frac, scale):
        """ Calculate value of pdf for a given momentum faction and scale. """
        if self.pdf is None:
            return mom_frac
        return self.xfxq2(pid, mom_frac, scale**2)/mom_frac

    def fxq2(self, pid, mom_frac, scale2):
        """ Calculate value of pdf for a given momentum faction and scale. """
        return self.xfxq2(pid, mom_frac, scale2)/mom_frac

    def xfxq2(self, pid, mom_frac, scale2):
        """ Calculate value of xpdf for a given momentum faction and scale. """
        return self.pdf.xfxQ2(pid, mom_frac, scale2)

    def xfxq(self, pid, mom_frac, scale):
        """ Calculate value of xpdf for a given momentum faction and scale. """
        return self.xfxq2(pid, mom_frac, scale**2)

    def alphasq2(self, scale2):
        """ Calculate the value of alpha_s at a given scale. """
        return self.pdf.alphasQ2(scale2)

    def alphasq(self, scale):
        """ Calculate the value of alpha_s at a given scale. """
        return self.alphasq2(scale**2)
