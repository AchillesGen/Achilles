#import lhapdf
#
#class PDF:
#    def __init__(self,name,iset):
#        self.name = name
#        self.iset = iset
#        self.pdf = lhapdf.mkPDF(self.name,self.iset)
#
#    def fxQ(self,pid,x,Q):
#        return self.xfxQ2(pid,x,Q*Q)/x
#
#    def fxQ2(self,pid,x,Q2):
#        return self.xfxQ2(pid,x,Q2)/x
#
#    def xfxQ2(self,pid,x,Q2):
#        return self.pdf.xfxQ2(pid,x,Q2)
#
#    def xfxQ(self,pid,x,Q):
#        return self.xfxQ2(pid,x,Q*Q)
#
#    def alphasQ2(self,Q2):
#        return self.pdf.alphasQ2(Q2)
#
#    def alphasQ(self,Q):
#        return self.alphasQ2(Q*Q)
#
#if __name__ == '__main__':
#    import numpy as np
#    pdf = PDF('CT10nlo',0)
#
#    xs = [x for x in np.logspace(-7,0,5)]
#    qs = [q for q in np.logspace(1,4,4)]
#    gluon_xfs = np.empty([len(xs), len(qs)])
#    for ix, x in enumerate(xs):
#        for iq, q in enumerate(qs):
#            gluon_xfs[ix,iq] = pdf.xfxQ(21,x,q)
#    print(gluon_xfs)
