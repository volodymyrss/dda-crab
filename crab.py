from numpy import *

class Crab:
    def __init__(self):
         self.Crab=load("/Integral/throng/common/energy_calibration/crab_spectrum/example/crab_example.npy")
         self.Crab_ecs=genfromtxt("/Integral/throng/common/energy_calibration/crab_spectrum/example/crab_unfolded_example.txt").transpose()

    def counts_in(self,e1,e2):
        _e1,_e2,_c=self.Crab
        return _c[logical_and(_e1>e1,_e2<e2)].sum()

    def to_mcrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000

    def to_mucrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000*1000

    def to_ecs(self,c,e1,e2):
        _ec,_de,_f,_fe=self.Crab_ecs
        _e1=_ec-_de
        _e2=_ec+_de

        s=logical_and(_e1>e1,_e2<e2)
        
        crab_ecs=(_f*_de*2)[s].sum()/(e2-e1)
        return  self.to_mcrab(c,e1,e2)/1000*crab_ecs


