from numpy import *
import os

class Crab:
    def __init__(self):
         self.Crab=load(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_example.npy")
         self.Crab_ecs=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_unfolded_example.txt").transpose()
         
         self.Crab_jemx=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_jemx_counts.qdp")
         self.Crab_jemx_ecs=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_jemx_uf.qdp").transpose()
    
    def counts_in(self,e1,e2,kind="ISGRI"):
        if kind=="ISGRI":
            _e1,_e2,_c=self.Crab
        else:
            _ec,_de,_c,_ce,_m=self.Crab_jemx.transpose()
            _e1=_ec-_de
            _e2=_ec+_de
            _c*=_de
        return _c[logical_and(_e1>e1,_e2<e2)].sum()

    def to_mcrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000

    def to_mucrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000*1000
    
    def ecs_in(self,e1,e2,kind="ISGRI"):
        if kind=="ISGRI":
            _ec,_de,_f,_fe=self.Crab_ecs
            _e1=_ec-_de
            _e2=_ec+_de
            _f*=1/_ec
        else:
            _ec,_de,_f,_fe,_m=self.Crab_jemx_ecs
            _e1=_ec-_de
            _e2=_ec+_de
            _f*=1.6e-9/_ec


        s=logical_and(_e1>e1,_e2<e2)
        
        return (_f*_de*2)[s].sum() #/sum(_de[s])/2 #(e2-e1)

    def to_ecs(self,e1,e2):
        _ec,_de,_f,_fe=self.Crab_ecs
        _e1=_ec-_de
        _e2=_ec+_de

        s=logical_and(_e1>e1,_e2<e2)
        
        crab_ecs=(_f*_de*2)[s].sum()/(e2-e1)
        return  self.to_mcrab(c,e1,e2)/1000*crab_ecs


    def to_ecs(self,c,e1,e2,kind="ISGRI"):
        return c/self.counts_in(e1,e2,kind)*self.ecs_in(e1,e2,kind)
