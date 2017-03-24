from numpy import *
import os

class Crab:
    def __init__(self,kind="ISGRI"):
        self.kind=kind

        if kind=="ISGRI":
            self.Crab=load(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_example.npy")
            self.Crab_ecs=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_unfolded_example.txt").transpose()

            self.e1, self.e2, self.c=self.Crab
            self.ecs_ec, self.ecs_de, self.ecs_f, self.ecs_fe=self.Crab_ecs
            self.ecs_e1=self.ecs_ec-self.ecs_de
            self.ecs_e2=self.ecs_ec+self.ecs_de

        elif kind=="JEMX":
            self.Crab_jemx=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_jemx_counts.qdp")
            self.Crab_jemx_ecs=genfromtxt(os.environ['HOME']+"/soft/energy_calibration/crab_spectrum/example/crab_jemx_uf.qdp").transpose()

            self.e1, self.e2, self.c, self.ce, self.m=self.Crab_jemx.transpose()
            self.ecs_ec, self.ecs_de, self.ecs_f, self.ecs_fe,self.ecs_m=self.Crab_jemx_ecs
            self.ecs_f*=1.6e-9
            self.ecs_fe*=1.6e-9
            self.ecs_e1=self.ecs_ec-self.ecs_de
            self.ecs_e2=self.ecs_ec+self.ecs_de
        else:
            raise Exception("unknown Crab kind!")
    
    def counts_in(self,e1,e2):
        _c=self.c*self.de
        return _c[logical_and(self.e1>e1,self.e2<e2)].sum()

    def to_mcrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000

    def to_mucrab(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*1000*1000
    
    def ecs_in(self,e1,e2,kind="ISGRI"):
        _f=self.ecs_f/self.ecs_ec

        s=logical_and(self.ecs_e1>e1,self.ecs_e2<e2)
        
        return (_f*self.ecs_de*2)[s].sum() #/sum(_de[s])/2 #(e2-e1)

    def to_ecs(self,e1,e2):
        _ec,_de,_f,_fe=self.Crab_ecs
        _e1=_ec-_de
        _e2=_ec+_de

        s=logical_and(_e1>e1,_e2<e2)
        
        crab_ecs=(_f*_de*2)[s].sum()/(e2-e1)
        return  self.to_mcrab(c,e1,e2)/1000*crab_ecs


    def to_ecs(self,c,e1,e2):
        return c/self.counts_in(e1,e2)*self.ecs_in(e1,e2)
