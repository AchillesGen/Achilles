from .nuChic.Cascade import FSI
from .nuChic.Nucleus import Nucleus
from .nuChic.Constants import MeV, fm, mqe as mN
from .nuChic.FourVector import Vec4 
from .nuChic.ThreeVector import Vec3
from .nuChic.Particle import Particle
import numpy as np
from numpy import linalg as LA
import pytest

nuc = Nucleus(6,12,10,225)
energy_transfer=500*MeV
dt=1*fm
cascade = FSI(nuc, dt)
cascade_b = FSI(nuc, dt)

def test_Cascade_init():
#    assert len(cascade.outgoing_particles) == 0
    number_of_protons = 0
    number_of_neutrons = 0
    for i in range(len(cascade.nucleons)) :
        if cascade.nucleons[i].pid == 2212 :
            number_of_protons += 1
        if cascade.nucleons[i].pid == 2212 :
            number_of_neutrons += 1
    assert number_of_protons == 6 and number_of_neutrons == 6

def test_kick():
    cascade.kick(energy_transfer)
    assert len(cascade.kicked_idxs)==1
    assert abs(cascade.nucleons[cascade.kicked_idxs[0]].mom.E - (mN+energy_transfer))<1E-10
    assert not(np.isnan(cascade.nucleons[cascade.kicked_idxs[0]].mom.px))
    assert not(np.isnan(cascade.nucleons[cascade.kicked_idxs[0]].mom.py))
    assert not(np.isnan(cascade.nucleons[cascade.kicked_idxs[0]].mom.pz))

def test_reset():
    cascade.kick(energy_transfer)
    cascade.reset()
    assert cascade.cylinder_pt1 ==0
    assert cascade.cylinder_pt2 ==0
    assert set([abs(cascade.nucleons[i].mom.E-mN)<1E-10 for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.px) for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.py) for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.pz) for i in range(len(cascade.nucleons))])

def test_points_in_cylinder():    
    assert cascade.points_in_cylinder(pt1=[0,0,0], pt2=[0,0,1], r=1, q=[0.2,0.2,0.2])
    assert not(cascade.points_in_cylinder(pt1=[0,0,0], pt2=[0,0,1], r=1, q=[0,1.1,0.5]))
    
def test_to_cartesian():
    coords = [[1,0,0], [1,np.pi/2,0], [1,np.pi/2,np.pi/2]]
    results = [[0,0,1], [0,1,0], [1,0,0]]
    for i in range(3):
        assert LA.norm((cascade.to_cartesian(coords[i])-results[i])) < 1E-10
        
def test_pauli_blocking():
    kf = cascade.nucleus.kf
    fv1 = Vec4(0,0,0,0.99*kf)
    fv2 = Vec4(0,0,0,1.01*kf)
    assert cascade.pauli_blocking(fv1)
    assert not(cascade.pauli_blocking(fv2))
        
def test_interacted():
    Ep = energy_transfer # massless proton for easy propagation
    pp = np.sqrt(Ep**2)

    positions = [Vec3(100,100,100), Vec3(0,0,dt/2), Vec3(0,0,1.5*dt), Vec3(0,0.5,0.5*dt), Vec3(0,1.5,0.5*dt)]
    tests = [False, True, False, True, False]
    for i in range(len(positions)) :
        # Put all other nucleaons everyone far
        cascade_b.nucleons = [Particle(pid=2212, mom=Vec4(mN,0,0,0), pos=Vec3(100,100,100)) for i in range(len(cascade.nucleons))]
        # Kicked particle
        cascade_b.nucleons[0] = Particle(pid=2212, mom=Vec4(Ep,0,0,pp), pos=Vec3(0,0,0))
        cascade_b.kicked_idxs = []
        cascade_b.kicked_idxs.append(0)
        cascade_b.nucleons[cascade_b.kicked_idxs[0]].status = -1 # propagating nucleon
        cascade_b.nucleons[10] = Particle(pid=2212, mom=Vec4(mN,0,0,0), pos=positions[i])
        print(cascade_b.nucleons[0].pos, cascade_b.nucleons[10].pos)
        in_cylinder, idx = cascade_b.interacted(0, sigma=1)
#        assert in_cylinder==tests[i]
        print(in_cylinder, tests[i], cascade_b.nucleons[0].pos, cascade_b.nucleons[10].pos)
    
    
    
    
    
    
    
    
    
    
    
    
