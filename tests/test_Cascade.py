from nuChic.Cascade import FSI
from nuChic.Nucleus import Nucleus
from nuChic.Constants import MeV, fm, mqe as mN
from nuChic.FourVector import Vec4 
from nuChic.ThreeVector import Vec3
from nuChic.Particle import Particle
import numpy as np
from numpy import linalg as LA
import pytest
from mock import patch

nuc = Nucleus(6,12,10,225)
dt=1*fm
cascade_b = FSI(nuc, dt)

nprotons = 6
protons = np.random.random((nprotons,3))
nneutrons = 6
neutrons = np.random.random((nneutrons,3))

@patch('nuChic.Cascade.Nucleus')
def test_Cascade_init(MockNucleus):
    MockNucleus.generate_config.return_value = (protons,neutrons) 
    cascade = FSI(MockNucleus, dt)
    assert cascade.number_nucleons == nprotons + nneutrons
    assert cascade.number_protons == nprotons 
    assert cascade.number_neutrons == nneutrons
    assert MockNucleus.generate_config.called_once()

@patch('nuChic.Cascade.Nucleus')
def test_kick(MockNucleus):
    MockNucleus.generate_config.return_value = (protons,neutrons) 
    cascade = FSI(MockNucleus, dt)

    Ep = mN + 500*MeV
    pp = np.sqrt(Ep**2-mN**2)
    energy_transfer=Vec4(Ep,0,0,pp)
    cascade.kick(energy_transfer)

    assert len(cascade.kicked_idxs)==1
    assert cascade.nucleons[cascade.kicked_idxs[0]].status == -1
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.E == Ep
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.px == 0
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.py == 0
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.pz == pp

@patch('nuChic.Cascade.Nucleus')
def test_reset(MockNucleus):
    MockNucleus.generate_config.return_value = (protons,neutrons) 
    cascade = FSI(MockNucleus, dt)

    Ep = mN + 500*MeV
    pp = np.sqrt(Ep**2-mN**2)
    energy_transfer=Vec4(Ep,0,0,pp)
    cascade.kick(energy_transfer)
    cascade.reset()

    assert cascade.number_nucleons == nprotons + nneutrons
    assert cascade.number_protons == nprotons 
    assert cascade.number_neutrons == nneutrons
    assert cascade.cylinder_pt1 ==0
    assert cascade.cylinder_pt2 ==0
    assert set([abs(cascade.nucleons[i].mom.E-mN)<1E-10 for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.px) for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.py) for i in range(len(cascade.nucleons))])
    assert set([np.isnan(cascade.nucleons[i].mom.pz) for i in range(len(cascade.nucleons))])

def test_points_in_cylinder():    
    assert FSI.points_in_cylinder(pt1=[0,0,0], pt2=[0,0,1],
                                      r=1, q=[0.2,0.2,0.2])
    assert not(FSI.points_in_cylinder(pt1=[0,0,0], pt2=[0,0,1],
                                          r=1, q=[0,1.1,0.5]))
    
def test_to_cartesian():
    coords = [[1,0,0], [1,np.pi/2,0], [1,np.pi/2,np.pi/2]]
    results = [[0,0,1], [0,1,0], [1,0,0]]
    for i in range(3):
        assert LA.norm((FSI.to_cartesian(coords[i])-results[i])) < 1E-10
        
@patch('nuChic.Cascade.Nucleus')
def test_pauli_blocking(MockNucleus):
    kf = 225 * MeV
    MockNucleus.generate_config.return_value = (protons,neutrons) 
    MockNucleus.kf = kf
    cascade = FSI(MockNucleus, dt)

    fv1 = Vec4(0,0,0,0.99*kf)
    fv2 = Vec4(0,0,0,1.01*kf)
    assert cascade.pauli_blocking(fv1)
    assert not cascade.pauli_blocking(fv2)
        
@patch('nuChic.Cascade.Nucleus')
def test_interacted(MockNucleus):
    MockNucleus.generate_config.return_value = (protons,neutrons) 
    cascade = FSI(MockNucleus, dt)

    Ep = 500 * MeV # massless proton for easy propagation
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
    
    
    
    
    
    
    
    
    
    
    
    
