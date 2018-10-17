# Program Layout

## Beams:
 - The nuetrino energy is given by ![equation](http://latex.codecogs.com/gif.latex?E_%5Cnu%20%3D%20%5Cfrac%7Bm_%5Cpi%5E2-m_%5Cmu%5E2%7D%7B2%28E_%5Cpi-p_%5Cpi%5Ccos%5Ctheta_%5Cnu%29%7D)
 - The target can be defined by the number of protons and the number of neutrons
 - Is there a way to make this more general?

## particle class:
 - PID of the particle
 - Four momentum of the particle
 - etc.

## NucModel: 
 - Should contain some selector based on the input Q value to determine if it should use QE, DIS, etc. 
 - Maybe calculate all and chose based off the ratio of a given process to the total

## Pseudocode

```python
b1 = Beam(name='nu',Epi=Epi,Ppi=Ppi,cosTheta=cosTheta)

# The target beam is given by the number of protons and the number of neutrons only
b2 = Beam(name='nucleus',nProtons=nProtons,nNeutrons=nNeutrons)

# Generate events (Calls MakeEvent through an integrator and a phase space generator)
GenerateEvents(b1,b2)

# DectectorSim if needed
if dectectorSim:
    event = ToGeant4(wgt,event)
 
# Write event to file for analysis
ToFile(wgt,event)

def MakeEvent(b1,b2):
    # Integrate over the following:
    # Need to generate the incoming and outgoing momenta for the leptonic sector; and type of interaction (W,\gamma)
    event_leptonic = GenerateLeptonic(b1)
    Q = -pow(event_leptonic[0]+event_leptonic[1],2)
    interaction = 'W' if event_leptonic[1].pid is lepton else 'A'

    wgt, particle_hadronic = NucModel(b2,Q,interaction)

    event_hadronic = Cascade(particle_hadronic)  

    if event_leptonic[1] is lepton:
        event_leptonic = PhotonShower(event_leptonic)

    event = event_leptonic + event_hadronic

    # Selection criteria and cuts??
    cut = Cuts(event)

    if cut:
        return 0

    return wgt #to integrator
```

# Remaining Questions:
1. How to determine where in detector event occurs?
2. Should the neutrino beam be defined through the pion kinematics that produced it, or should we generate a profile in some other way?
