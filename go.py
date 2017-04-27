from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys
#from mdtraj.reporters import NetCDFReporter, HDF5Reporter

pdb = PDBFile('input.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('./MDFTFF.xml','./spce.xml')

print('Solvating by SPCE...')
modeller.addSolvent(forcefield, model='spce', padding=1.5*nanometers, neutralize=False)

print('Configuring simulation...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer)
#integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.01)
integrator.setRandomNumberSeed(0)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('solvatedBox.pdb', 'w'))

print('Minimizing')
simulation.minimizeEnergy()

print('Printing minimizedSolvatedBox')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimizedSolvatedBox.pdb', 'w'))

printevery = 50

#simulation.reporters.append(PDBReporter('output.pdb', printevery))
simulation.reporters.append(DCDReporter('output.dcd', printevery))
steps = 3500000 # 500000 == 1ns

simulation.reporters.append(StateDataReporter(sys.stdout, printevery, progress=True, remainingTime=True, step=True, time=True, totalSteps=steps, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True, density=True, speed=True, separator='  '))

simulation.reporters.append(StateDataReporter('state.out', printevery, progress=True, remainingTime=True, step=True, time=True, totalSteps=steps, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True, density=True, speed=True, separator='  '))

simulation.reporters.append(CheckpointReporter('checkpnt.chk', 50))

if os.path.isfile('checkpnt.chk'):
    print('Reading checkpoint...')
    simulation.loadCheckpoint('checkpnt.chk')

simulation.step(steps)
