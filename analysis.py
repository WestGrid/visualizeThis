# Edited version of WestGrid's 2018 Visualize This! submission by Usman Alim and Roberta Cabral Ramos
# Mota (Dept. of Computer Science, University of Calgary).

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.density import density_from_Universe
from MDAnalysis.analysis.distances import distance_array
import tvtk.api as tvtk

base = '/Users/razoumov/Documents/visualizeThis/2018/data/scientific/'
u = mda.Universe(base+'sim.part0001_protein_lipid_popc.gro', [base+'5us_traj_popc_protein_lipid_popc.xtc',base+'15us_traj_protein_lipid_popc.xtc'])

# define the contact condition in ångströms

contact_dist = 5
nearby = "name PO4 and (around " + str(contact_dist) + " name SC1 or around " + str(contact_dist) + " name SC2 or around " + str(contact_dist) + " name SC3 or around " + str(contact_dist) + " name SC4)"

# write out 3D density of PO4 everywhere and near the protein in VTK Image Data

Dall = density_from_Universe(u, delta=0.5, atomselection="name PO4", update_selection=True)
Dnear = density_from_Universe(u, delta=1, atomselection=nearby, update_selection=True)

print('near:', Dnear.grid.shape)   # should be (60, 57, 51)
grid = tvtk.tvtk.ImageData(spacing=Dnear.delta, origin=Dnear.origin, dimensions=Dnear.grid.shape)
grid.point_data.scalars = Dnear.grid.ravel(order='F')
grid.point_data.scalars.name = 'POPC near SC1_4'
tvtk.write_data(grid, 'Dnear.vtk')   # omit `vtk` extension to write newer XML VTK, otherwise legacy VTK

print('all:', Dall.grid.shape)   # should be (302, 298, 124)
grid = tvtk.tvtk.ImageData(spacing=Dall.delta, origin=Dall.origin, dimensions=Dall.grid.shape)
grid.point_data.scalars = Dall.grid.ravel(order='F')
grid.point_data.scalars.name = 'POPC all'
tvtk.write_data(grid, 'Dall.vtk')   # omit `vtk` extension to write newer XML VTK, otherwise legacy VTK

### identify membrane's PO4 beads that interact with the protein, and classify these into high-, medium-,
### and low-interacting

selection = u.select_atoms(nearby, updating=True)
clist = []   # will contain a list of 1102 frames, each frame listing PO4 beads (IDs) close to protein's SC{1..4}
for timestep in u.trajectory:
    clist.append(set(selection.ids))

itrajs = dict();   # a dictionary of PO4 beads (IDs) that stay close to protein's SC{1..4} for 5, 6, ... frames
for wsize in range(5,25):
    bindList = []
    for i in range(0,len(clist) - wsize):
        bindList.append(set.intersection(*clist[i:i+wsize]))
    itrajs[wsize] = set.union(*bindList)

# partition the trajectories into groups
A = set.intersection(itrajs[5],  itrajs[6],  itrajs[7])
B = set.intersection(itrajs[8],  itrajs[9],  itrajs[10])
C = set.intersection(itrajs[11], itrajs[12], itrajs[13])

LOW = A - B - C
MID = B - C
HIGH = C

print("Number of HIGH PO4 beads: " + str(len(HIGH)))
print("Number of MEDIUM PO4 beads: " + str(len(MID)))
print("Number of LOW PO4 beads: " + str(len(LOW)))

# for each MEDIUM PO4 bead, write a VTK file with all its positions and the shortest distance to a
# protein's SC{1..4} bead

prefix = "mid"           # change to 'high' for HIGH beads
for id in sorted(MID):   # change to HIGH for HIGH beads
    atom_selection = u.select_atoms("bynum " + str(id), updating=True)
    protein_selection = u.select_atoms("name SC1 or name SC2 or name SC3 or name SC4", updating=True)
    points = tvtk.tvtk.Points()
    distances = tvtk.tvtk.FloatArray()
    distances.name = 'distances'
    times = tvtk.tvtk.FloatArray()
    times.name = 'TIME'
    for i, timestep in enumerate(u.trajectory):
        d = distance_array(atom_selection.positions, protein_selection.positions)
        if d.min() <= contact_dist:
            points.insert_next_point(tuple(timestep.positions[id]))
            distances.insert_next_value( d.min() )
            times.insert_next_value( i )

    polyline = tvtk.tvtk.PolyData()
    polyline.trait_set(points=points)
    polyline.point_data.scalars = distances
    polyline.point_data.add_array(times)
    tvtk.write_data(polyline, prefix + str(id) + ".vtk")

# for each of 1102 frames, write a VTK file `high/HIGH?????.vtk` that contains the positions of all HIGH
# PO4 beads at that timestep

for i, timestep in enumerate(u.trajectory):
    points = tvtk.Points()
    polyData = tvtk.PolyData()
    for id in HIGH:
        atom_selection = u.select_atoms("bynum " + str(id), updating=True)
        points.insert_next_point(tuple(timestep.positions[id]))
        
    polyData.trait_set(points=points)
    write_data(polyData, "high/HIGH" + "{:0>5d}".format(i) + ".vtk")
