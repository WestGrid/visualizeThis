---
layout: page
title: Scientific Dataset Details
---

**You will be able to download the dataset on October 1, 2018.**

<!-- Click here to download the dataset as a 293MB compressed ZIP file. If you have any issues accessing the -->
<!-- dataset using that link, you can try this alternative link or email alex.razoumov@westgrid.ca for -->
<!-- assistance. -->

<!-- The dataset is taken from a simulation of the air flow around multiple counter-rotating vertical-axis -->
<!-- wind turbines. It represents a single-time snapshot and contains the following five scalar variables on a -->
<!-- 3D unstructured mesh: -->

<!-- * pressure p relative to the reference pressure for normal conditions (can be positive or negative), -->
<!-- * Q-criterion for vorticity q, -->
<!-- * three Cartesian velocity components u, v, w. -->

<!-- The original data produced by the simulation code was stored in Tecplot file format with explicit x,y,z -->
<!-- coordinates for each of 5,307,199 mesh points and a list of 8 vertices for each of 16,679,253 cells. To -->
<!-- compress the dataset and speed up I/O, we converted the data into a VTK file format, with the vertex -->
<!-- coordinates and connections defining the cells encoded in the mesh itself. If for some reason you want -->
<!-- data in the original Tecplot format, please let us know. -->

<!-- The main VTK file air.vtu stores field data (the five scalar variables) on a volumetric unstructured -->
<!-- grid, whereas 10 smaller files blade{21,22,23,24,25,26,43,44,45,46}.vtp contain the field data on -->
<!-- polygonal surface meshes defining the turbine blades. -->

<!-- The main challenge in visualizing this dataset is dealing with the range of scales from the entire -->
<!-- simulation volume to small-scale turbulence near the blades. -->

## How to read the data

This information will be posted shortly.

<!-- Both ParaView and VisIt can open standard VTK unstructured grid and polygonal data files. These files can -->
<!-- also be read in Python and C++ using the VTK library (http://www.vtk.org). Note that the volumetric file -->
<!-- air.vtu contains over 16 million cells so it might take some time to read it, depending on your -->
<!-- computer's speed. -->

<!-- As a result, rendering might also take some time, and can be automated with scripting, so that you could -->
<!-- leave a script running for a few hours and come back to let's say several hundred frames of a movie, or -->
<!-- it could be run in parallel on a cluster. -->

<!-- We hosted a kickoff webinar on Sep-27 that gave a walk-through of the dataset. Click here to view the -->
<!-- archive recording (to get right to the dataset tour, skip ahead to 4:04 in the video). -->

## Sample Visualizations

We will provide sample visualizations shortly.

<!-- To give you an idea of the type of data in these files and to help you with actual visualizations, we -->
<!-- provide two sample visualizations, one done with ParaView and the other one with VisIt. Both workflows -->
<!-- demonstrate loading of all 11 VTK files. -->

<!-- The ParaView state file bladesWithLines.pvsm stores the pipeline to visualize the blades (coloured by the -->
<!-- pressure on their surfaces) and the airflow around them with uniform-colour streamlines. You can point -->
<!-- ParaView to this state file with File - Load State..., or start ParaView from the command line with -->
<!-- "paraview --state=bladesWithLines.pvsm". The resulting image bladesWithLines.png is shown below. -->

<!-- /files/webfm/Communications/bladesWithLines.png -->

<!-- The VisIt Python script positiveNegativePressure.py renders semi-transparent isosurfaces of positive -->
<!-- (blue) and negative (turquoise) pressure around the blades. You can run this script in VisIt either from -->
<!-- Controls - Launch CLI... or from Controls - Command..., or from the command line with "visit -nowin -cli -->
<!-- -s positiveNegativePressure.py". The resulting image positiveNegativePressure0000.png is shown below. -->

<!-- /files/webfm/Communications/positiveNegativePressure0000.png -->

<!-- We are looking for innovative visualizations of this dataset. For example, one could enhance these -->
<!-- renderings by drawing streamlines around the isosurfaces and producing some animations such as spinning -->
<!-- the visualization around the vertical axis or gradually turning on/off various visualization -->
<!-- elements. Speaking more generally, a nice animation would help us explore the spatial range and values of -->
<!-- multiple variables and show how various elements of the simulation are tied together. -->
