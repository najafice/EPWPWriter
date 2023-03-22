# EPWPWriter
Python script for extrapolating stabilized single-point (SSP) brick element stresses to equivalent nodal stress.
Wrriten by Abolfazl Najafi, M.Sc. - Copyright 2023, All right reserved.
Mail: abolfazlbox@gmail.com - Webpage: najafice.github.io
-----------------------------------------------------------------------------------------------------
Python script: EPWPWriter, v1.0
-----------------------------------------------------------------------------------------------------
1. put node data into 'data/nodesInfo.dat' file with the following format:
	{node id}	{X}	{Y}	{Z}

2. put element data into 'data/elementInfo.dat' file with the following format:
	{elem id}	{n1}	{n2}	{n3}	{n4}	{n5}	{n6}	{n7}	{n8}

3. write gravity stresses into 'data/Gstress.out' file, it should have 7 component per element:
	{elem id}	{Sxx}	{Syy}	{Szz}	{Sxy}	{Syz}	{Szx}	{hr}

4. write pore water pressure values into 'data/porePressure.out' file, it should have 2 column per node, 
   which the first column for the initial pore pressure and the second for the current time step, as below:
	{PP0}	{PPn}

5. Run EPWPWriter.py > it will generate 'EPPR.out' with Nodal Excess Pore Water Pressure.
-----------------------------------------------------------------------------------------------------
Note: This script only works for SSPBrickUP element in OpenSees which only has 1 integration point.
Hexahedron element nodal order:
       v
6----------7
|\     ^   |\
| \    |   | \
|  \   |   |  \
|   5------+---8
|   |  +-- |-- | -> u
1---+---\--2   |
 \  |    \  \  |
  \ |     \  \ |
   \|      w  \|
    3----------4
