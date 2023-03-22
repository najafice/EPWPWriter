# EPWPWriter
 Python script for extrapolating stabilized single-point (SSP) brick element stresses to equivalent nodal stress.
> Created by:  **Chris McGann**, **Pedro Arduino**: University of Washington <br/> Re-Written in python by: **Abolfazl Najafi**: Iran Univ. of Science and Tech.

**Mail:** abolfazlbox [at] gmail - **Webpage:** [najafice.github.io](https://najafice.github.io)
## How to use

1. put node data into '*data/nodesInfo.dat*' file with the following format:
	{node id}	{X}	{Y}	{Z}
2. put element data into '*data/elementInfo.dat*' file with the following format:
	{elem id}	{n1}	{n2}	{n3}	{n4}	{n5}	{n6}	{n7}	{n8}
3. write gravity stresses into '*data/Gstress.out*' file, it should have 7 component per element:
	{elem id}	{Sxx}	{Syy}	{Szz}	{Sxy}	{Syz}	{Szx}	{hr}
4. write pore water pressure values into '*data/porePressure.out*' file, it should have 2 column per node, which the first column for the initial pore pressure and the second for the current time step, as below:
	{PP0}	{PPn}
5. Run **EPWPWriter.py** > it will generate '*EPPR.out*' with Nodal Excess Pore Water Pressure.
> **Notes:** This script only works for [SSPBrickUP](https://openseespydoc.readthedocs.io/en/stable/src/SSPbrickUP.html) element in OpenSees which only has 1 integration point.
