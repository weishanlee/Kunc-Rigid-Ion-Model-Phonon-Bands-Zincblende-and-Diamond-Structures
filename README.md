# PDR of Zincblende and Diamond Crystals based on RIM
The codes in the repository intends to reproduce the classic work of phonon dispersion relation (PDR) of zincblende and diamond crystals by K. Kunc and O.H. Nielsen[1],[2].
The main idea of the rigid-ion model(RIM) is that in the Zincblende and Diamond structure, (up to the second neighbour) the interatomic force constants between two atoms are described by the 11 parameters, either obtained by fitting experimental data or from ab initio calculations.

The 11 parameters are given in RIM.DAT, including GaAs, ZnO (in Zincblende), ZnS and AlAs.

Descriptions of the fortran code:

(1)The main program is RIM_MAIN.for.
(2)The RIMCOUL.for calculates the Coulomb interaction with the Ewald sum.
(3)The RIMSR.for calculates the short-range interactions of atoms up to the second neighbour. 
(4)Some regular eispack packages are required, including
	ch.for; 
	epslon.for; 
	htribk.for; 
	htridi.for; 
	pythang.for; 
	tql2.for; 
	tqlrat.for, which may be found via the website: http://www.netlib.org/eispack/index.html
(5)The output phonon dispersion relation is recorded in PDR.DAT, which may be plotted by any other softwares. The plots of dispersion relations of GaAs, ZnO (in Zincblende), ZnS and AlAs are given in the .jpg format.

References:

[1]K. Kunc and O.H. Hielsen, "Lattice Dynamics of Zincblende Structure Compounds Using Deformation-Dipole Model and Rigid Ion Model", Computer PHysics Communications 16 (1979) 181-197.

[2]K. Kunc et at., "Lattice Dynamics of Several $A^{n}B^{8-N}$ Having the Zincblende Structure. I. Deformable-Bond Approximation", Phys. Stat. Sol. (b) 71, 341(1975).
