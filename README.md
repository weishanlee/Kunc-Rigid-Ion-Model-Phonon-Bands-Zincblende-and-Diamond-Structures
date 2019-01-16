# PDR of Zincblende and Diamond Crystals based on RIM
The codes in the repository intends to reproduce the classic work of phonon dispersion relation (PDR) of zincblende and diamond crystals by K. Kunc and O.H. Nielsen[1],[2].
The main idea of the rigid-ion model(RIM) is that in the Zincblende and Diamond structures, (up to the second neighbour) after considering the symmetry properties of atoms, the interatomic force constants between two atoms are described by the 11 parameters, either obtained by fitting experimental data or from ab initio calculations.

The 11 parameters are given in RIM.DAT, including GaAs, ZnO (in Zincblende), ZnS and AlAs.

Descriptions of the fortran codes:

(1)The main program is RIM_MAIN.for.
(2)The RIMCOUL.for calculates the Coulomb interactions of atoms by the Ewald sum.
(3)The RIMSR.for calculates the short-range interactions up to the second neighbour atoms. 
(4)Some regular eispack packages are required, including
	ch.for; 
	epslon.for; 
	htribk.for; 
	htridi.for; 
	pythang.for; 
	tql2.for; and
	tqlrat.for, which may be found via the website: http://www.netlib.org/eispack/. For copyright reasons the above fortran files are NOT uploaded and are NOT included in this repository.

The output phonon dispersion relation is in PDR.DAT, and may be plotted by any other softwares. The plots of phonon dispersion relations of GaAs, ZnO (in Zincblende), ZnS and AlAs are given in the .jpg format.

The executable file is also included. To calculate other materials, one may try to change the 11 parameters in the RIM.DAT file by just coping other parameters in the lists and pasting to replace the first two lines, and run the exe file. The codes are designed to read the RIM.DAT file only for the first two lines.

References:

[1]K. Kunc and O.H. Hielsen, "Lattice Dynamics of Zincblende Structure Compounds Using Deformation-Dipole Model and Rigid Ion Model", Computer PHysics Communications 16 (1979) 181-197.

[2]K. Kunc et at., "Lattice Dynamics of Several $A^{N}B^{8-N}$ Having the Zincblende Structure. I. Deformable-Bond Approximation", Phys. Stat. Sol. (b) 71, 341(1975).
