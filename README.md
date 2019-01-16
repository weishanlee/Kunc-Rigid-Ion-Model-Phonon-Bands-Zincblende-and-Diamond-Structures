# Phonon Bands of Zincblende and Diamond Crystals based on RIM
The codes in the repository intends to reproduce the classic work of phonon dispersion relations (PDR, or phonon bands) of zincblende and diamond crystals by K. Kunc et al[1],[2].
The main idea of the rigid-ion model(RIM) is that in the Zincblende and Diamond structures, (up to the second neighbour) after considering the symmetry properties of atoms, the interatomic force constants between two atoms are described by the 10 parameters (A, B, C1, D1, E1, F1, C2, D2, E2, F2,  either obtained by fitting experimental data or from ab initio calculations), together with Z2, the effective charge of two atoms (zero for diamond structure), as well as A (here representing for the lattice constant), and M1 and M2 (the atomic weights of the two atoms in a unit cell). The 14 parameters are given in RIM.DAT, including GaAs, ZnO (in Zincblende), ZnS, AlAs, and C (in diamond), Si, Ge, and Sn. 

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

Notice that parameters in E1 and E2 are generally NOT zero even for Group-IV elements. However, it is not easy to find out the appropriate values for E1 and E2. In Ref [2], Kunc et al. indicated that for central force assumptions we may set E1 and E2 to be zero.  

The executable file is also included. To calculate other materials, one may try to change the 11 parameters in the RIM.DAT file by just coping other parameters in the lists and pasting to replace the original parameters. Later, run the exe file. The codes are designed to read the RIM.DAT file only for the 11 parameters.

References:

[1]K. Kunc and O.H. Hielsen, "Lattice Dynamics of Zincblende Structure Compounds Using Deformation-Dipole Model and Rigid Ion Model", Computer PHysics Communications 16 (1979) 181-197.

[2]K. Kunc et at., "Lattice Dynamics of Several $A^{N}B^{8-N}$ Having the Zincblende Structure. I. Deformable-Bond Approximation", Phys. Stat. Sol. (b) 71, 341(1975).
