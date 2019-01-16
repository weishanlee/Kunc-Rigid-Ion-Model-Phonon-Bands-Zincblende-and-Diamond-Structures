	  PROGRAM Rigid_Ion_Model
	  implicit none
	  integer, parameter :: ndsr=14, nk=801, na=2,nb=3*na
	  ! ndsr : number of input parameters
          ! nk   : number of points chosen along k direction
	  !	   in the plot of dispersion curve
	  ! na   : number of atoms in a unit cell
	  ! nb   : number of coordinates
	  INTEGER IQ , J !, k , i
	  real*8 FBD(nb,nk),pshr(ndsr)   !FBD: returned values
	                                 !pshr: input parameters
	  real*8 PI
          OPEN(UNIT=62,FILE='RIM.DAT',STATUS='old')
	  OPEN(UNIT=20,FILE='PDR.DAT',STATUS='REPLACE') ! output
	  OPEN(UNIT=25,FILE='NKL_801.DAT',STATUS='old')
	  read(62,*) pshr ! read RIM.DAT input file for the 11 parameters.
          PI = dacos(-1.d0)
          ! Note: even for group IV, E1 and E2 are NOT zero!
          call RIMPDR(PSHR,FBD)
	  DO IQ=1,nk
		WRITE (20,8002) float( iQ-1 ),(FBD(J,IQ),J=1,nb)
          END DO
 8002     FORMAT(8F10.5)
	  END PROGRAM Rigid_Ion_Model
