C         SUBROUTINE FOR
C         COULOMB MATRIX ( IN UNIT OF e**2/Va ) by Ewald sum
C         IN GENERAL DIRECTIONS
C	  (PHASE USED: E(I*K*X(LK))
          SUBROUTINE RIMCOL(JN,JH,QX,QY,QZ,RI,CC,CC1,A)
          IMPLICIT real*8 (a-h,o-z)
          INTEGER , PARAMETER :: na=2
          DIMENSION T(3),VQ(3),Q(3),TQ(3),RI(na,3),PL(3,3),PRL(3,3)
          REAL*8 RL(3),LAT(3),R(3,na),R1(3),X(3),CC1(na,3,3)
          REAL*8 erf
          COMPLEX*16 CI,PHI1,PHI2,HAB,CC(na,na,3,3) !,SUM,PHI
C
          CI = (0.,1.)
	  PI=dacos(-1.d0)
	  C=A
	  VA=A**3/4.
	  B = 2*PI/A
	  BC= 2*PI/C
C---  ETA=ALPHA**2 ---------
      ETA=1.*SQRT(B*BC)
C---  fcc PRIMITIVE LATTICE VECTORS  (in a/2)
	  PL(1,1)=0.   ! (xyz, # atoms)
	  PL(2,1)=1.
	  PL(3,1)=1.

	  PL(1,2)=1.
	  PL(2,2)=0.
	  PL(3,2)=1.

	  PL(1,3)=1.
	  PL(2,3)=1.
	  PL(3,3)=0.

C---  RECIPROCAL L.V.  (in 2pi/a)
	  PRL(1,1)=-1.
	  PRL(2,1)=1.
	  PRL(3,1)=1.

	  PRL(1,2)=1.
	  PRL(2,2)=-1.
	  PRL(3,2)=1.

	  PRL(1,3)=1.
	  PRL(2,3)=1.
	  PRL(3,3)=-1.
C	
C     GENERATE POSITIONS
C
	  DO 5 II=1,na
	  DO 5  M=1,3
	  R(M,II) = RI(II,M)
  5       CONTINUE
CC		
	  Q(1) = B* QX
	  Q(2) = B* QY
	  Q(3) = BC* QZ
	  QQ = SQRT(QX*QX+QY*QY+QZ*QZ)
	  IF (QQ.NE.0.0) THEN
	  VQ(1) = QX/QQ
	  VQ(2) = QY/QQ
	  VQ(3) = QZ/QQ
	  ELSE
	  VQ(1) = 0.0
	  VQ(2) = 0.0
	  VQ(3) = 1.0
	  END IF
C	
	  CC=0.
C****************************************************
C      
  	  DO 4000 K1=1,na
       	  Q1=1.
          IF(K1.GT.1)Q1=-1.0
       	  DO 4000 K2=K1,na
       	  Q2=1.
      	  IF(K2.GT.1)Q2=-1.0
C      
	  DO 20 M=1,3
  20	  R1(M) = R(M,K2)-R(M,K1)
C
	  DO 3000 M1=1,3
	  DO 3000 M2=1,M1
CC
C     CALCULATE THE SUM IN RECIPROCAL SPACE
C
	  PHI1 = 0
	  DO 200 N1 = -JN ,JN
	  DO 200 N2 = -JN ,JN
	  DO 200 N3 = -JN ,JN
	  T(1)=B*(N1*PRL(1,1)+N2*PRL(1,2)+N3*PRL(1,3))
	  T(2)=B*(N1*PRL(2,1)+N2*PRL(2,2)+N3*PRL(2,3))
	  T(3)=B*(N1*PRL(3,1)+N2*PRL(3,2)+N3*PRL(3,3))
	  DO 105 M=1,3
 105      TQ(M)=T(M)+Q(M)
	  TQ2 = TQ(1)*TQ(1) + TQ(2)*TQ(2) + TQ(3)*TQ(3)
	  IF (TQ2.GT.1E-10) THEN
	  S  = T(1)*R1(1)  +  T(2)*R1(2)  +  T(3)*R1(3)
	  PHI1 = PHI1 +TQ(M1)*TQ(M2)*4*PI/TQ2*EXP(-TQ2/4/ETA)*CDEXP(-CI*S)
	  ELSE
          VVQ = VQ(1)**2 +VQ(2)**2 + VQ(3)**2
	  PHI1 = PHI1 + VQ(M1)*VQ(M2)*4.*PI/VVQ
	  END IF
 200      CONTINUE
C
C     EVALUATE CORRECTION SUM IN REAL SPACE
C
	  PHI2 = 0
	  DO 300 JH1= -JH,JH
	  DO 300 JH2= -JH,JH
	  DO 300 JH3= -JH,JH
	  LAT(1) = A/2.0*(JH1*PL(1,1)+JH2*PL(1,2)+JH3*PL(1,3))
	  LAT(2) = A/2.0*(JH1*PL(2,1)+JH2*PL(2,2)+JH3*PL(2,3))
	  LAT(3) = A/2.0*(JH1*PL(3,1)+JH2*PL(3,2)+JH3*PL(3,3))
	  DO 205 M=1,3
	  RL(M) = LAT(M)+R1(M)
 205      X(M) = SQRT(ETA)*RL(M)
	  X2 = X(1)**2. + X(2)**2. +X(3)**2.
	  X1 = SQRT(X2)
	  IF ( X1 .gt. 0.0 ) then
	  H1 = 3./(X1**3)*(1.-ERF(X1))+ 2./SQRT(PI)*(3./X2+2.)*EXP(-X2)
	  H1 = X(M1)*X(M2)*H1/X2*ETA ! add *ETA
	     IF (M1.EQ.M2) THEN
	     H2 = 1./(X1**3.0)*(1.-ERF(X1)) + 2./SQRT(PI)/X2*EXP(-X2)
	     ELSE
	        H2 = 0
	     END IF
	        H = H1-H2
	  ELSE IF (M1.EQ.M2 ) THEN
	        H = 4./3./SQRT(PI)
	  ELSE
	        H = 0.0 ! omit R at the center.
	  END IF
	  S1 = Q(1)*RL(1) +Q(2)*RL(2) + Q(3)*RL(3)
	  HAB = -ETA**(1.5) *H *CDEXP(CI*S1)
 300      PHI2 = PHI2 +HAB

	  CC(K1,K2,M1,M2) = (PHI1/va + PHI2)*Q1*Q2
C
 3000 CONTINUE
C
	  CC(K1,K2,1,2)=CC(K1,K2,2,1)
	  CC(K1,K2,1,3)=CC(K1,K2,3,1)
	  CC(K1,K2,2,3)=CC(K1,K2,3,2)
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
 4000 CONTINUE
CC
CC***********************************************************
C
        DO 5000 K2=1,na-1
        KK2= K2+1
       	DO 5000 K1=KK2,na
       	DO 5000 M1=1,3
	DO 5000 M2=1,3
 5000   CC(K1,K2,M1,M2) = CONJG( CC(K2,K1,M2,M1))

CC**** CHECK TRANSLATIONAL INVARIANCE ***************
	IF(QQ.LT.0.00001)THEN
       	DO 1200 K1=1,na
       	DO 1200 I=1,3
	DO 1200 J=1,3
	DD=0.
       	   DO 1100 K2=1,na
              DD=DD-CC(K1,K2,I,J)
 1100      CONTINUE
	CC1(K1,I,J)=DD
 1200   CONTINUE
	END IF
C
C  CORRECT SELF-INTERACTION TO ENSURE TRANSLATIONAL INVARIANCE
        DO 1300 I=1,3
	DO 1300 J=1,3
       	DO 1300 K1=1,na
      	   CC(K1,K1,I,J)=CC(K1,K1,I,J)+CC1(K1,I,J)
 1300   CONTINUE
   	RETURN
	END SUBROUTINE RIMCOL
