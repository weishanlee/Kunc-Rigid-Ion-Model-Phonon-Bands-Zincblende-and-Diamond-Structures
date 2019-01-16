c------------------------------------------------------
C       PHONON DISPERSION IN rim IN GENERAL DIRECTIONS
C
C       PHASE USED IS e^(i*q*(x(l'k')-x(lk))
C
      SUBROUTINE RIMPDR(PSHR,FBD)
      IMPLICIT NONE
      INTEGER, PARAMETER :: ndsr=14, nk=801, na=2,nb=3*na
      REAL*8 DDR(na,3,3)
      REAL*8 Q1, PI, Z2, C, A, Va, QX, QY, QZ, IER, THz
      REAL*8 M(na),HR(nb,nb),HI(nb,nb),FREQ(nb),RI(na,3)
     &  ,W2(nb),ZR(nb,nb),ZI(nb,nb),FV1(nb),FV2(nb),FM1(2,nb)
     &  ,FBD(nb,nk),NKL(nk,3),PSHR(ndsr), CC1(na,3,3)
      INTEGER nkn, i, JN, JH,  IQ, J, K1, K2, M1, M2, IS
      INTEGER J1, J2, I1, I2
      COMPLEX*16 CI,CC(na,na,3,3),DSR(na,na,3,3)
      COMPLEX*16 DR(nb,nb)
      DO nkn=1,nk,1
         READ(25,*) ( NKL(nkn,i), i=1,3 )
      END DO

      PI=dacos(-1.d0)
      CI=(0.d0,1.d0)  ! complex I
      Z2=pshr(11)
      A=pshr(12)
      C=A
      Va= A**3/4.0 ! Va: volume of primitive unit cell
                   ! t1= A/2*(0,1,1); t2= A/2*(1,0,1); t3= A/2*(1,1,0).
c                  ! Va= t1 dot (t2 cross t3)
      M(1)=pshr(13)
      M(2)=pshr(14)
c--------to install ddr------------
      QX=0.
      QY=0.
      QZ=0.
      Q1=2.*PI/A*SQRT(QX**2+QY**2+QZ**2)
      CALL RIMSR(A,C,QX,QY,QZ,DSR,DDR,pshr)
c------ ddr installed-------------
      JN=6
      JH=6
      DATA RI/0.,0.5,0.,0.5,0.,0.5/  ! ionic coordinates in a/2

      RI=RI*A/2.
      CALL RIMCOL(JN,JH,QX,QY,QZ,RI,CC,CC1,A)
      DO IQ=1,nk
         QX=NKL(IQ,1)/20.
         QY=NKL(IQ,2)/20.
         QZ=NKL(IQ,3)/20.
         Q1=2.*PI/A*SQRT(QX**2+QY**2+QZ**2)

         CALL RIMCOL(JN,JH,QX,QY,QZ,RI,CC,CC1,A)

         CALL RIMSR(A,C,QX,QY,QZ,DSR,DDR,pshr)

         DO K1= 1,na
            DO K2= 1,na
               DO M1=1,3
                  DO M2=1,3
                     J1=M1+(K1-1)*3
                     J2=M2+(K2-1)*3
                     DR(J1,J2)= DSR(K1,K2,M1,M2)  + Z2*CC(K1,K2,M1,M2)
                     IF (J1 <= 3 .and. J2 <= 3) THEN
                     DR(J1,J2)= DR(J1,J2)/SQRT(M(K1)*M(K1))
                     ELSE IF (J1 > 3 .and. J2 > 3 ) THEN
                     DR(J1,J2)= DR(J1,J2)/SQRT(M(K2)*M(K2))
                     ELSE
                     DR(J1,J2)= DR(J1,J2)/SQRT(M(K1)*M(K2))
                     END IF
                  END DO
               END DO
            END DO
         END DO
C
         DO I1= 1,nb
            DO I2= 1,nb
               J1=I1
               J2=I2
               HR(I1,I2) = DR(J1,J2)
               HI(I1,I2) = -CI*DR(J1,J2)
            END DO
         END DO

         IS=0 ! IS = 0 gives eigenvalues only;
              ! IS = 1 gives both eigenvalues and eigenvectors.

         CALL CH(nb,nb,HR,HI,W2,IS,ZR,ZI,FV1,FV2,FM1,IER)
              !give HR and HI, then return W2 (frequency^2)
c         IF(IER.NE.0.0)WRITE(6,*) '  EISPACK ERROR, IER=',IER
C        CONVERT sqrt(e^2/Va/M)/2pi into tera Hz
C                e: esu, Va: amstrong cube, M: amu
         THZ= SQRT(4.803**2/(1.6605*VA))*100.0/(2.0*PI)
         DO J=1,6
            FREQ(J) = SQRT(ABS(W2(J)))* THZ
            FBD(J,IQ)=FREQ(J)        ! returned frequency
         END DO
      END DO
      RETURN
      END SUBROUTINE RIMPDR
