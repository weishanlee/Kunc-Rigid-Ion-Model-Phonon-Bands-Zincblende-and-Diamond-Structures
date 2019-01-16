C     SUBROUTINE TO CALCULATE DYNAMICAL MATRIX FOR
C     SHORT RANGE PART OF THE INTERACTION FOR RIM
C
      SUBROUTINE RIMSR_paper(A,C,QX,QY,QZ,DSR,DDR,pshr)
c---    for 2nd n.n. int.
      implicit none
      integer, parameter :: ndsr=14,na=2
      real*8 PI, A1, B1, C1, D1, E1, F1, C2, D2, E2, F2
      real*8 Q1, Q2, Q, PI2, DOT, S1, S2, C, A, QX, QY, QZ
      REAL*8 X1(3),DDR(na,3,3),pshr(ndsr)
      DIMENSION IRN(3,4,3),IR(3,4) ! IRN: coordinates for
                                 !      second nearest neighbour atoms
                                 ! IR : coordinates for
                                 !      the first nearest neighbour atoms
      integer IRN, IR, I, J, IB, JB, II, IH, K1, K2
      COMPLEX*16 CI,CPH1,CPH(4,3),DSR(na,na,3,3),DD
C     UNITS IN L.V. or position vectors : (A,A,C)
C     UNITS IN R.L.V. or wave vectors : 2 pi (1/A,1/A,1/C)
      data irn/2,2,0, -2,-2,0, 2,-2,0, -2,2,0
     &,2,0,2, -2,0,-2, 2,0,-2, -2,0,2
     &,0,2,2, 0,-2,-2, 0,2,-2, 0,-2,2/
                ! IRN: coordinates for
                !      second nearest neighbour atoms
                !     in unit of A/4
      data ir/1,1,1, -1,1,-1, 1,-1,-1, -1,-1,1/
                ! IR : coordinates for
                !      the first nearest neighbour atoms
                !      in unit of A/4
      A=pshr(12)
      CI= (0.d0,1.d0)
      PI=dacos(-1.d0)
      PI2=2*PI

      a1=pshr(1)
      b1=pshr(2)
      c1=pshr(3)
      d1=pshr(4)
      e1=pshr(5)
      f1=pshr(6)
      C2=pshr(7)
      D2=pshr(8)
      E2=pshr(9)
      F2=pshr(10)
! Note: for group IV, we must have E1=E2=0 to preserve symmetry
C
C     V1 = MAGNITUDE OF WAVEVECTOR
C     VQ = DIRECTION OF UNIT WAVEVECTOR
C
      Q2 = QX*QX + QY*QY + QZ*QZ
      Q1 = 2.*PI/A*SQRT(Q2)
C--- INITIALIZE THE DSR MATRICES ------
      DO 12 I=1,3
      DO 12 J=1,3
      DO 2 IB=1,na
      DO 2 JB=1,na
      DSR(IB,JB,I,J)=0.
 2    CONTINUE
 12   CONTINUE
C
C======  Dynamic Matrix for Short Ranges ======================
!!!!!!!!!!!!!!!!!!!!!!!!!1111111111111!!!!!!!!!!!!!!!!!!!!!!!!!
      DSR(1,1,1,1)=-4.*( A1+
     & C1*( 2. - cos(PI*QX) * ( cos(PI*QY) + cos(PI*QZ) ) ) +
     & F1*( 1. - cos(PI*QY) * cos(PI*QZ) )
     & )

      DSR(1,1,2,2)=-4.*( A1+
     & C1*( 2. - cos(PI*QY) * ( cos(PI*QZ) + cos(PI*QX) ) ) +
     & F1*( 1. - cos(PI*QZ) * cos(PI*QX) )
     & )

      DSR(1,1,3,3)=-4.*( A1+
     & C1*( 2. - cos(PI*QZ) * ( cos(PI*QX) + cos(PI*QY) ) ) +
     & F1*( 1. - cos(PI*QX) * cos(PI*QY) )
     & )

      DSR(2,2,1,1)=-4.*( A1+
     & C2*( 2. - cos(PI*QX) * ( cos(PI*QY) + cos(PI*QZ) ) ) +
     & F2*( 1. - cos(PI*QY) * cos(PI*QZ) )
     & )

      DSR(2,2,2,2)=-4.*( A1+
     & C2*( 2. - cos(PI*QY) * ( cos(PI*QZ) + cos(PI*QX) ) ) +
     & F2*( 1. - cos(PI*QZ) * cos(PI*QX) )
     & )

      DSR(2,2,3,3)=-4.*( A1+
     & C2*( 2. - cos(PI*QZ) * ( cos(PI*QX) + cos(PI*QY) ) ) +
     & F2*( 1. - cos(PI*QX) * cos(PI*QY) )
     & )
!!!!!!!!!!!!!!!!!!!!!!!!!!22222222222222!!!!!!!!!!!!!!!!!!!!!!!!!
      DSR(1,1,1,2)=-4.*( D1 * sin(PI*QX) * sin(PI*QY) -
     & CI * E1 *sin(PI*QZ) * ( cos(PI*QX)-cos(PI*QY) )
     & )

      DSR(1,1,1,3)=-4.*( D1 * sin(PI*QX) * sin(PI*QZ) -
     & CI * E1 *sin(PI*QY) * ( cos(PI*QX)-cos(PI*QZ) )
     & )

      DSR(1,1,2,3)=-4.*( D1 * sin(PI*QY) * sin(PI*QZ) -
     & CI * E1 *sin(PI*QX) * ( cos(PI*QY)-cos(PI*QZ) )
     & )

      DSR(1,1,2,1)=-4.*( D1 * sin(PI*QY) * sin(PI*QX) -
     & CI * E1 *sin(PI*QZ) * ( cos(PI*QY)-cos(PI*QX) )
     & )

      DSR(1,1,3,1)=-4.*( D1 * sin(PI*QZ) * sin(PI*QX) -
     & CI * E1 *sin(PI*QY) * ( cos(PI*QZ)-cos(PI*QX) )
     & )

      DSR(1,1,3,2)=-4.*( D1 * sin(PI*QZ) * sin(PI*QY) -
     & CI * E1 *sin(PI*QX) * ( cos(PI*QZ)-cos(PI*QY) )
     & )

      DSR(2,2,1,2)=-4.*( D2 * sin(PI*QX) * sin(PI*QY) +
     & CI * E2 *sin(PI*QZ) * ( cos(PI*QX)-cos(PI*QY) )
     & )

      DSR(2,2,1,3)=-4.*( D2 * sin(PI*QX) * sin(PI*QZ) +
     & CI * E2 *sin(PI*QY) * ( cos(PI*QX)-cos(PI*QZ) )
     & )

      DSR(2,2,2,3)=-4.*( D2 * sin(PI*QY) * sin(PI*QZ) +
     & CI * E2 *sin(PI*QX) * ( cos(PI*QY)-cos(PI*QZ) )
     & )

      DSR(2,2,2,1)=-4.*( D2 * sin(PI*QY) * sin(PI*QX) +
     & CI * E2 *sin(PI*QZ) * ( cos(PI*QY)-cos(PI*QX) )
     & )

      DSR(2,2,3,1)=-4.*( D2 * sin(PI*QZ) * sin(PI*QX) +
     & CI * E2 *sin(PI*QY) * ( cos(PI*QZ)-cos(PI*QX) )
     & )

      DSR(2,2,3,2)=-4.*( D2 * sin(PI*QZ) * sin(PI*QY) +
     & CI * E2 *sin(PI*QX) * ( cos(PI*QZ)-cos(PI*QY) )
     & )
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!33333333333333333!!!!!!!!!!!!!!!!!!
      DSR(1,2,1,1)=4.*A1*(
     & cos(PI/2.*QX) * cos(PI/2.*QY) * cos(PI/2.*QZ) - CI*
     & sin(PI/2.*QX) * sin(PI/2.*QY) * sin(PI/2.*QZ)
     & )

      DSR(1,2,2,2)=4.*A1*(
     & cos(PI/2.*QY) * cos(PI/2.*QZ) * cos(PI/2.*QX) - CI*
     & sin(PI/2.*QY) * sin(PI/2.*QZ) * sin(PI/2.*QX)
     & )

      DSR(1,2,3,3)=4.*A1*(
     & cos(PI/2.*QZ) * cos(PI/2.*QX) * cos(PI/2.*QY) - CI*
     & sin(PI/2.*QZ) * sin(PI/2.*QX) * sin(PI/2.*QY)
     & )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!44444444444444444444444!!!!!!!!!!
      DSR(1,2,1,2)=4.*B1*(
     & -sin(PI/2.*QX) * sin(PI/2.*QY) * cos(PI/2.*QZ) + CI*
     &  cos(PI/2.*QX) * cos(PI/2.*QY) * sin(PI/2.*QZ)
     &  )

      DSR(1,2,1,3)=4.*B1*(
     & -sin(PI/2.*QX) * sin(PI/2.*QZ) * cos(PI/2.*QY) + CI*
     &  cos(PI/2.*QX) * cos(PI/2.*QZ) * sin(PI/2.*QY)
     &  )

      DSR(1,2,2,1)=4.*B1*(
     & -sin(PI/2.*QY) * sin(PI/2.*QX) * cos(PI/2.*QZ) + CI*
     &  cos(PI/2.*QY) * cos(PI/2.*QX) * sin(PI/2.*QZ)
     &  )

      DSR(1,2,2,3)=4.*B1*(
     & -sin(PI/2.*QY) * sin(PI/2.*QZ) * cos(PI/2.*QX) + CI*
     &  cos(PI/2.*QY) * cos(PI/2.*QZ) * sin(PI/2.*QX)
     &  )

      DSR(1,2,3,1)=4.*B1*(
     & -sin(PI/2.*QZ) * sin(PI/2.*QX) * cos(PI/2.*QY) + CI*
     &  cos(PI/2.*QZ) * cos(PI/2.*QX) * sin(PI/2.*QY)
     &  )

      DSR(1,2,3,2)=4.*B1*(
     & -sin(PI/2.*QZ) * sin(PI/2.*QY) * cos(PI/2.*QX) + CI*
     &  cos(PI/2.*QZ) * cos(PI/2.*QY) * sin(PI/2.*QX)
     &  )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 140 K1=2,na
      DO 140 K2=1,K1-1
      DO 140 I=1,3
      DO 140 J=1,3
         DSR(K1,K2,I,J)=CONJG(DSR(K2,K1,J,I))
 140  CONTINUE
C --- DEFINE SELF-INTERACTIONS FOR Q1=0. -----
      IF(Q1 .LT. 0.00001)THEN
      DO 200 IB=1,na
      DO 200 I=1,3
      DO 200 J=1,3
      DD=0.
      DO 190 JB=1,na
c   IF(JB.EQ.IB) GO TO 190  ! redundant
      DD=DD-DSR(IB,JB,I,J)    ! why negative?
                            ! sum over two atoms for each (x,y,z)
      !WRITE (*,*) DD
 190  CONTINUE
      DDR(IB,I,J)=DD
 200  CONTINUE
      END IF
c--- INCLUDE SELF-INTERACTION
      DO 300 IB=1,na
      DO 300 I=1,3
      DO 300 J=1,3
      DSR(IB,IB,I,J)=DDR(IB,I,J)+DSR(IB,IB,I,J)  ! DDR(1,1,1)=154.2296
 300  CONTINUE
! 7007  format(12f6.1)
c  440   CONTINUE
c        CLOSE(33)
      RETURN
      END SUBROUTINE RIMSR_paper
