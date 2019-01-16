C     SUBROUTINE TO CALCULATE DYNAMICAL MATRIX FOR
C     SHORT RANGE PART OF THE INTERACTION FOR RIM
C
      SUBROUTINE RIMSR(A,C,QX,QY,QZ,DSR,DDR,pshr)
c---  for 2nd n.n. int.
      IMPLICIT NONE
      INTEGER, PARAMETER :: ndsr=14, na=2
      REAL*8 PI, A1, B1, C1, D1, E1, F1, C2, D2, E2, F2
      REAL*8 Q1, Q2, Q, PI2, DOT, S1, S2, C, A, QX, QY, QZ
      REAL*8 X1(3),DDR(na,3,3),pshr(ndsr)

      DIMENSION IRN(3,4,3),IR(3,4) ! IRN: coordinates for
                                   !      second nearest neighbour atoms
                                   ! IR : coordinates for
                                   !      the first nearest neighbour atoms
      INTEGER IRN, IR, I, J, IB, JB, II, IH, K1, K2
      COMPLEX*16 CI,CPH1,CPH(4,3),DSR(na,na,3,3),DD
C     UNITS IN L.V. or position vectors : (A,A,C)
C     UNITS IN R.L.V. or wave vectors : 2*pi*(1/A,1/A,1/C)
      DATA irn/2,2,0, -2,-2,0, 2,-2,0, -2,2,0
     &,2,0,2, -2,0,-2, 2,0,-2, -2,0,2
     &,0,2,2, 0,-2,-2, 0,2,-2, 0,-2,2/
      ! IRN: coordinates for the second nearest neighbour atoms
      !      in unit of A/4
      DATA ir/1,1,1, -1,1,-1, 1,-1,-1, -1,-1,1/
      ! IR : coordinates for the first nearest neighbour atoms
      ! in unit of A/4
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
!     Note: for group IV, E1 and E2 are NOT zero!
C     V1 = MAGNITUDE OF WAVEVECTOR
C     VQ = DIRECTION OF UNIT WAVEVECTOR
C
      Q2 = QX*QX + QY*QY + QZ*QZ
      Q1 = 2.*PI/A*SQRT(Q2)
C---  INITIALIZE THE DSR MATRICES ------
      DO 12 I=1,3
      DO 12 J=1,3
      DO 2 IB=1,na
      DO 2 JB=1,na
      DSR(IB,JB,I,J)=0.
2     CONTINUE
12    CONTINUE
C====on-ion Interaction ================

C-------  C-A  bonds ----
      DO 110 II=1,4
         DO 100 I=1,3
         X1(I)=IR(I,II) ! real space coordinates for the nearest atom
 100     CONTINUE
      DOT=(QX*X1(1)+QY*X1(2)+QZ*X1(3))*PI/2.0
      CPH1=CDEXP(CI*DOT)
        IF(II.EQ.1) THEN
c---  (1,1,1) --
          DSR(1,2,1,1)=DSR(1,2,1,1)+A1*CPH1 ! first 2 elements refer to atoms
          DSR(1,2,2,2)=DSR(1,2,2,2)+A1*CPH1 ! last 2 elements refer to positions
          DSR(1,2,3,3)=DSR(1,2,3,3)+A1*CPH1 ! x, y, z
          DSR(1,2,1,3)=DSR(1,2,1,3)+B1*CPH1
          DSR(1,2,3,1)=DSR(1,2,3,1)+B1*CPH1
          DSR(1,2,2,3)=DSR(1,2,2,3)+B1*CPH1
          DSR(1,2,3,2)=DSR(1,2,3,2)+B1*CPH1
          DSR(1,2,1,2)=DSR(1,2,1,2)+B1*CPH1
          DSR(1,2,2,1)=DSR(1,2,2,1)+B1*CPH1
        END IF
c---  (-1,1,-1) --
        IF(II.EQ.2) THEN
          DSR(1,2,1,1)=DSR(1,2,1,1)+A1*CPH1
          DSR(1,2,2,2)=DSR(1,2,2,2)+A1*CPH1
          DSR(1,2,3,3)=DSR(1,2,3,3)+A1*CPH1
          DSR(1,2,1,3)=DSR(1,2,1,3)+B1*CPH1
          DSR(1,2,3,1)=DSR(1,2,3,1)+B1*CPH1
          DSR(1,2,2,3)=DSR(1,2,2,3)-B1*CPH1
          DSR(1,2,3,2)=DSR(1,2,3,2)-B1*CPH1
          DSR(1,2,1,2)=DSR(1,2,1,2)-B1*CPH1
          DSR(1,2,2,1)=DSR(1,2,2,1)-B1*CPH1
        END IF
c---  (1,-1,-1) --
        IF(II.EQ.3) THEN
          DSR(1,2,1,1)=DSR(1,2,1,1)+A1*CPH1
          DSR(1,2,2,2)=DSR(1,2,2,2)+A1*CPH1
          DSR(1,2,3,3)=DSR(1,2,3,3)+A1*CPH1
          DSR(1,2,1,3)=DSR(1,2,1,3)-B1*CPH1
          DSR(1,2,3,1)=DSR(1,2,3,1)-B1*CPH1
          DSR(1,2,2,3)=DSR(1,2,2,3)+B1*CPH1
          DSR(1,2,3,2)=DSR(1,2,3,2)+B1*CPH1
          DSR(1,2,1,2)=DSR(1,2,1,2)-B1*CPH1
          DSR(1,2,2,1)=DSR(1,2,2,1)-B1*CPH1
        END IF
c---  (-1,-1,1) --
        IF(II.EQ.4) THEN
          DSR(1,2,1,1)=DSR(1,2,1,1)+A1*CPH1
          DSR(1,2,2,2)=DSR(1,2,2,2)+A1*CPH1
          DSR(1,2,3,3)=DSR(1,2,3,3)+A1*CPH1
          DSR(1,2,1,3)=DSR(1,2,1,3)-B1*CPH1
          DSR(1,2,3,1)=DSR(1,2,3,1)-B1*CPH1
          DSR(1,2,2,3)=DSR(1,2,2,3)-B1*CPH1
          DSR(1,2,3,2)=DSR(1,2,3,2)-B1*CPH1
          DSR(1,2,1,2)=DSR(1,2,1,2)+B1*CPH1
          DSR(1,2,2,1)=DSR(1,2,2,1)+B1*CPH1
        END IF
 110  CONTINUE
C-------  C-C bonds ----
      DO 120 IH=1,3
        DO II=1,4
           DO 102 I=1,3
           X1(I)=IRN(I,II,IH)
 102       CONTINUE
        DOT=QX*X1(1)+QY*X1(2)+QZ*X1(3)
        CPH(II,IH)=CDEXP(CI*DOT*PI/2.0)
        END DO
        IF(IH==1)THEN
c---  (110), (-1-10), (1-10), (-110)
c-------------correct---------------
          DO II=1,4
             IF (II==1) THEN
               s1=1.
               s2=1.
             ELSE IF (II==2) THEN
               s1=-1.
               s2=-1.
             ELSE IF (II==3) THEN
               s1=-1.
               s2= 1.
             ELSE
               s1= 1.
               s2=-1.
             END IF
c--------correct---------
c--------wrong----------
c     DO II=1,4
c     s2=1.
c     IF(II.eq.2.or.II.eq.3)s2=-1.
c     s1=1.
c--------wrong----------
             DSR(1,1,1,1)=DSR(1,1,1,1)+C1*CPH(II,IH)
             DSR(1,1,2,2)=DSR(1,1,2,2)+C1*CPH(II,IH)
             DSR(1,1,3,3)=DSR(1,1,3,3)+F1*CPH(II,IH)
             DSR(1,1,1,2)=DSR(1,1,1,2)+s1*s2*D1*CPH(II,IH)
             DSR(1,1,2,1)=DSR(1,1,2,1)+s1*s2*D1*CPH(II,IH)
             DSR(1,1,1,3)=DSR(1,1,1,3)+s1*E1*CPH(II,IH)
             DSR(1,1,3,1)=DSR(1,1,3,1)-s1*E1*CPH(II,IH)
             DSR(1,1,2,3)=DSR(1,1,2,3)+s2*E1*CPH(II,IH)
             DSR(1,1,3,2)=DSR(1,1,3,2)-s2*E1*CPH(II,IH)
          END DO
        END IF
        IF(IH==2)THEN
c---  (101), (-10-1), (10-1), (-101)
c----------correct----------------------
          DO II=1,4
             if (II==1) then
              s1=1.
              s2=1.
             else if (II==2) then
              s1=-1.
              s2=-1.
             else if (II==3) then
              s1=-1.
              s2= 1.
             else
              s1= 1.
              s2=-1.
             end if
c----------correct----------------------
c----------wrong----------------------
c     DO II=1,4
c     s2=1.
c     IF(II.eq.2.or.II.eq.3)s2=-1.
c     s1=1.
c     if(II.eq.2.or.II.eq.4)s1=-1.
c----------wrong----------------------
             DSR(1,1,1,1)=DSR(1,1,1,1)+C1*CPH(II,IH)
             DSR(1,1,2,2)=DSR(1,1,2,2)+F1*CPH(II,IH)
             DSR(1,1,3,3)=DSR(1,1,3,3)+C1*CPH(II,IH)
             DSR(1,1,1,3)=DSR(1,1,1,3)+s1*s2*D1*CPH(II,IH)
             DSR(1,1,3,1)=DSR(1,1,3,1)+s1*s2*D1*CPH(II,IH)
             DSR(1,1,1,2)=DSR(1,1,1,2)+s1*E1*CPH(II,IH)
             DSR(1,1,2,1)=DSR(1,1,2,1)-s1*E1*CPH(II,IH)
             DSR(1,1,2,3)=DSR(1,1,2,3)-s2*E1*CPH(II,IH)
             DSR(1,1,3,2)=DSR(1,1,3,2)+s2*E1*CPH(II,IH)
          END DO
        END IF
      if(IH==3)then
c--- (011), (0-1-1), (01-1), (0-11)
c----------correct----------------------
        DO II=1,4
        if (II==1) then
        s1=1.
        s2=1.
      else if (II==2) then
        s1=-1.
        s2=-1.
      else if (II==3) then
        s1=-1.
        s2= 1.
      else
        s1= 1.
        s2=-1.
        end if
c----------correct----------------------
c----------wrong----------------------
c   DO II=1,4
c   s2=1.
c   IF(II.eq.2.or.II.eq.3)s2=-1.
c   s1=1.
c   if(II.eq.2.or.II.eq.4)s1=-1.
c----------wrong----------------------
        DSR(1,1,1,1)=DSR(1,1,1,1)+F1*CPH(II,IH)
        DSR(1,1,2,2)=DSR(1,1,2,2)+C1*CPH(II,IH)
        DSR(1,1,3,3)=DSR(1,1,3,3)+C1*CPH(II,IH)
        DSR(1,1,2,3)=DSR(1,1,2,3)+s1*s2*D1*CPH(II,IH)
        DSR(1,1,3,2)=DSR(1,1,3,2)+s1*s2*D1*CPH(II,IH)
        DSR(1,1,1,3)=DSR(1,1,1,3)-s2*E1*CPH(II,IH)
        DSR(1,1,3,1)=DSR(1,1,3,1)+s2*E1*CPH(II,IH)
        DSR(1,1,1,2)=DSR(1,1,1,2)-s1*E1*CPH(II,IH)
        DSR(1,1,2,1)=DSR(1,1,2,1)+s1*E1*CPH(II,IH)
        END DO
      endif
C--------  A-A  bonds ------
C------ DSR(3,3,1,1)=DSR(3,3',1,1),...etc------
c        DOTX=QX*R1(1)/A
c        CPH1X=CDEXP(CI*DOTX*PI2)
c        DOTY=QY*R2(2)/A
c        CPH1Y=EXP(CI*DOTY*PI2)
      if(IH==1)then
c--- (110), (-1-10), (1-10), (-110)
c----------correct----------------------
        DO II=1,4
      if (II==1) then
        s1=1.
        s2=1.
      else if (II==2) then
        s1=-1.
        s2=-1.
      else if (II==3) then
        s1=-1.
        s2= 1.
      else
        s1= 1.
        s2=-1.
      end if
c----------correct----------------------
c----------wrong----------------------
c      DO II=1,4
c   s2=1.
c        IF(II.eq.2.or.II.eq.3)s2=-1.
c   s1=1.
c   if(II.eq.2.or.II.eq.4)s1=-1.
c----------wrong----------------------
        DSR(2,2,1,1)=DSR(2,2,1,1)+C2*CPH(II,IH)
        DSR(2,2,2,2)=DSR(2,2,2,2)+C2*CPH(II,IH)
        DSR(2,2,3,3)=DSR(2,2,3,3)+F2*CPH(II,IH)
        DSR(2,2,1,2)=DSR(2,2,1,2)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,2,1)=DSR(2,2,2,1)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,1,3)=DSR(2,2,1,3)-s1*E2*CPH(II,IH)
        DSR(2,2,3,1)=DSR(2,2,3,1)+s1*E2*CPH(II,IH)
        DSR(2,2,2,3)=DSR(2,2,2,3)-s2*E2*CPH(II,IH)
        DSR(2,2,3,2)=DSR(2,2,3,2)+s2*E2*CPH(II,IH)
        END DO
      endif
      if(IH==2)then
c---  (101), (-10-1), (10-1), (-101)
c----------correct----------------------
        DO II=1,4
        if (II==1) then
        s1=1.
        s2=1.
      else if (II==2) then
        s1=-1.
        s2=-1.
      else if (II==3) then
        s1=-1.
        s2= 1.
      else
        s1= 1.
        s2=-1.
        end if
c----------correct----------------------
c----------wrong----------------------
c      DO II=1,4
c     s2=1.
c        IF(II.eq.2.or.II.eq.3)s2=-1.
c     s1=1.
c      if(II.eq.2.or.II.eq.4)s1=-1.
c----------wrong----------------------
        DSR(2,2,1,1)=DSR(2,2,1,1)+C2*CPH(II,IH)
        DSR(2,2,2,2)=DSR(2,2,2,2)+F2*CPH(II,IH)
        DSR(2,2,3,3)=DSR(2,2,3,3)+C2*CPH(II,IH)
        DSR(2,2,1,3)=DSR(2,2,1,3)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,3,1)=DSR(2,2,3,1)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,1,2)=DSR(2,2,1,2)-s1*E2*CPH(II,IH)
        DSR(2,2,2,1)=DSR(2,2,2,1)+s1*E2*CPH(II,IH)
        DSR(2,2,2,3)=DSR(2,2,2,3)+s2*E2*CPH(II,IH)
        DSR(2,2,3,2)=DSR(2,2,3,2)-s2*E2*CPH(II,IH)

        END DO
      endif
      if(IH==3)then
c---  (011), (0-1-1), (01-1), (0-11)
c----------correct----------------------
       DO II=1,4
        if (II==1) then
        s1=1.
        s2=1.
      else if (II==2) then
        s1=-1.
        s2=-1.
      else if (II==3) then
        s1=-1.
        s2= 1.
      else
        s1= 1.
        s2=-1.
      end if
c----------correct----------------------
c----------wrong----------------------
c        DO II=1,4
c   s2=1.
c        IF(II.eq.2.or.II.eq.3)s2=-1.
c   s1=1.
c   if(II.eq.2.or.II.eq.4)s1=-1.
c----------wrong----------------------
        DSR(2,2,1,1)=DSR(2,2,1,1)+F2*CPH(II,IH)
        DSR(2,2,2,2)=DSR(2,2,2,2)+C2*CPH(II,IH)
        DSR(2,2,3,3)=DSR(2,2,3,3)+C2*CPH(II,IH)
        DSR(2,2,2,3)=DSR(2,2,2,3)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,3,2)=DSR(2,2,3,2)+s1*s2*D2*CPH(II,IH)
        DSR(2,2,1,3)=DSR(2,2,1,3)+s2*E2*CPH(II,IH)
        DSR(2,2,3,1)=DSR(2,2,3,1)-s2*E2*CPH(II,IH)
        DSR(2,2,1,2)=DSR(2,2,1,2)+s1*E2*CPH(II,IH)
        DSR(2,2,2,1)=DSR(2,2,2,1)-s1*E2*CPH(II,IH)
        END DO
      endif
 120  continue

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
      DD=DD-DSR(IB,JB,I,J)
 190  CONTINUE
      DDR(IB,I,J)=DD
 200  CONTINUE
      END IF
c--- INCLUDE SELF-INTERACTION
      DO 300 IB=1,na
      DO 300 I=1,3
      DO 300 J=1,3
      DSR(IB,IB,I,J)=DDR(IB,I,J)+DSR(IB,IB,I,J)
 300  CONTINUE

        RETURN
      END SUBROUTINE RIMSR
