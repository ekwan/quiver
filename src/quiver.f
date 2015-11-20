
C
C      QUIVER (KWIV'ER), V., N.
C       1: SHAKE; SHIVER; TREMBLE. [CF. OE CWIFERLICE ACTIVELY]
C       2: TO SHAKE WITH A SLIGHT BUT RAPID MOTION.
C
C                  KEITH E LAIDIG    YALE UNIVERSITY 1988
C
C       November 2015: modified to compile correctly with Singleton
C                      modifications, Eugene Kwan
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (MNA = 150)
        PARAMETER (MQ = 3*MNA)
c       PARAMETER (MMOL = 6)
c       PARAMETER (MMOL = 15)
        PARAMETER (MMOL = 42)
        PARAMETER (MTEM = 8)
        PARAMETER (NNREL = 10)
        CHARACTER*80 FILE, FILE1, TITLE
        DIMENSION C(MQ),E(MQ,MQ),F(MQ),A(MMOL,MQ),V(MQ,MQ),AI(3,3)
        DIMENSION VM(MQ,MQ,MMOL),X(MMOL,MQ),Y(MMOL,MQ),Z(MMOL,MQ)
        dimension xmem(mna),ymem(mna),zmem(mna),atw(600)
        integer*2 mwt(mna)
        DIMENSION W(3,3)
        DIMENSION RMASFA(NNREL),TRFA(NNREL),IREL(2,MMOL)
C
C    TEMPERATURE DEPENDENT VALUES
C
        DIMENSION TEMS(MTEM),TRRULE(NNREL,MTEM),PRODFA(NNREL,MTEM)
        DIMENSION TMASS(MMOL),PTL(MMOL),TMASFA(MMOL),PMU(MMOL,MTEM)
        DIMENSION SUMNU(MMOL),SSUMNU(MMOL),SUMU(MMOL,MTEM)
        DIMENSION BEXC(MMOL,MTEM),U(MQ),ONIUM(MMOL)
        DIMENSION RBEXC(NNREL,MTEM),DSUMU(NNREL,MTEM)
        DIMENSION QRAT(NNREL,MTEM),EFFR(NNREL,MTEM)
        character*2 sym
C
C    VARIOUS CONSTANTS TO BE USED
C    CM IS THE SPEED OF LIGHT IN CM/SEC
C    EM IS THE CONVERSION FACTOR FROM ERG TO KCAL/MOL
C    HM IS PLANCK'S CONSTANT
C    HBC IS PLANCK'S CONSTANT/BOLTZMANN'S CONSTANT (CM*K)
C
        DATA HM/6.625E-27/,CM/2.998E10/,EM/1.440E13/,HBC/1.4387/
        atw(1) = 1.00783
        atw(2) = 2.0141
        atw(3) = 3.016049
        atw(7) = 6.94
        atw(10) = 10.81
        atw(11) = 11.0
        atw(12) = 12.0
        atw(13) = 13.00335
        atw(14) = 14.0031
        atw(16) = 15.9949
        atw(17) = 16.9991
        atw(18) = 17.9992
        atw(19) = 18.9984
        atw(28) = 28.05
        atw(31) = 30.9738
        atw(32) = 31.9721
        atw(35) = 34.9689
        atw(37) = 36.9659
        atw(48) = 47.867
        atw(63) = 62.93
        atw(80) = 79.9
        atw(103) = 102.9055
        atw(106) = 106.4
        atw(190) = 190.2
        atw(192) = 192.2
        atw(101) = 101.0
        atw(127) = 126.904
        atw(15) = 15.0001
        atw(19) = 18.9984
        atw(48) = 47.947
C
C    READ IN INPUT AND OUTPUT FILE NAMES
C
        open(31,file=' freqs.',status='unknown')
        open(30,file=' ratios.',status='unknown')
C        WRITE (6,10)
C10      FORMAT(' NAME OF INPUT FILE ',$)
C        READ (5,20) FILE
        FILE='temp.q1'
C        WRITE (6,11)
C11      FORMAT(' NAME OF FORCE CONST. FILE ',$)
C        READ (5,20) FILE1
20      FORMAT(A80)
        FILE1='temp.q2'
        OPEN (10,STATUS='OLD',FILE=FILE)
        OPEN (40,STATUS='OLD',FILE=FILE1)
C
C        WRITE (6,30)
C30      FORMAT(' NAME OF OUTPUT FILE ',$)
C        READ (5,20) FILE
        FILE='temp.qout'
        OPEN (20,STATUS='NEW',FILE=FILE)
C
C    SHORT OR LONG OUTPUT?
C
C        WRITE (6,40)
C40      FORMAT(' SHORT OR LONG OUTPUT? [0/1] ',$)
C        READ (5,*) IPR
        IPR=1
C
        READ (10,20) TITLE
        WRITE(20,20) TITLE
        WRITE(6,20) TITLE
C
C   READ IN NUMBER OF INPUT MOLECULES
C
        READ (10,*)  NMOL
C
C    LOOP OVER ISOTOPOMERS
C
        DO 330 IMOL = 1,NMOL
C
C    READ IN TITLE, NUMBER OF ATOMS, MASS AND CARTESIAN COORDINATES
C    FOR EACH ATOM
C
c        WRITE (20,50)
50       FORMAT(/,'==================================================')
         READ (10,20) TITLE
         WRITE (20,60) TITLE(1:60)
60       FORMAT(/,'ISOTOPOMER NAME ',A60)
         IF (IPR.EQ.1.and.imol.eq.1) WRITE (20,70)
70       FORMAT(/21x,'  X             Y             Z ')
C
C    LOOP OVER NUMBER OF ATOMS
C
          if(imol.eq.1) then
         READ (10,*)  NAT
         DO 90 I = 1,NAT
           READ (10,80) sym, Xmem(I),Ymem(I),Zmem(I)
80        FORMAT(a2,3F14.0)
c       type *, Xmem(I),Ymem(I),Zmem(I)
          x(imol,i) = xmem(i)
          y(imol,i) = ymem(i)
          z(imol,i) = zmem(i)
          IF (IPR .EQ. 1) WRITE (20,81)
     1     sym, X(IMOL,I),Y(IMOL,I),Z(IMOL,I)
81        FORMAT(1x,a2,10x,3F14.6)
90       CONTINUE
        endif
        read (10,92)(mwt(j),j=1,nat)
        write(6,92)(mwt(j),j=1,nat)
92      format(150i6)
        IF (IPR.EQ.1) WRITE (20,96)
96      format(/,'          Weights')
        do 95 i=1,nat
          x(imol,i) = xmem(i)
          y(imol,i) = ymem(i)
          z(imol,i) = zmem(i)
          II = 3*I
        A(IMOL,II) = atw(mwt(i))
          IF (IPR .EQ. 1) WRITE (20,97) i,atw(mwt(i))
97      format(i5,f12.5)
95      continue
C
C    LOAD THE MASSES INTO A VECTOR OF 3*NAT IN ORDER
C    TO FACILITATE MASS WEIGHTING OF POTENTIAL ENERGY
C    MATRIX
C
c       type *,' NAT = ',nat
         DO 100 I = 3,3*NAT,3
          A(IMOL,I-1) = A(IMOL,I)
          A(IMOL,I-2) = A(IMOL,I)
100      CONTINUE
C
C    READ IN SCALING FACTOR FOR FREQUENCIES
C
       IF (IMOL .EQ. 1) THEN
          READ (10,106) SCALEOFFREQ
106       FORMAT(8F10.0)

       WRITE (20,108) SCALEOFFREQ
       WRITE (6,108) SCALEOFFREQ
108       FORMAT(' SCALING FACTOR FOR FREQUENCIES ',F9.6)

       
          SCALE = SCALEOFFREQ*SCALEOFFREQ
       WRITE (20,110) SCALE
       WRITE (6,110) SCALE
110       FORMAT(' SCALING FACTOR FOR FORCE CONSTANT MATRIX ',F9.6)
       END IF

C
C    READ IN POTENTIAL ENERGY MATRIX OF SECOND DERIVATIVE FORCE
C    CONSTANTS
C
         IF (IMOL .EQ. 1) THEN
c        iloop = INT((3*NAT+4)/5)
c       type *,' Loop: ',iloop
c          DO 130 K=1, iloop
c           KSTART=(5*K)-4
c       write(17,116) 3*nat,kstart+4
c116    format(' 3*nat kstart+4 ',2i8)
           DO 136 I1 = 1, 3*NAT
              DO 135 I2 = 1, I1
                 READ (40,*) V(I1, I2)
135           CONTINUE
136        CONTINUE
c           DO 130 I=KSTART,3*NAT
c            READ (10,120,ERR=130) IN,(V(I,J),J=KSTART,KSTART+4)
c120         FORMAT(I4,5d14.0)
c       write(17,121) i,v(i,kstart)
c121    format(i3,f22.8)
c       write(17,122) IN,(V(I,J),J=KSTART,KSTART+4)
c122        FORMAT(I4,5d14.6)
c130       CONTINUE
C
C    FILL IN SYMMETRIC MATRIX
C
          DO 140 I=1,3*NAT
           DO 140 J=1,I
            V(J,I)=V(I,J)
140       CONTINUE
C
C     WE CONVERT HARTREES PER SQUARE BOHR TO MILLIDYNES PER ANGSTROM
C     AND SCALE USING INPUT SCALING FACTOR
C
C       GIVEN:  1 BOHR IS 5.29167E-01 ANGSTROM
C               1 HARTREE IS 4.35942E-11 ERG
C               1 ANGSTROM IS 1.0E-08 CENTIMETER
C               1 ERG IS 1 DYNE * 1 CM
C
          GDF = 435.942 / ( 5.29167 * 5.29167 )
          DO 150 I = 1, 3*NAT
           DO 150 J = 1, I
            V(I,J) = V(I,J) * GDF * SCALE
            V(J,I) = V(I,J)
150       CONTINUE
         END IF
C
C    READ IN TEMPERATURES TO BE USED IN THIS RUN
C
         IF (IMOL .EQ. 1) THEN
          READ (10,160) NOTEM
160       FORMAT(I4)
c       type *,' number of temperatures',notem
          READ (10,170) (TEMS(IT),IT=1,NOTEM)
170       FORMAT(8F10.0)
          WRITE(6,172) (TEMS(IT),IT=1,NOTEM)
172       FORMAT(' TEMPERATURE                     ',F9.4)

         END IF
C
C    TRANSFER FORCE CONSTANTS TO MATRIX VM AND MASS WEIGHT THEM
C
         DO 180 I = 1,3*NAT
          DO 180 J = 1,3*NAT
           VM(J,I,IMOL) = V(J,I)/(DSQRT(A(IMOL,I))*DSQRT(A(IMOL,J)))
180      CONTINUE
C
C    DIAGONaLIZE THE MASS WEIGHTED CARTESIAN COORDINATES
C    THIS PRODUCES THE EIGENVALUES AND EIGENVECTORS
C    WHICH CORRESPOND TO THE FREQUENCIES (AFTER CORRECTION)
C    AND THE NORMAL MODE VIBRATIONAL DISPLACEMENTS, RESPECTIVELY
C
         CALL HDIAG (VM(1,1,IMOL),3*NAT,MQ,0,E)
C
C    CONVERT EIGENVALUES TO WAVENUMBERS
C
C    THE CONVERSION FACTOR USED IS
C
         CNV = SQRT(5.8890E-07)
C
         DO 200 I = 1,3*NAT
          F(I) = SIGN((DSQRT(ABS(VM(I,I,IMOL)))/CNV),VM(I,I,IMOL))
          DO 190 J = 1,3*NAT
           C(J) = E(J,I)/DSQRT(A(IMOL,J))
190       CONTINUE
C
C    USE A CUTOFF OF 50 CM-1 TO SEPARATE THE SIX NEAR ZERO
C    FREQUENCIES FROM THE OTHER 3N-6
C
        if(imol.eq.1)write(31,193)i,f(i)
193     format(i6,f14.3)
          IF (F(I) .GT. 40.0 .OR. F(I) .LT. -40.0) THEN
           IF (IPR .EQ. 1) WRITE (20,210) F(I)
           IF (imol .eq. 1) WRITE (6,210) F(I)
           IF (IPR .EQ. 1) WRITE (20,220) (C(K),K=1,3*NAT)
          END IF
200      CONTINUE
210      FORMAT(' FREQUENCY OF ', F15.5, ' 1/CM ')
220      FORMAT(' NORMAL COORDINATE ',/,3(3X,3F7.4))
C
C    GET MOMENTS OF INERTIA FOR THE MOLECULE
C
C    ZERO OUT ARRAY AI
C
         DO 230 I = 1,3
          DO 230 J = 1,3
           AI(J,I) = 0.0D0
230      CONTINUE
C
C    CONVERT CARTESIAN COORDINATES INTO ATOMIC UNITS
C    CONVERSION FACTOR IS 1.889727 BOHR/ANGSTROM
C
         ATB = 1.889727
C
         DO 240 J = 1,NAT
          X(IMOL,J) = X(IMOL,J)*ATB
          Y(IMOL,J) = Y(IMOL,J)*ATB
          Z(IMOL,J) = Z(IMOL,J)*ATB
240      CONTINUE
C
C    LOAD UP ARRAY AI FOR DIAGONALIZATION
C
         DO 250 J=1,NAT
          AI(1,1) = AI(1,1) + A(IMOL,3*J)*(Y(IMOL,J)**2+Z(IMOL,J)**2)
          AI(2,2) = AI(2,2) + A(IMOL,3*J)*(X(IMOL,J)**2+Z(IMOL,J)**2)
          AI(3,3) = AI(3,3) + A(IMOL,3*J)*(X(IMOL,J)**2+Y(IMOL,J)**2)
          AI(1,2) = AI(1,2) - A(IMOL,3*J)*X(IMOL,J)*Y(IMOL,J)
          AI(1,3) = AI(1,3) - A(IMOL,3*J)*X(IMOL,J)*Z(IMOL,J)
          AI(2,3) = AI(2,3) - A(IMOL,3*J)*Y(IMOL,J)*Z(IMOL,J)
250      CONTINUE
C
C    FILL IN SYMMETRIC OFF-DIAGONAL ELEMENTS
C
         AI(2,1) = AI(1,2)
         AI(3,1) = AI(1,3)
         AI(3,2) = AI(2,3)
C
C    DIAGONALIZE AI
C    EIGENVALUES RETURNED ARE THE MOMEMNTS OF INTERTIA AND
C    THE EIGENVECTORS REPRESENT THE INTERTIAL FRAME
C
         CALL HDIAG(AI,3,3,0,W)
C
c        WRITE (20,260) (AI(J,J), J=1,3)
260      FORMAT (/,' MOMENTS OF INERTIA'/' IX = ', 1PE12.4,
     1                               '  IY = ', 1PE12.4,
     2                               '  IZ = ', 1PE12.4 )
C
C    DETERMINE THE TOTAL MASS
C
C    SET PTL EQUAL TO ONE
C
         PTL(IMOL) = 1.0D0
         DO 270 I = 1,NAT
          TMASS(IMOL) = TMASS(IMOL) + A(IMOL,3*I)
          PTL(IMOL) = PTL(IMOL) + 1.5*LOG(A(IMOL,3*I))
270      CONTINUE
C%%%
c        WRITE(20,280) TMASS(IMOL)
280      FORMAT(/,' TOTAL MASS ',F8.4)
C
C   DETERMINE THE MASS FACTOR FOR EACH MOLECULE
C   !!!AVOID DIVIDE BY ZERO FOR LINEAR MOLECULES!!!
C
         IF (AI(1,1) .EQ. 0.0) AI(1,1) = 1.0D0
         IF (AI(2,2) .EQ. 0.0) AI(2,2) = 1.0D0
         IF (AI(3,3) .EQ. 0.0) AI(3,3) = 1.0D0
C
         TMASFA(IMOL) = (AI(1,1)*AI(2,2)*AI(3,3))**0.5
     1              *(TMASS(IMOL))**1.5
C%%%
C
C    TRAP OUT SMALL AND NEGATIVE FREQUENCIES HERE
C
         DO 290 I = 1,3*NAT
          IF (F(I) .GT.  50.0) SUMNU(IMOL) = SUMNU(IMOL) + F(I)
          IF (F(I) .LT. -50.0) ONIUM(IMOL) = -F(I)
          SSUMNU(IMOL) = SSUMNU(IMOL) + SIGN(F(I)**2,F(I))
290      CONTINUE
C
C    NOW FACTOR IN APPROPRIATE CONVERSION FACTORS TO THE
C    ZERO POINT ENERGY
C
         SUMNU(IMOL) = SUMNU(IMOL)*0.5*HM*CM*EM
C
c        WRITE (20,300) SUMNU(IMOL)
300      FORMAT(/,' THE ZERO-POINT ENERGY IS ',F8.4,' KCAL/MOLE')
C
C    FOR EACH TEMPERATURE INPUT,
C    1) CALCULATE THE PRODUCT SUM OF ONE MINUS EXP(-MU(I)) OVER THE
C       NORMAL MODES; GIVEN MU = H*DV(I)/K*T OR 1.4387*DV(I)/T.
C    2) CALCULATE THE SUM OF MU/2.0 OVER THE NORMAL MODES
C
         DO 320 L = 1,NOTEM
          BEXC(IMOL,L) = 1.0D0
          PMU(IMOL,L) = 1.0D0
          SUMU(IMOL,L) = 1.0D0
          DO 310 I = 1,3*NAT
           U(I) = HBC*F(I)/TEMS(L)
           IF (F(I) .GT. 50.0) BEXC(IMOL,L)=BEXC(IMOL,L)*(1.0-EXP(-U(I))
     1)
           IF (F(I) .GT. 50.0) PMU(IMOL,L) = PMU(IMOL,L)+LOG(U(I))
           IF (F(I) .GT. 50.0) SUMU(IMOL,L) = SUMU(IMOL,L) + U(I)
310       CONTINUE
C%%%
C
320      CONTINUE
C
C    READ IN FURTHER ISOTOPOMERS
C
330     CONTINUE
C
C    LOOP OVER NUMBER OF ISOTOPIC RELATIONSHIPS
C
340     READ (10,160) NIREL
C
C    LOOP OVER THE RELATIONSHIPS OF INTEREST
C
        DO 1440 JJ = 1,NIREL
C
         READ (10,350) (IREL(J,JJ),J=1,2)
350      FORMAT(2I4)
         WRITE(20,50)
         WRITE(20,360) IREL(2,JJ),IREL(1,JJ)
360      FORMAT(/,' RELATIONSHIPS BETWEEN ISOTOPOMERS ',I4,'    /',I4)
C
C    FIND TELLER-REDLICH PRODUCT RULE RATIO, RATIO OF MASS FACTORS,
C    FUNCTION OF SUMMED LOGS OF ATOMIC WEIGHTS AND SUMMED LOGS OF
C    FREQUENCIES
C
         RMASFA(JJ) = TMASFA(IREL(2,JJ))/TMASFA(IREL(1,JJ))
         TRFA(JJ) = RMASFA(JJ)*EXP(PTL(IREL(1,JJ))-PTL(IREL(2,JJ)))
C%%@@
C
C   LOOP OVER TEMPERATURES OF INTEREST
C
         DO 440 L = 1,NOTEM
C
C     RBEXC IS THE EXCITATION FACTOR
C
          PRODFA(JJ,L) = EXP(PMU(IREL(2,JJ),L)-PMU(IREL(1,JJ),L))
          RBEXC(JJ,L) = BEXC(IREL(1,JJ),L)/BEXC(IREL(2,JJ),L)
          DSUMU(JJ,L) = 0.5*(SUMU(IREL(1,JJ),L) - SUMU(IREL(2,JJ),L))
C
C    TELLER-REDLICH PRODUCT RULE
C
          TRRULE(JJ,L) = PRODFA(JJ,L)/TRFA(JJ)
          WRITE(20,365) L
365       FORMAT(/,' I am here now                ',I4)
C
C     EXP(DSUMU) IS THE ZERO-POINT FACTOR
C
          QRAT(JJ,L) = RMASFA(JJ)*RBEXC(JJ,L)*EXP(DSUMU(JJ,L))
C
C     EFFER IS THE BIGELEISEN-MAYER FUNCTION : (S2/S1)F
C
          EFFR(JJ,L) = PRODFA(JJ,L)*RBEXC(JJ,L)*EXP(DSUMU(JJ,L))
C%%
          WRITE(20,370) TEMS(L)
370       FORMAT(/,' TEMPERATURE                     ',F9.4)
          WRITE(20,380) PRODFA(JJ,L)
380       FORMAT(' RATIO OF PRODUCT OF FREQUENCIES ',E15.8)
          WRITE(20,390) BEXC(JJ,L)
390       FORMAT(' EXCITATION FACTOR               ',E15.8)
          WRITE(20,395) RBEXC(JJ,L)
395       FORMAT(' RBEXC                           ',E15.8)
          WRITE(20,400) EXP(DSUMU(JJ,L))
400       FORMAT(' ZERO-POINT DIFFERENCE FACTOR    ',E15.8)
          WRITE(20,410) TRRULE(JJ,L)
410       FORMAT(' TELLER-REDLICH PRODUCT RATIO    ',E15.8)
          WRITE(20,420) EFFR(JJ,L)
420       FORMAT(' (S2/S1)F                        ',E15.8)
c         WRITE(30,422) JJ,EFFR(JJ,L)
422       FORMAT(i5,E15.8)
C     WRITE(20,114) QRAT(JJ,L)
430       FORMAT(' RATIO OF PARTITION FUNCTIONS ',E15.8)
C
440     CONTINUE
1440    CONTINUE
        sum = 0.0
        sum1 = 0.0
        sum2 = 0.0
        do 500 i=1,nirel
        if(i.le.nirel/2) then
        sum1 = sum1 + effr(i,1)
        else
        sum2 = sum2 + effr(i,1)
        endif
500     sum = sum + effr(i,1)
        sum = sum/nirel
        do 510 i=1,nirel
        fnorm = effr(i,1)/sum
        write(30,522)i,fnorm
522       FORMAT(i5,f15.8)
510     continue
        ratio = sum1/sum2
        write(30,523)ratio
523     format(/,'  Ratio = ',f12.7)
C
450     STOP ' END OF QUIVER ANALYSIS '
C
        END
C
        SUBROUTINE HDIAG ( A, N, NDIM, IEGEN, EIVR )
        IMPLICIT NONE
C*****
C     THIS ROUTINE IS A MODIFICATION OF HDIAG, WRITTEN AT MIT,
C     THAT USES THE JACOBI METHOD.  THIS ROUTINE USES A VARIABLE
C     THRESHOLD JACOBI METHOD.  IT GIVES VERY GOOD EIGENVALUES
C     AND EIGENVECTORS AND IS MUCH FASTER.  THIS IS JACVAT.
C*****
        INTEGER I, I2, IEGEN, IFLAG, II, IROW, J, JCOL, JCOL1, JJ
        INTEGER K, N, NDIM, NU
        DOUBLE PRECISION A(NDIM,NDIM), AII, AIJ, AJJ, AMIN, ATOP, AVGF
        DOUBLE PRECISION C
        DOUBLE PRECISION D, DSTOP, EIVR(NDIM,NDIM), S, T, TEMP, THRSH, U
C*****
        IF ( N .LE. 1 ) THEN
         EIVR(1,1) = 1.0
         RETURN
        ENDIF
        IF ( IEGEN .EQ. 0 ) THEN
         DO 20 J = 1, N
          DO 10 I = 1, N
           EIVR(I,J) = 0.0
10        CONTINUE
          EIVR(J,J) = 1.0
20       CONTINUE
        ENDIF
C*****
C     FIND THE ABSOLUTELY LARGEST ELEMENT OF A
C*****
        ATOP = 0.0
        DO 30 I = 1, N
         DO 30 J = I, N
          IF ( ATOP .LT. ABS(A(I,J)) ) ATOP = ABS(A(I,J))
30      CONTINUE
        IF ( ATOP .LE. 0.0 ) RETURN
C*****
C     CALCULATE THE STOPPING CRITERION -- DSTOP
C*****
        AVGF = FLOAT ( N * ( N - 1 )) * 0.55
        D = 0.0
        DO 40 JJ = 2, N
         DO 40 II = 2, JJ
          S = A(II-1,JJ) / ATOP
          D = S * S + D
40      CONTINUE
        DSTOP = 1.0E-06 * D
C*****
C     CALCULATE THE THRESHOLD, THRSH
C*****
        THRSH = SQRT ( D / AVGF ) * ATOP
C*****
C     START A SWEEP
C*****
50      IFLAG = 0
        DO 110 JCOL = 2, N
         JCOL1 = JCOL - 1
         DO 110 IROW = 1, JCOL1
          AIJ = A(IROW,JCOL)
C*****
C     COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH
C*****
          IF ( ABS(AIJ) .LE. THRSH ) GO TO 110
          AII = A(IROW,IROW)
          AJJ = A(JCOL,JCOL)
          S = AJJ - AII
C*****
C     CHECK TO SEE IF  THE CHOSEN ROTATION IS LESS THAN THE ROUNDING
C     IF SO, THEN DO NOT ROTATE.
C*****
          IF ( ABS(AIJ) .LE. 1.0E-09 * ABS(S) ) GO TO 110
          IFLAG = 1
C*****
C     IF THE ROTATION IS VERY CLOSE TO 45 DEGREES, SET SIN AND COS
C        TO 1/(ROOT 2).
C*****
          IF ( 1.0E-10 * ABS(AIJ) .GE. ABS(S) ) THEN
           S = .707106781
           C = S
          ELSE
C*****
C     CALCULATION OF SIN AND COS FOR A ROTATION THAT IS NOT VERY CLOSE
C     TO 45 DEGREES.
C     COS = C, SIN = S
C*****
           T = AIJ / S
           S = 0.25 / SQRT ( 0.25 + T * T )
           C = SQRT ( 0.5 + S )
           S = 2.0 * T * S / C
          ENDIF
C*****
C     CALCULATION OF THE NEW ELEMENTS OF MATRIX A
C*****
          DO 60 I = 1, IROW
           T = A(I,IROW)
           U = A(I,JCOL)
           A(I,IROW) = C * T - S * U
           A(I,JCOL) = S * T + C * U
60        CONTINUE
          I2 = IROW + 2
          IF ( I2 .LE. JCOL ) THEN
           DO 70 I = I2, JCOL
            T = A(I-1,JCOL)
            U = A(IROW,I-1)
            A(I-1,JCOL) = S * U + C * T
            A(IROW,I-1) = C * U - S * T
70         CONTINUE
          ENDIF
          A(JCOL,JCOL) = S * AIJ + C * AJJ
          A(IROW,IROW) = C * A(IROW,IROW) -
     1                 S * ( C * AIJ - S * AJJ )
          DO 80 J = JCOL, N
           T = A(IROW,J)
           U = A(JCOL,J)
           A(IROW,J) = C * T - S * U
           A(JCOL,J) = S * T + C * U
80        CONTINUE
C*****
C     ROTATION COMPLETED.
C     SEE IF EIGENVECTORS ARE WANTED BY USER.
C*****
          IF ( IEGEN .EQ. 0 ) THEN
           DO 90 I = 1, N
            T = EIVR(I,IROW)
            EIVR(I,IROW) = C * T - EIVR(I,JCOL) * S
            EIVR(I,JCOL) = S * T + EIVR(I,JCOL) * C
90         CONTINUE
          ENDIF
C*****
C     CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP
C*****
          S = AIJ / ATOP
          D = D - S * S
          IF ( D .LT. DSTOP ) THEN
C*****
C     RECALCULATE DSTOP AND THRSH TO DISCARD ROUNDING ERRORS
C*****
           D = 0.0
           DO 100 JJ = 2, N
            DO 100 II = 2, JJ
             S = A(II-1,JJ) / ATOP
             D = S * S + D
100        CONTINUE
           DSTOP = 1.0E-06 * D
          ENDIF
          THRSH = SQRT ( D / AVGF ) * ATOP
110     CONTINUE
        IF ( IFLAG .NE. 0 ) GO TO 50
C*****
C     ARRANGE THE EIGENVALUES IN THE ORDER OF INCREASING ENERGY.
C     ARRANGE THE EIGENVECTORS IN THE SAME ORDER.
C*****
        NU = N
        DO 140 I = 1, N
         IF ( I .GE. NU ) RETURN
         AMIN = A(I,I)
         DO 130 J = I, NU
          IF ( A(J,J) .GE. AMIN ) GO TO 130
C*****
C     IF IEGEN IS -1, EXCLUDE UNCONVERGED EIGENVALUES FROM THIS ORDERING
C*****
          TEMP = ABS (EIVR(N,J)) + ABS (EIVR(N-1,J))
          IF ( TEMP .GT. 0.05 .AND. IEGEN .EQ. -1 ) THEN
           TEMP = A(J,J)
           A(J,J) = A(NU,NU)
           A(NU,NU) = TEMP
           II = NU
           NU = NU - 1
          ELSE
           II = I
           AMIN = A(J,J)
           A(J,J) = A(I,I)
           A(I,I) = AMIN
          ENDIF
          DO 120 K = 1, N
           TEMP = EIVR(K,II)
           EIVR(K,II) = EIVR(K,J)
           EIVR(K,J) = TEMP
120       CONTINUE
130      CONTINUE
140     CONTINUE
        RETURN
        END
