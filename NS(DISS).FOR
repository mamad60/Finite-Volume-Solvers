*     FINITE VOLUME METHOD OF SOLUTION OF THE TWO DIMENTIONAL                        
*             INCOMPRESSIBLE NAVIER-STOKES EQUATIONS FOR                                                
*                      FLOW AROUND CYLENDER                            
*                                                                    
*           BY ARTIFICIAL COMPRESSIBILITY METHOD. 
*..................................................................................................................................

         
*        LEGEND

*
*        X,Y                    CARTESIAN COORDINATES  
*        A                      AREA OF THE PRIMARY CELLS
*        U                      HORIZONTAL VELOCITY
*        V                      VERTICAL VELOCITY  
*        P                      PRESSURE
*        U0,V0,P0               INITIAL VALUES (INITIAL GUESS)
*        M                      NO. OF POINTS IN THE NORMAL TO THE BODY
*        N                      NO. OF POINTS AROUND CYLENDER 
*        CFL                    CORANT FACTOR  
*        DT                     DELTA(TIME)
*        ITR                    NO. OF ITERATION (INPUT)
*        ITER                   NO. OF RUNNING ITERATION  
*        R                      RADIUS OF CIRCLES
*        R(1)                   RADIUS OF CYLENDER 
*        BETA                   CLUSTERING PARAMETER
*        RADI                   RATIO OF THE OUTER BOUNDARY TO SOLID BOUNDARY
*        REN                    REYNOLDS NUMBER
*        ACOM                   ARTIFITIAL COMPRESSIBILITY
*        VOR                    VORTICITY
*        EU,EV,EP               ARTIFITIAL DISSIPATION PARAMETERS
*        CP                     PRESSURE COEF.
*        UT                     OVERAL VELOCITY
*        CD                     DRAG COEF.
*        PS                     PRESSURE ON THE SOLID BOUNDARY 
*        ERR,ENORM              ACCURCY PARAMETERS


*..................MAIN PROG.
         
         PARAMETER(NI=100,NJ=100)

         REAL*4     X(NI,NJ),Y(NI,NJ),A(NI,NJ) 
         REAL*4     U(NI,NJ),V(NI,NJ),P(NI,NJ)
         REAL*4     UT(NI,NJ),CP(NI,NJ),R(NJ)
         REAL*4     PS(NI),VOR(NI,NJ),REN,CFL

         COMMON/B0 /X,Y,A
         COMMON/B1 /U,V,P
         COMMON/B2/U0,V0,P0
         COMMON/B3 /M,JM,N
         COMMON/B4 /DT,CFL
         COMMON/B6 /ITER
         COMMON/B7 /R
         COMMON/B8 /BETA,RADI
         COMMON/B10/REN
         COMMON/B11/ACOM
         COMMON/B13/VOR
         COMMON/B14/EU,EV,EP
         COMMON/B17/CP,UT,CD
         COMMON/B19/PS

         OPEN(1,FILE='INPUT1.DAT',STATUS='UNKNOWN')
         OPEN(2,FILE= 'HIST.DAT',STATUS='UNKNOWN')

         READ(1,*)
         READ(1,*) BETA,RADI,ACOM
         READ(1,*)
         READ(1,*) EU,EV,EP
         READ(1,*)
         READ(1,*) U0,V0,P0
         READ(1,*)
         READ(1,*) M,N,R(1)
         READ(1,*)
         READ(1,*) CFL,ERR
         READ(1,*)
         !READ(1,*) REN,ITR
	   REN=40
	   ITR=1000
         
         JM=M-2
         
         CALL GRIDGEN
         CALL INITIALVAL
         
         WRITE(2,*)'    ITER','    ENORM'

*........START  ITERATION

         ITER=0
1        ITER=ITER+1

         WRITE(2,2)ITER, ENORM
2        FORMAT(2x,I4,2x,F10.7)

         CALL CFLTEST
         CALL RK4
         CALL ACCUR(ENORM)

         PI=4*ATAN(1.)
         SG=ABS(SIN(ITER*PI/2.)) 

         IF(SG.EQ.1.)THEN

             WRITE(*,*)ITER, ENORM

         ENDIF
	   write(*,*)iter,enorm
         IF((ITER.LT.ITR).OR.(ENORM.GT.ERR))GO TO 1

*..................PRESS. COEF.         

         DO I=1,N
            DO J=1,JM

               CP(I,J)=2*P(I,J)
            
            ENDDO
         ENDDO

*..................DRAG COEF.

         DTETA=2*PI/N
         DF=0.0
         
         DO I=1,N

            TETA=(I-1)*DTETA
            DF=DF+(PS(I)*COS(TETA)+2*VOR(I,1)*SIN(TETA)/REN)*DTETA
         
         ENDDO

         CD=DF/1.

         CALL OUTP

         END

*....................................................
*                                                   :
*       CELL GENERATION         :
*                                                   :
*...................................................:

         SUBROUTINE GRIDGEN

         PARAMETER(NI=100,NJ=100)

         REAL*4  X(NI,NJ),Y(NI,NJ),A(NI,NJ)
         REAL*4  R(NJ),ETA(NJ)

         COMMON/B0/X,Y,A
         COMMON/B3/M,JM,N
         COMMON/B7/R
         COMMON/B8/BETA,RADI

*        ALREF       REFRENCE LENGTH
         
         ALREF=R(1)*2
         R(1)=R(1)/ALREF
         BT=(BETA+1)/(BETA-1)
         R(M)=RADI*R(1)
         DR=R(M)-R(1)
         D=.5*LOG(1/BT)
         PI=4.*ATAN(1.)
         TET=2.*PI/FLOAT(N)

         DO J=1,M

            P=D*FLOAT(M-J)/FLOAT(M-1)
            ETA(J)=TANH(P)/TANH(D)
            PR=BT**ETA(J)
            R(J)=R(1)+DR*(1-BETA*(PR-1)/(PR+1))

         ENDDO

         DO J=1,M
            DO I=1,N+1

               TETA=PI-(I-1)*TET

               X(I,J)=R(J)*COS(TETA)
               Y(I,J)=R(J)*SIN(TETA)

            ENDDO
         ENDDO

         OPEN(7,FILE='GRID',STATUS='UNKNOWN')

         WRITE(7,*)'         X               Y     '

         DO J=1,JM
            DO I=1,N

               A(I,J)=.5*((X(I,J)-X(I+1,J+1))*(Y(I+1,J)-Y(I,J+1))
     &                   -(Y(I,J)-Y(I+1,J+1))*(X(I+1,J)-X(I,J+1)))

            ENDDO
         ENDDO

         DO J=1,JM
            DO I=1,N

               WRITE(7, 2 )X(I,J),Y(I,J),A(I,J)
 2             FORMAT(2X,3(6X,F12.6))

            ENDDO
         ENDDO

         CLOSE(7)

         RETURN
         END

*...................................................
*                                                  :
*   MAKING INITIAL VALUES :
*                                                  :
*..................................................:

         SUBROUTINE INITIALVAL

         PARAMETER(NI=100,NJ=100)

         REAL*4   U (NI,NJ),V (NI,NJ),P (NI,NJ)
         REAL*4   UN(NI,NJ),VN(NI,NJ),PN(NI,NJ)

         COMMON/B1/U,V,P
         COMMON/B2/U0,V0,P0
         COMMON/B3/M,JM,N
         COMMON/B9/UN,VN,PN
         
*..................INITIAL GUESS

         DO J=1,JM+1
            DO I=1,N

               U(I,J) = U0
               V(I,J) = V0
               P(I,J) = P0

               UN(I,J)= U(I,J)
               VN(I,J)= V(I,J)
               PN(I,J)= P(I,J)

            ENDDO
         ENDDO

         RETURN
         END

*...................................................
*                                                  :
*          CFL TEST                       :
*                                                  :
*..................................................:

         SUBROUTINE CFLTEST

         PARAMETER(NI=100,NJ=100)

         REAL*4     U(NI,NJ),V(NI,NJ),P(NI,NJ)
         REAL*4     X(NI,NJ),Y(NI,NJ),A(NI,NJ)
         REAL*4     C(NI,NJ),D(NI,NJ)

         COMMON/B0/X,Y,A
         COMMON/B1/U,V,P
         COMMON/B3/M,JM,N
         COMMON/B4/DT,CFL
         COMMON/B6/ITER
         COMMON/B11/ACOM

         IF(ITER.GT.1)GO TO 1

*        TOP,RIGHT,LEFT & BOTTOM DELTA Y & DELTA X (FOR PRIMARY CELLS)

         DO J=1,JM
            DO I=1,N

               DYT=Y(I+1,J+1)-Y(I,J+1)
               DYR=Y(I+1,J)-Y(I+1,J+1)
               DYB=Y(I,J)-Y(I+1,J)
               DYL=Y(I,J+1)-Y(I,J)

               DXT=X(I+1,J+1)-X(I,J+1)
               DXR=X(I+1,J)-X(I+1,J+1)
               DXB=X(I,J)-X(I+1,J)
               DXL=X(I,J+1)-X(I,J)

               DLT=SQRT(DYT**2+DXT**2)
               DLL=SQRT(DYL**2+DXL**2)
               DLB=SQRT(DYB**2+DXB**2)
               DLR=SQRT(DYR**2+DXR**2)

               D(I,J)=MIN(DLT,DLL,DLB,DLR)

            ENDDO
         ENDDO

         DXMIN = D(1,1)

         DO J=1,JM
            DO I=1,N

               IF(DXMIN.GT.D(I,J))DXMIN=D(I,J)

            ENDDO
         ENDDO

1        DO J=1,JM
            DO I=1,N

               C(I,J)=SQRT(U(I,J)**2+V(I,J)**2)+SQRT(ACOM)
                                                          
            ENDDO
         ENDDO

         CMAX =C (1,1)
         
         DO J=1,JM
            DO I=1,N

               IF (CMAX.LT.C(I,J))CMAX=C(I,J)
               
            ENDDO
         ENDDO

         DT=CFL*DXMIN/CMAX

         RETURN
         END

*...........................................................
*                                                          :
*      SOLUTION UPDATE                 :
*    4TH ORDER RUNGE-KUTTA     :
*..........................................................:

         SUBROUTINE RK4

         PARAMETER(NI=100,NJ=100)

         REAL*4     QPA(NI,NJ),QUA(NI,NJ),QVA(NI,NJ)
         REAL*4     QPB(NI,NJ),QUB(NI,NJ),QVB(NI,NJ)
         REAL*4     QPC(NI,NJ),QUC(NI,NJ),QVC(NI,NJ)
         REAL*4     QPD(NI,NJ),QUD(NI,NJ),QVD(NI,NJ)
         REAL*4     UN (NI,NJ),VN (NI,NJ),PN (NI,NJ)
         REAL*4     U  (NI,NJ),V  (NI,NJ),P  (NI,NJ)
         REAL*4     DIF(NI,NJ)
         
         COMMON/B1/U,V,P
         COMMON/B2/UND,VND,PND
         COMMON/B3/M,JM,N
         COMMON/B4/DT,CFL
         COMMON/B5/DIF
         COMMON/B9/UN,VN,PN
         
         DO J=1,JM
            DO I=1,N

               CALL FLUX(I,J,QPA(I,J),QUA(I,J),QVA(I,J))
               CALL DISS(I,J,DISU,DISV,DISP)

               U(I,J) =UN(I,J) - ( QUA(I,J)+ DISU ) * DT/2.
               V(I,J) =VN(I,J) - ( QVA(I,J)+ DISV ) * DT/2.
               P(I,J) =PN(I,J) - ( QPA(I,J)+ DISP ) * DT/2.

            ENDDO
         ENDDO

         DO J=1,JM
            DO I=1,N

               CALL FLUX(I,J,QPB(I,J),QUB(I,J),QVB(I,J))
               CALL DISS(I,J,DISU,DISV ,DISP)

               U(I,J) =UN(I,J) - ( QUB(I,J)+ DISU ) * DT/2.
               V(I,J) =VN(I,J) - ( QVB(I,J)+ DISV ) * DT/2.
               P(I,J) =PN(I,J) - ( QPB(I,J)+ DISP ) * DT/2.


            ENDDO
         ENDDO

         DO J=1,JM
            DO I=1,N

               CALL FLUX(I,J,QPC(I,J),QUC(I,J),QVC(I,J))
               CALL DISS(I,J,DISU,DISV ,DISP)

               U(I,J) =UN(I,J) - ( QUC(I,J)+ DISU ) * DT/2.
               V(I,J) =VN(I,J) - ( QVC(I,J)+ DISV ) * DT/2.
               P(I,J) =PN(I,J) - ( QPC(I,J)+ DISP ) * DT/2.

            ENDDO
         ENDDO

         DO J=1,JM
            DO I=1,N

               CALL FLUX(I,J,QPD(I,J),QUD(I,J),QVD(I,J))
               CALL DISS(I,J,DISU,DISV ,DISP)

               P(I,J)=PN(I,J)-(QPA(I,J)+2.*QPB(I,J)+2.*
     &                QPC(I,J)+QPD(I,J))*DT/6.-DT*DISP/6.

               U(I,J)=UN(I,J)-(QUA(I,J)+2.*QUB(I,J)+2.*
     &                QUC(I,J)+QUD(I,J))*DT/6.-DT*DISU/6.

               V(I,J)=VN(I,J)-(QVA(I,J)+2.*QVB(I,J)+2.*
     &                QVC(I,J)+QVD(I,J))*DT/6.-DT*DISV/6.


            ENDDO
         ENDDO

         CALL FARF
         
         DO J=1,JM
            DO I=1,N

               DIF(I,J)=(P(I,J)-PN(I,J))**2

               PN(I,J)=P(I,J)
               UN(I,J)=U(I,J)
               VN(I,J)=V(I,J)
               
            ENDDO
         ENDDO

         DO J=1,JM+1

            P(N+1,J)=P(1,J)
            U(N+1,J)=U(1,J)
            V(N+1,J)=V(1,J)
            
         ENDDO

         RETURN
         END
*.........................................................
*                                                        :
*                CELL FLUXES                :
*                                                        :
*........................................................:

         SUBROUTINE FLUX(I,J,QP,QU,QV)

         PARAMETER(NI=100,NJ=100)

         REAL*4    U(NI,NJ),V(NI,NJ),P(NI,NJ)
         REAL*4    X(NI,NJ),Y(NI,NJ),A(NI,NJ)
         REAL*4    VOR(NI,NJ),PS(NI)
         
         COMMON/B0/X,Y,A
         COMMON/B1/U,V,P
         COMMON/B3/M,JM,N
         COMMON/B10/REN
         COMMON/B11/ACOM
         COMMON/B13/VOR
         COMMON/B19/PS

*        CONTINUITY OF THE INTIAL & FINAL CELLS BY THE TWO (IF) STATMENT

         IF(I.EQ.1)THEN
            II=N+1
         ELSE
            II=I
         ENDIF

         IF(I.EQ.N)THEN
            IN=0
         ELSE
            IN=I
         ENDIF

*        TOP,RIGHT,BOTTOM & LEFT DELTA Y OF THE PRIMARY CELLS
         
         DYT=Y(IN+1,J+1)-Y(I,J+1)
         DYR=Y(IN+1,J)-Y(IN+1,J+1)
         DYB=Y(I,J)-Y(IN+1,J)
         DYL=Y(I,J+1)-Y(I,J)

*        TOP,RIGHT,BOTTOM & LEFT DELTA X OF THE PRIMARY CELLS

         DXT=X(IN+1,J+1)-X(I,J+1)
         DXR=X(IN+1,J)-X(IN+1,J+1)
         DXB=X(I,J)-X(IN+1,J)
         DXL=X(I,J+1)-X(I,J)

*        CENTER POINT OF THE TOP,RIGHT,BOTTOM & LEFT SECONDARY CELLS 

         YC=.25*(Y(I,J)+Y(I,J+1)+Y(IN+1,J+1)+Y(IN+1,J))
         YT=.25*(Y(I,J+1)+Y(IN+1,J+1)+Y(IN+1,J+2)+Y(I,J+2))
         YR=.25*(Y(IN+1,J)+Y(IN+1,J+1)+Y(IN+2,J)+Y(IN+2,J+1))
         YL=.25*(Y(I,J)+Y(II-1,J)+Y(I,J+1)+Y(II-1,J+1))
         YB=.25*(Y(I,J-1)+Y(IN+1,J-1)+Y(I,J)+Y(IN+1,J))

         XC=.25*(X(I,J)+X(I,J+1)+X(IN+1,J+1)+X(IN+1,J))
         XT=.25*(X(I,J+1)+X(IN+1,J+1)+X(IN+1,J+2)+X(I,J+2))
         XR=.25*(X(IN+1,J)+X(IN+1,J+1)+X(IN+2,J)+X(IN+2,J+1))
         XL=.25*(X(I,J)+X(II-1,J)+X(I,J+1)+X(II-1,J+1))
         XB=.25*(X(I,J-1)+X(IN+1,J-1)+X(I,J)+X(IN+1,J))

*..................CONVECTIVE TERMS

*..................TOP F FLUX

         FT1=.5*(U(I,J+1)+U(I,J))*DYT*ACOM
         FT2=.5*(U(I,J+1)*U(I,J+1)+P(I,J+1)+U(I,J)*U(I,J)+P(I,J))*DYT
         FT3=.5*(U(I,J+1)*V(I,J+1)+U(I,J)*V(I,J))*DYT

*..................BOTTOM F FLUX

         FB1=.5*(U(I,J)+U(I,J-1))*DYB*ACOM
         FB2=.5*(U(I,J)*U(I,J)+P(I,J)+U(I,J-1)*U(I,J-1)+P(I,J-1))*DYB
         FB3=.5*(U(I,J)*V(I,J)+U(I,J-1)*V(I,J-1))*DYB

*..................RIGHT F FLUX

         FR1=.5*(U(IN+1,J)+U(I,J))*DYR*ACOM
         FR2=.5*(U(IN+1,J)*U(IN+1,J)+P(IN+1,J)+U(I,J)*U(I,J)+P(I,J))*DYR
         FR3=.5*(U(IN+1,J)*V(IN+1,J)+U(I,J)*V(I,J))*DYR

*..................LEFT F FLUX

         FL1=.5*(U(I,J)+U(II-1,J))*DYL*ACOM
         FL2=.5*(U(I,J)*U(I,J)+P(I,J)+U(II-1,J)*U(II-1,J)+P(II-1,J))*DYL
         FL3=.5*(U(I,J)*V(I,J)+U(II-1,J)*V(II-1,J))*DYL

*..................TOP G FLUX

         GT1=.5*(V(I,J+1)+V(I,J))*DXT*ACOM
         GT2=.5*(U(I,J+1)*V(I,J+1)+U(I,J)*V(I,J))*DXT
         GT3=.5*(V(I,J+1)*V(I,J+1)+P(I,J+1)+V(I,J)*V(I,J)+P(I,J))*DXT

*..................BOTTOM G FLUX

         GB1=.5*(V(I,J)+V(I,J-1))*DXB*ACOM
         GB2=.5*(U(I,J)*V(I,J)+U(I,J-1)*V(I,J-1))*DXB
         GB3=.5*(V(I,J)*V(I,J)+P(I,J)+V(I,J-1)*V(I,J-1)+P(I,J-1))*DXB

*..................RIGHT G FLUX

         GR1=.5*(V(IN+1,J)+V(I,J))*DXR*ACOM
         GR2=.5*(U(IN+1,J)*V(IN+1,J)+U(I,J)*V(I,J))*DXR
         GR3=.5*(V(IN+1,J)*V(IN+1,J)+P(IN+1,J)+V(I,J)*V(I,J)+P(I,J))*DXR

*..................LEFT G  FLUX

         GL1=.5*(V(I,J)+V(II-1,J))*DXL*ACOM
         GL2=.5*(U(I,J)*V(I,J)+U(II-1,J)*V(II-1,J))*DXL
         GL3=.5*(V(I,J)*V(I,J)+P(I,J)+V(II-1,J)*V(II-1,J)+P(II-1,J))*DXL

*        TOP,RIGHT,BOTTOM & LEFT AREA OF THE SECONDARY CELLS

         AT=.5*((YT-YC)*DXT-(XT-XC)*DYT)
         AB=.5*((YB-YC)*DXB-(XB-XC)*DYB)
         AL=.5*((XC-XL)*DYL-(YC-YL)*DXL)
         AR=.5*((XC-XR)*DYR-(YC-YR)*DXR)
          
*        VELOCITY OF THE PRIMARY CELLS NODE(CORNER)  

         UNE=.25*(U(I,J)+U(IN+1,J)+U(I,J+1)+U(IN+1,J+1))
         USE=.25*(U(I,J-1)+U(IN+1,J-1)+U(I,J)+U(IN+1,J))
         USW=.25*(U(II-1,J-1)+U(I,J-1)+U(II-1,J)+U(I,J))
         UNW=.25*(U(II-1,J)+U(I,J)+U(II-1,J+1)+U(I,J+1))

         VNE=.25*(V(I,J)+V(IN+1,J)+V(I,J+1)+V(IN+1,J+1))
         VSE=.25*(V(I,J-1)+V(IN+1,J-1)+V(I,J)+V(IN+1,J))
         VSW=.25*(V(II-1,J-1)+V(I,J-1)+V(II-1,J)+V(I,J))
         VNW=.25*(V(II-1,J)+V(I,J)+V(II-1,J+1)+V(I,J+1))

*        THE SECONDARY CELLS DELTA X & DELTA Y 

         DYTAB=Y(I,J+1)-YT
         DYTBC=YT-Y(IN+1,J+1)
         DYTCD=Y(IN+1,J+1)-YC
         DYTDA=YC-Y(I,J+1)
         
         DXTAB=X(I,J+1)-XT
         DXTBC=XT-X(IN+1,J+1)
         DXTCD=X(IN+1,J+1)-XC
         DXTDA=XC-X(I,J+1)

         DYBAB=Y(I,J)-YC
         DYBBC=YC-Y(IN+1,J)
         DYBCD=Y(IN+1,J)-YB
         DYBDA=YB-Y(I,J)

         DXBAB=X(I,J)-XC
         DXBBC=XC-X(IN+1,J)
         DXBCD=X(IN+1,J)-XB
         DXBDA=XB-X(I,J)

         DYLAB=YL-Y(I,J+1)
         DYLBC=Y(I,J+1)-YC
         DYLCD=YC-Y(I,J)
         DYLDA=Y(I,J)-YL

         DXLAB=XL-X(I,J+1)
         DXLBC=X(I,J+1)-XC
         DXLCD=XC-X(I,J)
         DXLDA=X(I,J)-XL

         DYRAB=YC-Y(IN+1,J+1)
         DYRBC=Y(IN+1,J+1)-YR
         DYRCD=YR-Y(IN+1,J)
         DYRDA=Y(IN+1,J)-YC

         DXRAB=XC-X(IN+1,J+1)
         DXRBC=X(IN+1,J+1)-XR
         DXRCD=XR-X(IN+1,J)
         DXRDA=X(IN+1,J)-XC

*..................VISCOUS TERMS


*..................R - TOP

         RT1=0.0

         RT2=(.5*(UNE+U(I,J+1))*DYTAB
     &       +.5*(UNW+U(I,J+1))*DYTBC
     &       +.5*(UNW+U(I,J))  *DYTCD
     &       +.5*(UNE+U(I,J))  *DYTDA)*DYT/AT

         RT3=(.5*(VNE+V(I,J+1))*DYTAB
     &       +.5*(VNW+V(I,J+1))*DYTBC
     &       +.5*(VNW+V(I,J))  *DYTCD
     &       +.5*(VNE+V(I,J))  *DYTDA)*DYT/AT

*..................R - RIGHT
        
         RR1=0.0

         RR2=(.5*(UNE+U(I,J))   *DYRAB 
     &       +.5*(UNE+U(IN+1,J))*DYRBC
     &       +.5*(USE+U(IN+1,J))*DYRCD
     &       +.5*(USE+U(I,J))   *DYRDA)*DYR/AR

         RR3=(.5*(VNE+V(I,J))   *DYRAB
     &       +.5*(VNE+V(IN+1,J))*DYRBC
     &       +.5*(VSE+V(IN+1,J))*DYRCD
     &       +.5*(VSE+V(I,J))   *DYRDA)*DYR/AR

*..................R - LEFT
        
         RL1=0.0

         RL2=(.5*(UNW+U(II-1,J))*DYLAB 
     &       +.5*(UNW+U(I,J))   *DYLBC
     &       +.5*(USW+U(I,J))   *DYLCD
     &       +.5*(USW+U(II-1,J))*DYLDA)*DYL/AL

         RL3=(.5*(VNW+V(II-1,J))*DYLAB
     &       +.5*(VNW+V(I,J))   *DYLBC
     &       +.5*(VSW+V(I,J))   *DYLCD
     &       +.5*(VSW+V(II-1,J))*DYLDA)*DYL/AL

*..................R - BOTTOM
        
         RB1=0.0

         RB2=(.5*(USE+U(I,J))  *DYBAB
     &       +.5*(USW+U(I,J))  *DYBBC
     &       +.5*(USW+U(I,J-1))*DYBCD
     &       +.5*(USE+U(I,J-1))*DYBDA)*DYB/AB

         RB3=(.5*(VSE+V(I,J))  *DYBAB
     &       +.5*(VSW+V(I,J))  *DYBBC
     &       +.5*(VSW+V(I,J-1))*DYBCD
     &       +.5*(VSE+V(I,J-1))*DYBDA)*DYB/AB

*..................S - TOP

         ST1=0.0

         ST2=-(.5*(UNE+U(I,J+1))*DXTAB
     &        +.5*(UNW+U(I,J+1))*DXTBC
     &        +.5*(UNW+U(I,J))  *DXTCD
     &        +.5*(UNE+U(I,J))  *DXTDA)*DXT/AT

         ST3=-(.5*(VNE+V(I,J+1))*DXTAB
     &        +.5*(VNW+V(I,J+1))*DXTBC
     &        +.5*(VNW+V(I,J))  *DXTCD
     &        +.5*(VNE+V(I,J))  *DXTDA)*DXT/AT

*..................S - RIGHT
        
         SR1=0.0

         SR2=-(.5*(UNE+U(I,J))   *DXRAB 
     &        +.5*(UNE+U(IN+1,J))*DXRBC
     &        +.5*(USE+U(IN+1,J))*DXRCD
     &        +.5*(USE+U(I,J))   *DXRDA)*DXR/AR

         SR3=-(.5*(VNE+V(I,J))   *DXRAB
     &        +.5*(VNE+V(IN+1,J))*DXRBC
     &        +.5*(VSE+V(IN+1,J))*DXRCD
     &        +.5*(VSE+V(I,J))   *DXRDA)*DXR/AR

*..................S - LEFT
        
         SL1=0.0

         SL2=-(.5*(UNW+U(II-1,J))*DXLAB 
     &        +.5*(UNW+U(I,J))   *DXLBC
     &        +.5*(USW+U(I,J))   *DXLCD
     &        +.5*(USW+U(II-1,J))*DXLDA)*DXL/AL

         SL3=-(.5*(VNW+V(II-1,J))*DXLAB
     &        +.5*(VNW+V(I,J))   *DXLBC
     &        +.5*(VSW+V(I,J))   *DXLCD
     &        +.5*(VSW+V(II-1,J))*DXLDA)*DXL/AL

*..................S - BOTTOM
        
         SB1=0.0

         SB2=-(.5*(USE+U(I,J))  *DXBAB
     &        +.5*(USW+U(I,J))  *DXBBC
     &        +.5*(USW+U(I,J-1))*DXBCD
     &        +.5*(USE+U(I,J-1))*DXBDA)*DXB/AB

         SB3=-(.5*(VSE+V(I,J))  *DXBAB
     &        +.5*(VSW+V(I,J))  *DXBBC
     &        +.5*(VSW+V(I,J-1))*DXBCD
     &        +.5*(VSE+V(I,J-1))*DXBDA)*DXB/AB

*..................STREAM FUNCTION (VORTICITY)

         VT=.5*(V(I,J)+V(I,J+1)) 
         VB=.5*(V(I,J)+V(I,J-1)) 
         VL=.5*(V(I,J)+V(II-1,J))
         VR=.5*(V(I,J)+V(IN+1,J))

         UT=.5*(U(I,J)+U(I,J+1)) 
         UB=.5*(U(I,J)+U(I,J-1)) 
         UL=.5*(U(I,J)+U(II-1,J))
         UR=.5*(U(I,J)+U(IN+1,J))


         IF(J.EQ.1)THEN

            VB=0.0
            UB=0.0

         ENDIF

         VDY=VT*DYT+VR*DYR+VB*DYB+VL*DYL
         UDX=UT*DXT+UR*DXR+UB*DXB+UL*DXL
         
         VOR(I,J)=(VDY+UDX)/A(I,J)

*..................SOLID BOUNDARY

         IF(J.EQ.1) THEN

            AB=.5*(YC*(X(IN+1,J)-X(I,J))
     &            +Y(I,J)*(XC-X(IN+1,J))
     &            +Y(IN+1,J)*(X(I,J)-XC))


            RB1=0.0
            RB2=.5*U(I,J)*DYB*DYB/AB   
            RB3=.5*V(I,J)*DYB*DYB/AB

            SB1=0.0
            SB2=-.5*U(I,J)*DXB*DXB/AB
            SB3=-.5*V(I,J)*DXB*DXB/AB

C           PS(I)=((SR3+SL3+ST3+SB3)/REN
C     &            -.5*(P(II-1,J)+P(I,J))*DXL
C     &            -.5*(P(I,J)+P(I,J+1)) *DXT
C     &            -.5*(P(IN+1,J)+P(I,J))*DXR)/DXB

            PS(I)=1.5*P(I,J)-.5*P(I,J+1)

            FB1=0.0
            FB2=PS(I)*DYB
            FB3=0.0

            GB1=0.0
            GB2=0.0
            GB3=PS(I)*DXB

         ENDIF
  
         F1=FT1+FB1+FR1+FL1
         F2=FT2+FB2+FR2+FL2
         F3=FT3+FB3+FR3+FL3

         G1=GT1+GB1+GR1+GL1
         G2=GT2+GB2+GR2+GL2
         G3=GT3+GB3+GR3+GL3

         R1=RT1+RB1+RR1+RL1
         R2=RT2+RB2+RR2+RL2
         R3=RT3+RB3+RR3+RL3

         S1=ST1+SB1+SR1+SL1
         S2=ST2+SB2+SR2+SL2
         S3=ST3+SB3+SR3+SL3
       
         QP=((R1-S1)/REN-(F1-G1))/A(I,J)
         QU=((R2-S2)/REN-(F2-G2))/A(I,J)
         QV=((R3-S3)/REN-(F3-G3))/A(I,J)

         RETURN
         END

!...........................................................
!*                                                        *
!*   ARTIFITIAL DISSIPATION        *
!*........................................................*

         SUBROUTINE DISS(I,J,DISU,DISV,DISP)

         PARAMETER(NI=100,NJ=100)

         REAL*4   U(NI,NJ),V(NI,NJ),P(NI,NJ)
         REAL*4   X(NI,NJ),Y(NI,NJ),A(NI,NJ)

         COMMON/B0/X,Y,A
         COMMON/B1/U,V,P
         COMMON/B3/M,JM,N
         COMMON/B14/EU,EV,EP

         IF(J.LE.2)THEN
            II=N+1
         ELSE
            II=I
         ENDIF

         IF(I.EQ.N)THEN
            IN=0
         ELSE
            IN=I
         ENDIF

         IF((J.LT.2).OR.(J.EQ.JM))THEN

            DISU=0.0
            DISV=0.0
            DISP=0.0

         ELSE

            DISUI=EU*(U(IN+2,J)-4*U(IN+1,J)+6*U(I,J)
     &                         -4*U(II-1,J)+U(II-2,J))
            DISUJ=EU*(U(I,J+2) -4*U(I,J+1)+6*U(I,J)
     &                         -4*U(I,J-1)+U(I,J-2))

            DISVI=EV*(V(I+2,J)-4*V(IN+1,J)+6*V(I,J)
     &                        -4*V(II-1,J)+V(II-2,J))
            DISVJ=EV*(V(I,J+2)-4*V(I,J+1)+6*V(I,J)
     &                        -4*V(I,J-1)+V(I,J-2))

            DISPI=EP*(P(I+2,J)-4*P(IN+1,J)+6*P(I,J)
     &                        -4*P(II-1,J)+P(II-2,J))
            DISPJ=EP*(P(I,J+2)-4*P(I,J+1)+6*P(I,J)
     &                        -4*P(I,J-1)+P(I,J-2))

            DISU=(DISUI+DISUJ)/A(I,J)
            DISV=(DISVI+DISVJ)/A(I,J)
            DISP=(DISPI+DISPJ)/A(I,J)

         ENDIF

         RETURN
         END
     
!...................................................
!*                                                :
!*        ACCURACY TEST         :
!*                                                :
!*................................................:

         SUBROUTINE ACCUR(ENORM)

         PARAMETER(NI=100,NJ=100)

         REAL*4  DIF(NI,NJ)

         COMMON/B3/M,JM,N
         COMMON/B5/DIF

         SUM1=0.0

         DO J=1,JM
            DO I=1,N

               SUM1 =SUM1 +DIF(I,J)

            ENDDO
         ENDDO

         ENORM= SQRT(SUM1)/FLOAT(JM*N)

         RETURN
         END

!*.........................................................
!*                                                        :
!*               FAR-FIELD                      :
!*                                                        :
!*........................................................:

         SUBROUTINE FARF

         PARAMETER(NI=100,NJ=100)

         REAL*4    U(NI,NJ),V(NI,NJ),P(NI,NJ)

         COMMON/B1/U,V,P
         COMMON/B2/U0,V0,P0
         COMMON/B3/M,JM,N

!*..................INFLOW

         DO I=N/4+1,3*N/4

            U(I,JM+1)=U0
            V(I,JM+1)=V0
            P(I,JM+1)=1.5*P(I,JM)-.5*P(I,JM-1)

            ENDDO

!*..................OUTFLOW

         DO I=1,N/4

            U(I,JM+1)=U0
            V(I,JM+1)=V0
            P(I,JM+1)=P0
            
         ENDDO

         DO I=3*N/4+1,N

            U(I,JM+1)=U0
            V(I,JM+1)=V0
            P(I,JM+1)=P0
            
         ENDDO

         RETURN
         END
         
!.........................................................
!*                                                      :
!*            OUTPUTS                         :
!*                                                      :
!*......................................................:

         SUBROUTINE OUTP

         PARAMETER(NI=100,NJ=100)

         REAL*4    U(NI,NJ),V(NI,NJ),P(NI,NJ)
         REAL*4    X(NI,NJ),Y(NI,NJ),A(NI,NJ)
         REAL*4    VOR(NI,NJ),UT(NI,NJ),R(NJ)
         REAL*4    CP(NI,NJ)

         COMMON/B0/X,Y,A
         COMMON/B1/U,V,P
         COMMON/B3/M,JM,N
         COMMON/B4/DT,CFL
         COMMON/B6/ITER
         COMMON/B7/R
         COMMON/B13/VOR
         COMMON/B17/CP,CD

         TET=360/FLOAT(N)

         OPEN(UNIT=1 ,FILE='OUTP1.DAT',STATUS='UNKNOWN')
         OPEN(UNIT=2 ,FILE='OUTP2.DAT',STATUS='UNKNOWN')
         OPEN(UNIT=3 ,FILE=  'VOR.DAT',STATUS='UNKNOWN')
         OPEN(UNIT=4 ,FILE=  'PRS.DAT',STATUS='UNKNOWN')
         OPEN(UNIT=5 ,FILE=   'CP.DAT',STATUS='UNKNOWN')  
         OPEN(UNIT=6 ,FILE=    'U.DAT',STATUS='UNKNOWN')

         WRITE(1,1)N+1,M-1
 1       FORMAT(2X,'TITLE="OUTPUT"',/2X,'VARIABLES= X, Y, U, V, P,VOR',/2X,
     &          'ZONE T="ZONE-ONE"',5X,'I=',I3,3X,'J=',I3,3X,'F=POINT')

         WRITE(2,2 )ITER,DT,CD
 2       FORMAT(/,4X,'ITER=',I4,/,4X,'DT=',F7.3,/,4X,'D. COEF=',F7.3)
         
         DO J=1,JM+1
            DO I=1,N+1

               WRITE(1,5)X(I,J),Y(I,J),U(I,J),V(I,J),P(I,J),VOR(I,J)
            
            ENDDO
         ENDDO

         DO I=1,N

            TETA=(I-1)*TET

            WRITE(3,8)TETA,VOR(I,1)
            WRITE(4,8)TETA,  P(I,1)
            WRITE(5,8)TETA, CP(I,1)
            WRITE(6,8)TETA,  U(I,1)

         ENDDO

 5       FORMAT(5(3X,F10.3))         
 6       FORMAT(2(1X,F10.3))
 7       FORMAT(2X,3(2X,F10.3))
 8       FORMAT(2X,2(2X,F10.3))

         RETURN
         END



