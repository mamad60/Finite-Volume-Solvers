function [ QP,QU,QV ] = Flux( i,j )
%Return right hand side for each eqs.
global  X Y;
global XC YC;
global U V P;
global Beta;
global A;
global Ren;
global m n;

%Note:i,j in soloution shows cell centers while in X,Y refer to grid
%points

%*********************************************

%Begin Geometrical Values Calculations
%--------------------------------------------------------------
%Left, Bottom ,Right & Top DeltaY of Primary cells
DYT=Y(i,j)-Y(i+1,j);   %Deta Y of Top Face
DYL=Y(i,j+1)-Y(i,j);   %Deta Y of Left Face
DYB=Y(i+1,j+1)-Y(i,j+1); %Deta Y of Bottom Face
DYR=Y(i+1,j)-Y(i+1,j+1); %Deta Y of Right Face

%Left, Bottom ,Right & Top DeltaX of Primary cells
DXT=X(i,j)-X(i+1,j);  %Deta X of Top Face
DXL=X(i,j+1)-X(i,j);  %Deta X of Left Face
DXB=X(i+1,j+1)-X(i,j+1); %Deta X of Bottom Face
DXR=X(i+1,j)-X(i+1,j+1); %Deta X of Right Face 

%-------------------------------------------------------------


%Y of Center points of Left, Bottom ,Right & Top of Primary cells


    YL=0.25*( Y(i-1,j)+Y(i,j)+Y(i-1,j+1)+Y(i,j+1) );  %YC(i-1,j)
    YB=0.25*( Y(i,j+1)+Y(i+1,j+1)+Y(i,j+2)+Y(i+1,j+2) ); %YC(i,j+1)
    YR=0.25*( Y(i+1,j)+Y(i+2,j)+Y(i+1,j+1)+Y(i+2,j+1) ); %YC(i+1,j)
    YT=0.25*( Y(i,j-1)+Y(i+1,j-1)+Y(i,j)+Y(i+1,j) ); %YC(i,j+1)


%X of Center points of Left, Bottom ,Right & Top of Primary cells

    XL=0.25*( X(i-1,j)+X(i,j)+X(i-1,j+1)+X(i,j+1) ); %XC(i-1,j)
    XB=0.25*( X(i,j+1)+X(i+1,j+1)+X(i,j+2)+X(i+1,j+2) ); %XC(i,j+1)
    XR=0.25*( X(i+1,j)+X(i+2,j)+X(i+1,j+1)+X(i+2,j+1) ); %XC(i+1,j)
    XT=0.25*( X(i,j-1)+X(i+1,j-1)+X(i,j)+X(i+1,j) );  %XC(i,j+1)


%--------------------------------------------------------------------------


%Secondary cells DeltaY & DeltaX-----------------------------------------
    %DeltaY
            %L Secondary Cell      
            DYLDA=YC(i,j)-Y(i,j+1);
            DYLCD=Y(i,j)-YC(i,j);
            DYLBC=YL-Y(i,j);
            DYLAB=Y(i,j+1)-YL;
            
            %B Secondary Cell
            DYBDA=Y(i+1,j+1)-YB;
            DYBAB=YB-Y(i,j+1);
            DYBCD=YC(i,j)-Y(i+1,j+1);
            DYBBC=Y(i,j+1)-YC(i,j);
            
            %R Secondary Cell     
            DYRDA=YR-Y(i+1,j+1);
            DYRCD=Y(i+1,j)-YR;
            DYRBC=YC(i,j)-Y(i+1,j);
            DYRAB=Y(i+1,j+1)-YC(i,j);
            
            
            %T  Secondary Cell     
            DYTDA=Y(i+1,j)-YC(i,j);
            DYTCD=YT-Y(i+1,j);
            DYTBC=Y(i,j)-YT;
            DYTAB=YC(i,j)-Y(i,j);
            
    %DeltaX
            %L  Secondary Cell      
            DXLDA=XC(i,j)-X(i,j+1);
            DXLCD=X(i,j)-XC(i,j);
            DXLBC=XL-X(i,j);
            DXLAB=X(i,j+1)-XL;
            
            
            %B Secondary Cell 
            DXBDA=X(i+1,j+1)-XB;
            DXBAB=XB-X(i,j+1);
            DXBCD=XC(i,j)-X(i+1,j+1);
            DXBBC=X(i,j+1)-XC(i,j);
            
            %R Secondary Cell      
            DXRDA=XR-X(i+1,j+1);
            DXRCD=X(i+1,j)-XR;
            DXRBC=XC(i,j)-X(i+1,j);
            DXRAB=X(i+1,j+1)-XC(i,j);
            
            
            %T Secondary Cell      
            DXTDA=X(i+1,j)-XC(i,j);
            DXTCD=XT-X(i+1,j);
            DXTBC=X(i,j)-XT;
            DXTAB=XC(i,j)-X(i,j);

            
 %------------------------------------------------------------            
 %Area of Secondary cells-------------------------------------------------
 
              AL=0.5*( (XC(i,j)-XL)*DYL-(YC(i,j)-YL)*DXL ); %Area of L Secondary cell
              AB=0.5*( (YB-YC(i,j))*DXB-(XB-XC(i,j))*DYB ); %Area of B Secondary cell
              AR=0.5*( (XC(i,j)-XR)*DYR-(YC(i,j)-YR)*DXR ); %Area of R Secondary cell
              AT=0.5*( (YT-YC(i,j))*DXT-(XT-XC(i,j))*DYT ); %Area of T Secondary cell
 
 %End of Geometrical Values Calculation
%*********************************************
         
%--------------------------------------------------------------           
  %Velocity of the Primary cell nodes(Corners)-- Required for Second Order Terms Calculations
  %Interpolated From Cell Center Values
            %u Values at Corners
         Unw=0.25*( U(i-1,j-1)+U(i,j-1)+U(i-1,j)+U(i,j) );
         Une=0.25*(  U(i,j-1)+U(i+1,j-1)+U(i,j)+U(i+1,j) );
         Usw=0.25*( U(i-1,j)+U(i,j)+U(i-1,j+1)+U(i,j+1) );
         Use=0.25*( U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1) );
            %v Values at Corners
         Vnw=0.25*( V(i-1,j-1)+V(i,j-1)+V(i-1,j)+V(i,j) );
         Vne=0.25*(  V(i,j-1)+V(i+1,j-1)+V(i,j)+V(i+1,j) );
         Vsw=0.25*( V(i-1,j)+V(i,j)+V(i-1,j+1)+V(i,j+1) );
         Vse=0.25*( V(i,j)+V(i+1,j)+V(i,j+1)+V(i+1,j+1) );

        
     
%--------------------------------------------------------------------------

%***********************************************
%Calculations for First Order Terms
%Calaculations for internal(not next to boundary cells)

%----------------------------------------------------------------
    %Convective Terms--------Calculated From Averaging(Central Difference)
    %----F Fluxes------
        %Left F flux
    FL1=0.5*( U(i-1,j)+U(i,j) )*DYL*Beta ;   % U*DY*Beta on Left Face
    FL2=0.5*( U(i-1,j)*U(i-1,j)+P(i-1,j)+U(i,j)*U(i,j)+P(i,j) )*DYL ; % (U^2+P)*DY on Left Face
    FL3=0.5*( U(i-1,j)*V(i-1,j)+U(i,j)*V(i,j) )*DYL ; %U*V on Left Face
    
        %Bottom F flux
    FB1=0.5*( U(i,j)+U(i,j+1) )*DYB*Beta ;   % U*DY*Beta on Bottom Face
    FB2=0.5*( U(i,j)*U(i,j)+P(i,j)+U(i,j+1)*U(i,j+1)+P(i,j+1) )*DYB ;  %(U^2+P)*DY on Bottom Face
    FB3=0.5*( U(i,j)*V(i,j)+U(i,j+1)*V(i,j+1) )*DYB ;  %U*V on Bottom Face
    
        %Right F flux
    FR1=0.5*( U(i,j)+U(i+1,j) )*DYR*Beta ;  % U*DY*Beta on Right Face
    FR2=0.5*( U(i,j)*U(i,j)+P(i,j)+U(i+1,j)*U(i+1,j)+P(i+1,j) )*DYR ; % (U^2+P)*DY on Right Face
    FR3=0.5*( U(i,j)*V(i,j)+U(i+1,j)*V(i+1,j) )*DYR ;  %U*V on Right Face
    
        %Top F flux
    FT1=0.5*( U(i,j)+U(i,j-1) )*DYT*Beta ; % U*DY*Beta on Top Face
    FT2=0.5*( U(i,j)*U(i,j)+P(i,j)+U(i,j-1)*U(i,j-1)+P(i,j-1) )*DYT ; % (U^2+P)*DY on Top Face
    FT3=0.5*( U(i,j)*V(i,j)+U(i,j-1)*V(i,j-1) )*DYT ;   % (U^2+P)*DY on Top Face
    
   %----F Fluxes------
        %Left G flux
    GL1=0.5*( V(i,j)+V(i-1,j) )*DXL*Beta ;  % V*DX*Beta on Left Face
    GL2=0.5*( U(i,j)*V(i,j)+U(i-1,j)*V(i-1,j) )*DXL ; % U*V*DX on Left Face
    GL3=0.5*( V(i,j)*V(i,j)+P(i,j)+V(i-1,j)*V(i-1,j)+P(i-1,j) )*DXL ; % (V^2+p)*DX on Left Face

        %Bottom G flux
    GB1=0.5*( V(i,j)+V(i,j+1) )*DXB*Beta ; % V*DX*Beta on Bottom Face
    GB2=0.5*( U(i,j)*V(i,j)+U(i,j+1)*V(i,j+1) )*DXB ; % U*V*DX on Bottom Face
    GB3=0.5*( V(i,j)*V(i,j)+P(i,j)+V(i,j+1)*V(i,j+1)+P(i,j+1) )*DXB ; % (V^2+p)*DX on Bottom Face

        %Right G flux
    GR1=0.5*( V(i,j)+V(i+1,j) )*DXR*Beta ; % V*DX*Beta on Right Face
    GR2=0.5*( U(i,j)*V(i,j)+U(i+1,j)*V(i+1,j) )*DXR ; % U*V*DX on Right Face
    GR3=0.5*( V(i,j)*V(i,j)+P(i,j)+V(i+1,j)*V(i+1,j)+P(i+1,j) )*DXR ; % (V^2+p)*DX on Right Face
    

        %Top G flux
    GT1=0.5*( V(i,j)+V(i,j-1) )*DXT*Beta ; % V*DX*Beta on Top Face
    GT2=0.5*( U(i,j)*V(i,j)+U(i,j-1)*V(i,j-1) )*DXT ; % U*V*DX on Top Face
    GT3=0.5*( V(i,j)*V(i,j)+P(i,j)+V(i,j-1)*V(i,j-1)+P(i,j-1) )*DXT ; % (V^2+p)*DX on Top Face
       
 
%--------------------------------------------------------------------------
%-Viscous Terms---Second Order Term Computed from Secondary Cells
    %-------Compute R terms(Du/Dx)
        %R Term of Left(Computed on Left(L) Secondary Cell)
        RL1=0;
        RL2=( 0.5*(Usw+U(i,j))*DYLDA + 0.5*(U(i,j)+Unw)*DYLCD + 0.5*(Unw+U(i-1,j))*DYLBC + 0.5*(U(i-1,j)+Usw)*DYLAB) *DYL/AL; %Du/Dx on Left Face
        RL3=( 0.5*(Vsw+V(i,j))*DYLDA + 0.5*(V(i,j)+Vnw)*DYLCD + 0.5*(Vnw+V(i-1,j))*DYLBC + 0.5*(V(i-1,j)+Vsw)*DYLAB)*DYL/AL;  %Dv/x on Left Face

        %R Term of Bottom (Computed on Bottom(B) Secondary Cell))
        RB1=0;
        RB2=( 0.5*(U(i,j+1)+Use)*DYBDA + 0.5*(Use+U(i,j))*DYBCD + 0.5*(U(i,j)+Usw)*DYBBC + 0.5*(Usw+U(i,j+1))*DYBAB) *DYB/AB; %Du/Dx on Bottom Face
        RB3=( 0.5*(V(i,j+1)+Vse)*DYBDA + 0.5*(Vse+V(i,j))*DYBCD + 0.5*(V(i,j)+Vsw)*DYBBC + 0.5*(Vsw+V(i,j+1))*DYBAB)   *DYB/AB; %Dv/x on Bottom Face


        %R Term of Right(Computed on Right(R) Secondary Cell))
        RR1=0;
        RR2=( 0.5*(Use+U(i+1,j) )*DYRDA + 0.5*(U(i+1,j)+Une) *DYRCD + 0.5*(Une+U(i,j) )*DYRBC + 0.5*(U(i,j)+Use )*DYRAB ) *DYR/AR; %Du/Dx on Right Face
        RR3=( 0.5*(Vse+V(i+1,j))*DYRDA + 0.5*(V(i+1,j)+Vne)*DYRCD + 0.5*(Vne+V(i,j))*DYRBC + 0.5*(V(i,j)+Vse)*DYRAB) *DYR/AR; %Du/Dx on Right Face

        %R Term of Top(Computed on Top(R) Secondary Cell))
        RT1=0;
        RT2=( 0.5*(U(i,j)+Une)*DYTDA + 0.5*(Une+U(i,j-1))*DYTCD + 0.5*(U(i,j-1)+Unw)*DYTBC + 0.5*(Unw+U(i,j))*DYTAB) *DYT/AT; %Du/Dx on Top Face
        RT3=( 0.5*(V(i,j)+Vne)*DYTDA + 0.5*(Vne+V(i,j-1))*DYTCD + 0.5*(V(i,j-1)+Vnw)*DYTBC + 0.5*(Vnw+V(i,j))*DYTAB) *DYT/AT; %Du/Dx on Top Face

        
    %-------Compute s terms(Du/Dy)
        %S Term of Left(Computed on Left(L) Secondary Cell)
        SL1=0;
        SL2=-( 0.5*(Usw+U(i,j))*DXLDA + 0.5*(U(i,j)+Unw)*DXLCD + 0.5*(Unw+U(i-1,j))*DXLBC + 0.5*(U(i-1,j)+Usw)*DXLAB) *DXL/AL; %Du/Dy on Left Face
        SL3=-( 0.5*(Vsw+V(i,j))*DXLDA + 0.5*(V(i,j)+Vnw)*DXLCD + 0.5*(Vnw+V(i-1,j))*DXLBC + 0.5*(V(i-1,j)+Vsw)*DXLAB)*DXL/AL;  %Dv/Dy on Left Face


        %S Term of Bottom(Computed on Bottom(B) Secondary Cell)
        SB1=0;
        SB2=-( 0.5*(U(i,j+1)+Use)*DXBDA + 0.5*(Use+U(i,j))*DXBCD + 0.5*(U(i,j)+Usw)*DXBBC + 0.5*(Usw+U(i,j+1))*DXBAB) *DXB/AB;%Du/Dy on Bottom Face
        SB3=-( 0.5*(V(i,j+1)+Vse)*DXBDA + 0.5*(Vse+V(i,j))*DXBCD + 0.5*(V(i,j)+Vsw)*DXBBC + 0.5*(Vsw+V(i,j+1))*DXBAB)   *DXB/AB;%Dv/Dy on Bottom Face
 
 
        %S Term of Right(Computed on Right(R) Secondary Cell)
        SR1=0;
        SR2=-( 0.5*(Use+U(i+1,j))*DXRDA + 0.5*(U(i+1,j)+Une)*DXRCD + 0.5*(Une+U(i,j))*DXRBC + 0.5*(U(i,j)+Use)*DXRAB) *DXR/AR; %Du/Dy on Right Face
        SR3=-( 0.5*(Vse+V(i+1,j))*DXRDA + 0.5*(V(i+1,j)+Vne)*DXRCD + 0.5*(Vne+V(i,j))*DXRBC + 0.5*(V(i,j)+Vse)*DXRAB) *DXR/AR; %Dv/Dy on Right Face
        %S Term of Top(Computed on Top(T) Secondary Cell)
        ST1=0;
        ST2=-( 0.5*(U(i,j)+Une)*DXTDA + 0.5*(Une+U(i,j-1))*DXTCD + 0.5*(U(i,j-1)+Unw)*DXTBC + 0.5*(Unw+U(i,j))*DXTAB) *DXT/AT; %Du/Dy on Top Face
        ST3=-( 0.5*(V(i,j)+Vne)*DXTDA + 0.5*(Vse+V(i,j-1))*DXTCD + 0.5*(V(i,j-1)+Vnw)*DXTBC + 0.5*(Vnw+V(i,j))*DXTAB) *DXT/AT; %Dv/Dy on Top Face

             

%*************************************************************************
%----------Compute F and G for each Equ.------------------------------
F1=FL1+FB1+FR1+FT1;
F2=FL2+FB2+FR2+FT2;
F3=FL3+FB3+FR3+FT3;
%
G1=GL1+GB1+GR1+GT1;
G2=GL2+GB2+GR2+GT2;
G3=GL3+GB3+GR3+GT3;

%----------Compute R and S for each Equ.------------------------------
R1=RL1+RB1+RR1+RT1;
R2=RL2+RB2+RR2+RT2;
R3=RL3+RB3+RR3+RT3;
%
S1=SL1+SB1+SR1+ST1;
S2=SL2+SB2+SR2+ST2;
S3=SL3+SB3+SR3+ST3;

%----------Return RHS for each Equ.-------------------------------------
QP=((R1-S1)/Ren-(F1-G1))/A(i,j);
QU=((R2-S2)/Ren-(F2-G2))/A(i,j);
QV=((R3-S3)/Ren-(F3-G3))/A(i,j);


end

