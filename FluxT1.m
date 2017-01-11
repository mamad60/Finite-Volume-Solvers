function [ QT] = FluxT1( i,j )
%Return right hand side for each eqs.
% global n m;
global  X Y;
global XC YC;
global T;
global Alpha;
global A;
% global F0;


%Note:i,j in soloution shows cell centers while in X,Y refer to grid
%points

%*********************************************
%Begin Geometrical Values Calculations
%--------------------------------------------------------------
%Left, Bottom ,Right & Top DeltaY of Primary cells
DYT=Y(i,j)-Y(i+1,j);  %Deta Y of Top Face
DYL=Y(i,j+1)-Y(i,j);  %Deta Y of Left Face
DYB=Y(i+1,j+1)-Y(i,j+1); %Deta Y of Bottom Face
DYR=Y(i+1,j)-Y(i+1,j+1); %Deta Y of Right Face

%Left, Bottom ,Right & Top DeltaX of Primary cells
DXT=X(i,j)-X(i+1,j);  %Deta X of Top Face
DXL=X(i,j+1)-X(i,j);  %Deta X of Left Face
DXB=X(i+1,j+1)-X(i,j+1);  %Deta X of Bottom Face
DXR=X(i+1,j)-X(i+1,j+1);  %Deta X of Right Face 

%-------------------------------------------------------------


%Y of Center points of Left, Bottom ,Right & Top of Primary cells


    YL=0.25*( Y(i-1,j)+Y(i,j)+Y(i-1,j+1)+Y(i,j+1) );  %YC(i-1,j)

    YB=0.25*( Y(i,j+1)+Y(i+1,j+1)+Y(i,j+2)+Y(i+1,j+2) );  %YC(i,j+1)

    YR=0.25*( Y(i+1,j)+Y(i+2,j)+Y(i+1,j+1)+Y(i+2,j+1) );  %YC(i+1,j)

    YT=0.25*( Y(i,j-1)+Y(i+1,j-1)+Y(i,j)+Y(i+1,j) );  %YC(i,j+1


%X of Center points of Left, Bottom ,Right & Top of Primary cells


    XL=0.25*( X(i-1,j)+X(i,j)+X(i-1,j+1)+X(i,j+1) );  %XC(i-1,j)

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
            
            %R  Secondary Cell      
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
 
 AL=abs( 0.5*( (XC(i,j)-XL)*DYL-(YC(i,j)-YL)*DXL ) );  %Area of L Secondary cell
 AB=abs( 0.5*( (YB-YC(i,j))*DXB-(XB-XC(i,j))*DYB ) );  %Area of B Secondary cell
 AR=abs( 0.5*( (XC(i,j)-XR)*DYR-(YC(i,j)-YR)*DXR ) );  %Area of R Secondary cell
 AT=abs( 0.5*( (YT-YC(i,j))*DXT-(XT-XC(i,j))*DYT ) );  %Area of T Secondary cell
 
 %End of Geometrical Values Calculation
 %*********************************************       
%--------------------------------------------------------------           
  %Temparture of the Primary cell nodes(Corners) (internal cells Only)--------------------------
  %Required for Second Order Terms Calculations
  %Interpolated From Cell Center Values
         Tnw=0.25*( T(i-1,j-1)+T(i,j-1)+T(i-1,j)+T(i,j) );
         Tne=0.25*(  T(i,j-1)+T(i+1,j-1)+T(i,j)+T(i+1,j) );
         Tsw=0.25*( T(i-1,j)+T(i,j)+T(i-1,j+1)+T(i,j+1) );
         Tse=0.25*( T(i,j)+T(i+1,j)+T(i,j+1)+T(i+1,j+1) );
%--------------------------------------------------------------------------

%***********************************************
%Calaculations for internal(not next to boundaries)

%----------------------------------------------------------------
%-Second Order Terms

        %R of Left
        RL=( 0.5*(Tsw+T(i,j))*DYLDA + 0.5*(T(i,j)+Tnw)*DYLCD + 0.5*(Tnw+T(i-1,j))*DYLBC + 0.5*(T(i-1,j)+Tsw)*DYLAB) *DYL/AL;

        %R of Bottom
        RB=( 0.5*(T(i,j+1)+Tse)*DYBDA + 0.5*(Tse+T(i,j))*DYBCD + 0.5*(T(i,j)+Tsw)*DYBBC + 0.5*(Tsw+T(i,j+1))*DYBAB) *DYB/AB;

        


        %R of Right
        RR=( 0.5*(Tse+T(i+1,j) )*DYRDA + 0.5*(T(i+1,j)+Tne) *DYRCD + 0.5*(Tne+T(i,j) )*DYRBC + 0.5*(T(i,j)+Tse )*DYRAB ) *DYR/AR;

        %R of Top
        RT=( 0.5*(T(i,j)+Tne)*DYTDA + 0.5*(Tne+T(i,j-1))*DYTCD + 0.5*(T(i,j-1)+Tnw)*DYTBC + 0.5*(Tnw+T(i,j))*DYTAB) *DYT/AT;

        
%

        %S of Left
        SL=-( 0.5*(Tsw+T(i,j))*DXLDA + 0.5*(T(i,j)+Tnw)*DXLCD + 0.5*(Tnw+T(i-1,j))*DXLBC + 0.5*(T(i-1,j)+Tsw)*DXLAB) *DXL/AL;


        %S of Bottom
        SB=-( 0.5*(T(i,j+1)+Tse)*DXBDA + 0.5*(Tse+T(i,j))*DXBCD + 0.5*(T(i,j)+Tsw)*DXBBC + 0.5*(Tsw+T(i,j+1))*DXBAB) *DXB/AB;
 
 
        %S of Right
        SR=-( 0.5*(Tse+T(i+1,j))*DXRDA + 0.5*(T(i+1,j)+Tne)*DXRCD + 0.5*(Tne+T(i,j))*DXRBC + 0.5*(T(i,j)+Tse)*DXRAB) *DXR/AR;
        %S of Top
        ST=-( 0.5*(T(i,j)+Tne)*DXTDA + 0.5*(Tne+T(i,j-1))*DXTCD + 0.5*(T(i,j-1)+Tnw)*DXTBC + 0.5*(Tnw+T(i,j))*DXTAB) *DXT/AT;

             


%----------Compute R and S for each Equ.------------------------------
R=RL+RB+RR+RT;
%
S=SL+SB+SR+ST;

%----------Return RHS for each Equ.-------------------------------------
QT=(R-S)/A(i,j)*Alpha;

end

