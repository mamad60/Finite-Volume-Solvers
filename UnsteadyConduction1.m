clear
clc

%----------------Inputs----------------------
L=1;            %Lenght of domain (x direction)
H=1;            %Height of domain(y dirction0
global n m;
global Alpha;
Alpha=1e-3;
n=15 ;          %No. of cells in x direction
m=20;            %number of cell in y direction
MaxIT=10000;      %Maximum allowed iteration
ITN=1;              %current iteration number
eps=1e-6;      %Convergence criteria
t0=0;            %Initial value of T
Dt=.1;             %time step
global F0;
% F0=.1;              %fourier no.
%---------------------------------------------

%Definition of varibles-------------------------
            %Location of points
global  X Y;
global T;
global A;
global XC YC;

[X,Y]=deal(zeros(m+1,n+1)');                                             %position of grid points

[XC,YC]=deal(zeros(m,n)');                                             %position of cell centers(solution variables are stored 

            
[A]=deal(zeros(m,n)');                                             %Area of Primary cells

            %Solution variables
[T]=deal(zeros(m,n)');      
[Told]=deal(zeros(m,n)');   %old Iteration               
%----------------------------------------------------
%Grid gerneration
 Grid(L,H);
 %Initiate the Soloution
 InitiateT(t0);
 %Begin Iteration
  Dh=H/n;
 Dl=L/m; 
 F0=Alpha*Dt/(Dh*Dl)
 pause
%Velocities on each boundary
%*********************************************
    %left Boundry
Tl=1;
    %Bottom Boundry
Tb=0;
    %Right Boundry
Tr=0;

    %Top Boundry
Tt=1;

%*********************************************

 hold on
 errT=100;
 while(ITN<=MaxIT && errT>eps)
     %--------------------------------------
     %Appply Bcs 
     T(1,:)=(2/3)*T(2,:)+(1/3)*Tl;
     T(n,:)=(2/3)*T(n-1,:)+(1/3)*Tr;
     T(:,1)=(2/3)*T(:,2)+(1/3)*Tt;
%      T(:,1)=T(:,2);
     T(:,m)=(2/3)*T(:,m-1)+(1/3)*Tb;
%      T(:,m)=T(:,m-1);
    
     for j=2:m-1
         for i=2:n-1
             [QTA] = FluxT1( i,j );
             T(i,j)=Told(i,j)+( QTA )*Dt;
         end
     end
        %Calculate Error
    errT=max(max(abs((T-Told))));
    plot(ITN,log(errT));
    fprintf(1,'IT=%i  Error=%e\n',ITN,errT);
    ITN=ITN+1;

    
    Told=T;

 end
 hold off
 figure
 [C1,h1] = contourf(T',20);
 text_handle = clabel(C1,h1,'manual');
 colorbar
 figure
 plot(XC(:,n/2),T(:,n/2))
Dt*ITN
