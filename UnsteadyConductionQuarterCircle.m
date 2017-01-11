%Solves UnSteady Conduction in a Nonorthogonal quadilateral Grid
%Using Finite Volume Method, Easlily Customizable for Triangular Grids
%The Domain Could be any domain with Structed Grid in a Quarter Circle
%This version Uses Dummy Cells for Imposing Boundary Condition
%By Mohammad Aghakhani, 2012

clear
clc

%----------------Inputs----------------------
Ro=2; %Outer circle radius
Ri=1;  %Inner circle  radius
Theta0=pi/2; %Thata varies between 0 and Thata 0
global n m;
global Alpha;
Alpha=1e-3;
n=50 ;          %No. of cells Around circles (Thetha dirction)
m=20;            %No. of cell Normal to Circles (R direction)
MaxIT=100000;      %Maximum allowed iteration
eps=1e-6;      %Convergence criteria
t0=0;            %Initial value of T(Value of T at time t=0)
Dt=.6 ;             %time step
global F0;
% F0=.1;              %fourier no.
err=zeros(1,MaxIT);
err(1)=1000;    %Error in two con. time step
DisplayGrid=1;   %Show Grid in Contour Plot
%---------------------------------------------

%Definition of varibles-------------------------
%Location of points
global  X Y;
global XC YC;
global T;
global A;

[X,Y]=deal(zeros(m+3,n+3)');    %position of grid points(Two Extra points for Dummy Cells)

[XC,YC]=deal(zeros(m+2,n+2)');  %position of cell centers(solution variables are stored)


[A]=deal(zeros(m+2,n+2)');       %Area of Primary cells

%Solution variables
[T]=deal(zeros(m+2,n+2)');      %Solution Varible
[Told]=deal(zeros(m+2,n+2)');   %Solution Varible@old Iteration
%----------------------------------------------------
% Test for Fourier Number
dR=(Ro-Ri)/m; %Grid Spacing in Radial Direction
th=Theta0*Ri/n; %Grid Spacing in Circumcal Direction
F0=Alpha*Dt/(dR*th); %Output Fourier Number
fprintf(1,'The Time Step is:%2.6e Fourier Number is:%2.6e\n',Dt,F0);
disp('Press any key')
pause
%Grid gerneration
GridQuarterCircle(Ro,Ri,Theta0)
%Display Grid
figure;
ShowGrid(1,1,1,0);
xlabel('x');
ylabel('y');
title('Grid with Dummy cells(Grey)');
%Initiate the Soloution
InitiateT(t0);
%Temperatures on each boundary (Boundary Conditions)
%*********************************************
%left Boundry
Tl=100;
%Bottom Boundry
Tb=100;
%Right Boundry
Tr=0;

%Top Boundry
Tt=100;

%*********************************************
% Begin Euler Iteration
%  TempCountor=figure;
%  TCA = axes('parent',TempCountor);
%  errT=100;
%  TempError=figure;
%  TEA = axes('parent',TempError);
IT=1;              %current iteration number

while((IT<MaxIT)&&(err(IT)>eps))
    IT=IT+1;
    
    %shift solution from old iteration
    Told=T;
    %--------------------------------------
    
    %Appply Bcs
    BcsNeumann(Tr,Tl,Tt,Tb);
    for j=2:m+1
        for i=2:n+1
            [QTA] = FluxT1( i,j );  %Calculate Flux at
            T(i,j)=Told(i,j)+( QTA )*Dt;
        end
    end
    %Calculate Error
    errT=max(max(abs((T-Told))));
    err(IT)=errT;
    fprintf(1,'IT=%i   Error=%2.6e\n',IT,err(IT));
    
end

if(err(IT-1)>1000)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Time Step Or other Parmeters and Run the Code Again\n');
    disp('Press any key')
    pause
    return;
end
if(err(IT)<eps)
    fprintf(1,'Converged in %i Iterations\n',IT);
else
    disp('Maximum Iteration Number Reached');
    plot(1:IT,log10(error(1:IT)),'- r');
    xlabel('Iteration');
    ylabel('Log(error)');
    title('Convergence History');
    pause;
    return;
end

%Plot Convergence History
figure
plot( 1:IT,log10(err(1:IT)),'-. g');
xlabel('Iteration')
ylabel('Log10(Error)')
title('Error History')

%Exclude Ghost(Dummy) Cells from Calcualtions
T=T(2:n+1,2:m+1)';
XC=XC(2:n+1,2:m+1)';
YC=YC(2:n+1,2:m+1)';

%Plot Results
figure
if DisplayGrid
    hold on;
    ShowGrid(1,1,1,0);
end

[C1,h1] = contourf(XC,YC,T,20);
%text_handle = clabel(C1,h1,'manual');
colorbar
title('CONTOURES OF Temperature');
xlabel('x')
ylabel('y')
axis fill
if DisplayGrid
    hold on
    ShowGrid(1,1,1,0);
    hold off
end

figure
surf(XC,YC,T)
title('Surface of Temperature')
axis fill

fprintf(1,'Final Time is %2.6f\n',Dt*IT);