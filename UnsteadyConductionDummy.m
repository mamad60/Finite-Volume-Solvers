%Solves UnSteady Conduction in a Nonorthogonal quadilateral Grid
%Using Finite Volume Method, Easily Customizable for Triangular Grids
%The Domain Could be any domain with Structed Grid
%This version Uses Dummy Cells for Imposing Boundary Condition
%By Mohammad Aghakhani, 2012
clear
clc
%----------------Inputs----------------------
L=1;            %Lenght of domain (x direction)
H=1;            %Height of domain(y dirction)
global n m;
global Alpha;
Alpha=1e-3;
n=20 ;          %No. of cells in x direction
m=15;           %number of cell in y direction
MaxIT=100000;   %Maximum allowed iteration
eps=1e-6;       %Convergence criteria
t0=0;           %Initial value of T
Dt=0.1;          %time step
global F0;
% F0=.1;              %fourier no.
err=zeros(1,MaxIT);
err(1)=1000;    %Error in two con. time step
DisplayGrid=1;   %Show Grid in Contour Plot
%---------------------------------------------

%Definition of varibles-------------------------
%Location of points
global  X Y;
global T;
global A;
global XC YC;

[X,Y]=deal(zeros(m+3,n+3)');    %position of grid points(Two Extra points for Dummy Cells)

[XC,YC]=deal(zeros(m+2,n+2)');  %position of cell centers(solution variables are stored)


[A]=deal(zeros(m+2,n+2)');       %Area of Primary cells

%Solution variables
[T]=deal(zeros(m+2,n+2)');      %Solution Varible
[Told]=deal(zeros(m+2,n+2)');   %Solution Varible@old Iteration
%----------------------------------------------------
%Fourier Number & Time Step
Dh=H/n;
Dl=L/m;
F0=Alpha*Dt/(Dh*Dl);
fprintf(1,'The Time Step is:%2.6e Fourier Number is:%2.6e\n',Dt,F0);
disp('Press any key')
pause
%Grid gerneration
GridDummy(L,H);
%Display Grid
figure;
ShowGrid(1,1,1,0);
xlabel('x');
ylabel('y');
title('Grid with Dummy cells(Grey)');
%Initiate the Soloution
InitiateT(t0);
%Temperatures on each boundary
%*********************************************
%left Boundry
Tl=0;
%Bottom Boundry
Tb=0;
%Right Boundry
Tr=1;

%Top Boundry
Tt=0;
% For any Change in Bcs and Grid change Coresponding Files
%*********************************************
%Begin Iteration
IT=1;          %current iteration number
while((IT<MaxIT)&&(err(IT)>eps))
    IT=IT+1;
    %shift solution from old iteration
    Told=T;
    %--------------------------------------
    %Appply Bcs
    Bcs(Tr,Tl,Tt,Tb);
    %      T(1,:)=2*Tl-T(2,:);
    %      T(n+2,:)=2*Tr-T(n+1,:);
    %
    %      T(:,1)=2*Tt-T(:,2);
    %      T(:,1)=T(:,2);
    %      T(:,m+2)=2*Tb-T(:,m+1);
    %      T(:,m+1)=T(:,m);
    for j=2:m+1
        for i=2:n+1
            [QTA] = FluxT1( i,j );
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
T=T(2:n+1,2:m+1);
XC=XC(2:n+1,2:m+1);
YC=YC(2:n+1,2:m+1);

%Plot Results
figure
[C1,h1] = contourf(XC,YC,T,20);
%text_handle = clabel(C1,h1,'manual');
colorbar
title('CONTOURES OF Temperature');
xlabel('x')
ylabel('y')
axis fill
if DisplayGrid
    hold on
    ShowGrid(0,1,1,0);
    hold off
end

figure
surf(XC,YC,T)
title('Surface of Temperature')
axis fill

figure
hold on
plot(T(1,:),YC(1,:),'-* r')
plot(T(floor(n/5),:),YC(floor(n/5),:),'- g')
plot(T(floor(n/4),:),YC(floor(n/4),:),'-S m')
plot(T(floor(n/3),:),YC(floor(n/3),:),'- y')
plot(T(floor(n/2),:),YC(floor(n/2),:),'- b')
plot(T(floor(n*2/3),:),YC(floor(n*2/3),:),'-. c')
plot(T(floor(n*3/4),:),YC(floor(n*3/4),:),'-.. g')
plot(T(n,:),YC(n,:),'-.. g')
legend('x=0','x=n/5','x=n/4','x=n/3','x=n/2','x=2*n/3','x=3*n/4','x=n',1)
title('Profiles of the Temperature across the Domain')
xlabel('Temperature');
ylabel('y Coordinate');
axis fill
hold off;
fprintf(1,'Final Time is %2.6f\n',Dt*IT);
