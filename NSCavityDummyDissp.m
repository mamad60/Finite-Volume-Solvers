clear
clc
%Solves Steady Channnel Flow by Artificial Compressiblity Method
%Using Finite Volume Method
%and adds Artificial Dissipation for Numerical Stability and to overcome 
%Odd Even Coupling Issue
%Jamson 4th Order Dissipation
%BY Mohammad Aghakhani, 2012

%Definition of varibles-------------------------
%Location of points
global m n;
global  X Y;
global A;
global XC YC;
global U V P;
global Beta;
global Ren;
%----------------Inputs----------------------
L=1;            %Lenght of domain (x direction)
H=1;            %Height of domain(y dirction)
Beta=1.5;
Ren=100;        % Free Stream Reynolds Number
n=60 ;          %No. of cells in x direction
m=30;            %number of cell in y direction
MaxIT=50000;      %Maximum allowed iteration
%Dt=.000001;             %time step
eps=1e-3;      %Convergence criteria
u0=1;            %Initial value of u
v0=0;            %Initial value of v
p0=0.01;            %Initial value of p
CFL=0.25; %Courant Number
err=zeros(1,MaxIT);
err(1)=1000; %Error in two con. time step
epsP=0.005;  %Dissipation Coeficient P
epsU=0.005;  %Dissipation Coeficient U
epsV=0.005;  %Dissipation Coeficient in Y Direction 
%---------------------------------------------
[X,Y]=deal(zeros(m+3,n+3)');          %position of grid points

[XC,YC]=deal(zeros(m+2,n+2)');   %position of cell centers(solution variables are stored)


[A]=deal(zeros(m+2,n+2)');       %Area of Primary cells

%Solution variables
    %Solution varibles @ Current  iteration
[U,V,P]=deal(zeros(m+2,n+2)');
    %Solution varibles @ Previous iteration
[Uold,Vold,Pold]=deal(zeros(m+2,n+2)');   
%Right hand sites of Eqs
RHSP=0;
RHSU=0;
RHSV=0;
%----------------------------------------------------
%Grid gerneration
GridDummy(L,H);
%Initiate the Soloution
Initiate(u0,v0,p0);
[P,U,V ] = BcsCavity( n,m,P,U,V );
%Determine Maximum allowabe time step size
dL=L/n;
dH=H/m;
Dt=CFL_Test(Beta,CFL,dL,dH,U,V);
fprintf(1,'Maximum Allowable Time Step=%2.6e\n',Dt);
disp('Press any key')
pause
%Begin Iteration
IT=1;              %current iteration number
while((IT<MaxIT)&&(err(IT)>eps))
    IT=IT+1;
    %shift solution from old iteration
    Pold=P;
    Uold=U;
    Vold=V;
    %Apply Boundary Conditions
    [P,U,V]=BcsCavity(n,m,P,U,V);
    for j=2:m+1
        for i=2:n+1
           %Compute Fluxes From  values From P,V,T
            [ RHSP,RHSU,RHSV ] = Flux( i,j );
            %Calcualte Dissipation
            [ DISP,DISU,DISV ] = DISSP( i,j,epsP,epsU,epsV);
           %Add diissioation to RHS
            RHSP=RHSP-DISP;
            RHSU=RHSU-DISU;
            RHSV=RHSV-DISV;
           %Calculate Solution in new time step Using Explicit Euler
            [P(i,j)]=Euler(P(i,j),Dt,RHSP);
            [U(i,j)]=Euler(U(i,j),Dt,RHSU);
            [V(i,j)]=Euler(V(i,j),Dt,RHSV);

        end
    end
    
    %Calculate Error
    errP=max(max(abs((P-Pold))));
    errU=max(max(abs((U-Uold))));
    errV=max(max(abs((V-Vold))));
    erri=max(errU,errV);
    err(IT)=max(erri,errP);
    fprintf(1,'IT=%i   Error=%2.6e\n',IT,err(IT));
end

if(err(IT-1)>1000)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Time Step, Beta Or other Parmeters and Run the Code Again\n');
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
P=P(2:n+1,2:m+1);
U=U(2:n+1,2:m+1);
V=V(2:n+1,2:m+1);

XC=XC(2:n+1,2:m+1);
YC=YC(2:n+1,2:m+1);


disp('****************************************************************')
fprintf(1,'Beta=%2.2f       IT=%u       Dt=%2.2e        IT*Dt=%2.2f    Beta*Dt=%2.5f\n ',Beta, IT, Dt, IT*Dt, Beta*Dt);
Min=sum((U(1,:))*(dH));
Mout=sum((U(n,1))*(dH));
Mass_inbalance=Mout-Min
disp('Press any key')
pause


figure
[C1,h1] = contourf(XC,YC,U,20);
%text_handle = clabel(C1,h1);
colorbar
title('CONTOURES OF u(x,y)');
xlabel('x')
ylabel('y')
axis fill

figure
surf(XC,YC,P)
title('Surface of Pressure')
axis fill

figure
surf(XC,YC,U)
title('Surface of U(x) Velocity)')
axis fill


figure
surf(XC,YC,V)
title('Surface of V(y Velocity)')


% figure
% hold on
% plot(XC(1,1:m),U(1,1:m),'-* r')
% plot(Y(ceil(n/10),1:m),U(ceil(n/10),1:m),'- g')
% plot(Y(ceil(n/5),1:m),U(ceil(n/5),1:m),'-S m')
% plot(Y(ceil(n/3),1:m),U(ceil(n/3),1:m),'- y')
% plot(Y(ceil(n/2),1:m),U(ceil(n/2),1:m),'- b')
% plot(Y(ceil(3*n/2),1:m),U(ceil(3*n/2),1:m),'-. c')
% plot(Y(ceil(4*n/3),1:m),U(ceil(4*n/3),1:m),'- k')
% plot(Y(ceil(n),1:m),U(ceil(n),1:m),'-- r')
% legend('x=0','x=n/10','x=n/5','x=n/3','x=n/2','x=3*n/2','x=4*n/2','x=n',1)
% title('Profiles of the U velocity across Cavity Section')
% xlabel('H')
% ylabel('U Velocity')
% xlim([0 H]);
% ylim('auto')
% hold off;

figure
quiver(XC,YC,U*5,V*5)
xlim([0 L]);
ylim([0 H]);
time=IT*Dt;
axis fill;
% title('Velocity Vectors')
% streamline(XC,YC,U,V,L/2,H/2);
% title('Stream Lime Plot')

