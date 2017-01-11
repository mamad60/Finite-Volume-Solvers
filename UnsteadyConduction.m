clear
clc

%----------------Inputs----------------------
L=1;            %Lenght of domain (x direction)
H=1;            %Height of domain(y dirction0
global n m;
global Alpha;
Alpha=1e-3;
n=10 ;          %No. of cells in x direction
m=10;            %number of cell in y direction
MaxIT=10000;      %Maximum allowed iteration
ITN=1;              %current iteration number
eps=1e-10;      %Convergence criteria
t0=0;            %Initial value of T
Dt=1;             %time step
global F0;
% F0=.1;              %fourier no.
%---------------------------------------------

%Definition of varibles-------------------------
%Location of points
global  X Y;
global T;
global A;
global XC YC;

[X,Y]=deal(zeros(m+1,n+1)');          %position of grid points

[XC,YC]=deal(zeros(m,n)'); %position of cell centers(solution variables are stored)


[A]=deal(zeros(m,n)'); %Area of Primary cells

%Solution variables
[T]=deal(zeros(m,n)');
[Told]=deal(zeros(m,n)');   %The Varible For holding at old Iteration
%----------------------------------------------------
%Grid gerneration
[ X,Y,dL,dH ] = Grid(m,n,L,H );%Initiate the Soloution
InitiateT(t0);
%Define Fourier Number

F0=Alpha*Dt/(dL*dH)
pause
hold on
errT=1000;
%Begin Iteration
while(ITN<=MaxIT && errT>eps)
    %--------------------------------------
    %shift solution from old iteration
    
    
    for j=1:m
        for i=1:n
            if( (i==1 && j==1) || (i==n && j==1) || (i==1 && j==m) || (i==n && j==m))  %4 Corners
                break;
            end
            [QTA] = FluxT( i,j );     %Calculate Flux
            T(i,j)=Told(i,j)+( QTA )*Dt;
        end
    end
    %Calculate Error and Plot Error
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
