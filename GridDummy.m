function [] = GridDummy(L,H)
%Generates a uniform grid and returns the positions of Grid points and cell
%centers 
global m n;
global  X Y;
global  A;
global XC YC;


dL=L/n;
dH=H/m;


for j=1:m+3
        for i=1:n+3
            X(i,j)=(i-2)*dL;
            Y(i,j)=(j-2)*dH;
        end
end


%X and Y of Center points of of Primary cells

for j=1:m+2
        for i=1:n+2
            YC(i,j)=0.25*( Y(i,j)+Y(i+1,j)+Y(i,j+1)+Y(i+1,j+1) );
            XC(i,j)=0.25*( X(i,j)+X(i+1,j)+X(i,j+1)+X(i+1,j+1) );
        end
end

A(:)=dL*dH;
end

