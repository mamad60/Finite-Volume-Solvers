function [] = GridQuarterCircle(Ro,Ri,Theta0)
%Generates a uniform grid and returns the positions of Grid points and cell
%centers 
global m n;
global  X Y;
global  A;
global XC YC;

dR=(Ro-Ri)/m;
dTheta=Theta0/n;


for j=1:m+3
        for i=1:n+3
            R=Ri+(m-j+2)*dR;
            Thetha=Theta0-(i-2)*dTheta;
            X(i,j)=R*cos(Thetha);
            Y(i,j)=R*sin(Thetha);
        end
end


%X and Y of Center points of of Primary cells

for j=1:m+2
        for i=1:n+2
            YC(i,j)=0.25*( Y(i,j)+Y(i+1,j)+Y(i,j+1)+Y(i+1,j+1) );
            XC(i,j)=0.25*( X(i,j)+X(i+1,j)+X(i,j+1)+X(i+1,j+1) );
        end
end

%Areas of primary cells
for j=1:m+2
        for i=1:n+2
            A(i,j)=0.5*(  (X(i+1,j+1)-X(i,j) )*( Y(i+1,j)-Y(i,j+1)) - (X(i+1,j)-X(i,j+1))*(Y(i+1,j+1)-Y(i,j)) );
        end
end

end

