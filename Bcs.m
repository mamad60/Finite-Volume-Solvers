function [  ] = Bcs(Tr,Tl,Tt,Tb)
%sets bcs on each side of the domain

global n m;
global  T;            

     T(1,:)=2*Tl-T(2,:);
     T(n+2,:)=2*Tr-T(n+1,:);
%Neumann Boundary Condition  
%       T(1,:)=T(2,:);
%       T(n+2,:)=T(n+1,:);
%Dirichlet Boundary Condition     
     T(:,1)=2*Tt-T(:,2);
     T(:,m+2)=2*Tb-T(:,m+1);
%Neumann Boundary Condition  
%      T(:,m+2)=T(:,m+1);
%      T(:,1)=T(:,2);


end

