function [P,U,V] = BcsCavity( n,m,P,U,V )
%sets bcs on each side of the domain

%Velocities on each boundary
%*********************************************
     %Left Boundry(Wall)
 U(1,:)=-U(2,:);  % Zero U @ Wall
 V(1,:)=-V(2,:);  % Zero V @ Wall 
 P(1,:)=P(2,:);    %Zero Gradient @ Wall

      %Right Boundry(Wall)
 U(n+2,:)=-U(n+1,:);  % Zero U @ Wall
 V(n+2,:)=-V(n+1,:);  % Zero V @ Wall 
 P(n+2,:)=P(n+1,:);    %Zero Gradient @ Wall
Ut=1; %Top Wall Velocity
    %Top Boundry(Wall)
 U(:,m+2)=2*Ut-U(:,m+1);   % U=Ut @ Wall
 V(:,m+2)=-V(:,m+1);   % Zero V @ Wall 
 P(:,m+2)=P(:,m+1);    %Zero Gradient @ Wall

    %Bottom Boundry(Wall)
 U(:,1)=-U(:,2);  % Zero U @ Wall
 V(:,1)=-V(:,2);  % Zero V @ Wall 
 P(:,1)=P(:,2);    %Zero Gradient @ Wall





end

