function [ DISP,DISU,DISV ] = DISSP( I,J,EP,EU,EV)

global m n;
global A;
global U V P;

if ((I<4) || I>(n-3)|| (J<4) || J>(m-3))
    DISU=0.0;
    DISV=0.0;
    DISP=0.0;
else
    DISUI=EU*(U(I+2,J)-4.0*U(I+1,J)+6.0*U(I,J) -4.0*U(I-1,J)+U(I-2,J));
    
    DISUJ=EU*(U(I,J+2) -4.0*U(I,J+1)+6.0*U(I,J) -4.0*U(I,J-1)+U(I,J-2));
    
    
    DISVI=EV*(V(I+2,J)-4.0*V(I+1,J)+6.0*V(I,J) -4.0*V(I-1,J)+V(I-2,J));
    
    DISVJ=EV*(V(I,J+2)-4.0*V(I,J+1)+6.0*V(I,J) -4.0*V(I,J-1)+V(I,J-2));
    
    
    DISPI=EP*(P(I+2,J)-4.0*P(I+1,J)+6.0*P(I,J) -4.0*P(I-1,J)+P(I-2,J));
    
    DISPJ=EP*(P(I,J+2)-4.0*P(I,J+1)+6.0*P(I,J)-4.0*P(I,J-1)+P(I,J-2));
    
    DISU=(DISUI+DISUJ)/A(I,J);
    DISV=(DISVI+DISVJ)/A(I,J);
    DISP=(DISPI+DISPJ)/A(I,J);
    
    
end
% %----------------------------------------------
% 
% 
% if (J==2) || (J==3)
%     if (I<n-3)
%         DISUI=EU*(U(I+2,J)-4.0*U(I+1,J)+6.0*U(I,J) -4.0*U(I-1,J)+U(I-2,J));
%         DISVI=EV*(V(I+2,J)-4.0*V(I+1,J)+6.0*V(I,J) -4.0*V(I-1,J)+V(I-2,J));
%         DISPI=EP*(P(I+2,J)-4.0*P(I+1,J)+6.0*P(I,J) -4.0*P(I-1,J)+P(I-2,J));
%         
%         DISUJ=EU*(U(I,J+4) -4.0*U(I,J+3)+6.0*U(I,J+2) -4.0*U(I,J+1)+U(I,J));
%         DISVJ=EV*(V(I,J+4) -4.0*V(I,J+3)+6.0*V(I,J+2) -4.0*V(I,J+1)+V(I,J));
%         DISPJ=EP*(P(I,J+4)-4.0*P(I,J+3)+6.0*P(I,J+2)-4.0*P(I,J+1)+P(I,J));
%                 
%     end
%     if (I==n-2) || (I==n-1)
%         DISUI=EU*(U(I,J)-4.0*U(I-1,J)+6.0*U(I-2,J) -4.0*U(I-3,J)+U(I-4,J));
%         DISVI=EV*(V(I,J)-4.0*V(I-1,J)+6.0*V(I-2,J) -4.0*V(I-3,J)+V(I-4,J));
%         DISPI=EP*(P(I,J)-4.0*P(I-1,J)+6.0*P(I-2,J) -4.0*P(I-3,J)+P(I-4,J));
%         
%         DISUJ=EU*(U(I,J+4) -4.0*U(I,J+3)+6.0*U(I,J+2) -4.0*U(I,J+1)+U(I,J));
%         DISVJ=EV*(V(I,J+4) -4.0*V(I,J+3)+6.0*V(I,J+2) -4.0*V(I,J+1)+V(I,J));
%         DISPJ=EP*(P(I,J+4)-4.0*P(I,J+3)+6.0*P(I,J+2)-4.0*P(I,J+1)+P(I,J));
%               
%     end
%     DISU=(DISUI+DISUJ)/A(I,J);
%     DISV=(DISVI+DISVJ)/A(I,J);
%     DISP=(DISPI+DISPJ)/A(I,J);
% 		 
% end
%-------------------------------------------------


