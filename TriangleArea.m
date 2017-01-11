function [A] = TriangleArea(xa,ya,xb,yb,xc,yc)
A=abs(xa*(yc-yb)+xb*(ya-yc)+xc*(yb-ya))/2;

