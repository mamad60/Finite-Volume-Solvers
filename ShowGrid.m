function  ShowGrid(ShowDummy)
%Return right hand side for each eqs.
global n m;
global  X Y;
global XC YC;

%Horizental Lines
for j=2:m+2
    for i=1:n+2
        if (i==1 || i==n+2)
            if ShowDummy
                line([X(i,j) X(i+1,j)],[Y(i,j) Y(i+1,j)],'Color',[0.8 0.8 0.8]);
            end
        else
            line([X(i,j) X(i+1,j)],[Y(i,j) Y(i+1,j)],'Color','b');
        end
    end
    
end

%Vertical Lines
for i=2:n+2
    for j=1:m+2
        if (j==1 || j==m+2)
            if ShowDummy
                line([X(i,j) X(i,j+1)],[Y(i,j) Y(i,j+1)],'Color',[0.8 0.8 0.8]);
            end
        else
            line([X(i,j) X(i,j+1)],[Y(i,j) Y(i,j+1)],'Color','b');
        end
    end
    
end
%Show Cell Centers
if ShowDummy
    hold on
    plot(XC,YC,'.k')
end


if ShowDummy
    for i=1:n+2
        j=1;
        line([X(i,j) X(i+1,j)],[Y(i,j) Y(i+1,j)],'Color',[0.8 0.8 0.8]);
        j=m+3;
        line([X(i,j) X(i+1,j)],[Y(i,j) Y(i+1,j)],'Color',[0.8 0.8 0.8]);
    end
    for j=1:m+2
        i=1;
        line([X(i,j) X(i,j+1)],[Y(i,j) Y(i,j+1)],'Color',[0.8 0.8 0.8]);
        i=n+3;
        line([X(i,j) X(i,j+1)],[Y(i,j) Y(i,j+1)],'Color',[0.8 0.8 0.8]);
    end
end

