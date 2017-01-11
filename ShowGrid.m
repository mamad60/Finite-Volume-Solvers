function  ShowGrid(ShowDummy,ShowCenter,ShowLines,ShowNumbers)
%Return right hand side for each eqs.
global n m;
global  X Y;
global XC YC;

if ShowLines
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
end
%Show Cell Centers
if ShowCenter
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
if ShowNumbers
    for j=1:m+2
        for i=1:n+2
            text(XC(i,j),YC(i,j),strcat(num2str(i),',',num2str(j)),'FontSize',8)
         end
    end
    text(2*XC(1,1),0.5*(YC(n,m)+YC(1,1)),'I \uparrow','FontSize',18)
    text(0,2*YC(1,1),'J \rightarrow','FontSize',18);


end



