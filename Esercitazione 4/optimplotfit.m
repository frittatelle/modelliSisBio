function stop = optimplotfit(x,optimValues,state,data)
% OPTIMPLOTFIT Plot fitting at each iteration.
% Data:
% x,z,objFun,title,pauseTime

stop = false;
t = data.x;
z = data.z;
[~,y] = data.objFun(x);
disp('ciao')


switch state
    case 'iter'
        tLength = length(t);     
        if optimValues.iteration == 0  
            figure
            plotData = plot(t,y);
            hold on
            plot(t,z,'*')
            hold off
            title(data.title);
            set(gca,'xlim',[0,1 + tLength])
            set(plotData,'Tag','optimplotx');
        else
            plotData = findobj(get(gca,'Children'),'Tag','optimplotx');
            set(plotData,'Xdata',t,'Ydata',y);          
        end
        
pause(data.pauseTime)        
        
end







