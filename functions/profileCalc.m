function [ par] = profileCalc( par , dist, alignm, namestrc,cutths)


if alignm == 1;
    figure
    % 1 diameter above
    x = unique(par.x);
    y = par.avgMagnitude(par.centroPiv0(2)-dist*round(par.esd/(par.wdw*par.mmpropix)),:);
    
    findcutPos = min(find(x >= x(par.centroPiv0(1))+cutths));
    findcutNeg = max(find(x <= x(par.centroPiv0(1))-cutths));
    
    y(isnan(y)) = 0;
    yleft = y(findcutNeg:par.centroPiv0(1));
    %yright = fliplr(y(par.centroPiv0(1):end));
    yright = y(par.centroPiv0(1):findcutPos);
    xleft = x(findcutNeg:par.centroPiv0(1));
    %xright = fliplr(x(par.centroPiv0(1):end));
    xright = x(par.centroPiv0(1):findcutPos);
    
    if isempty(yleft) || isempty(yright)
        disp('Matlab Error: Out of bounds, set to max bounds')
        y(isnan(y)) = 0;
        yleft = y(1:par.centroPiv0(1));
        %yright = fliplr(y(par.centroPiv0(1):end));
        yright = y(par.centroPiv0(1):end);
        xleft = x(1:par.centroPiv0(1));
        %xright = fliplr(x(par.centroPiv0(1):end));
        xright = x(par.centroPiv0(1):end);
    end
    
    
    Umax = quantile(yleft,0.90);
    UmaxRight = quantile(yright,0.90);
    
    thsh = Umax * 0.95;
    thshRight = UmaxRight * 0.95;
    
    xvaltemp = xleft(yleft>thsh);
    ind1 = find(xleft == max(xvaltemp));
    ind2 = ind1+1;
    
    xvaltempRight = xright(yright>thshRight);
    ind1Right = find(xright == min(xvaltempRight));
    ind2Right = ind1Right-1;
    
    boundaryThick = interp1([yleft(ind1) yleft(ind2)],[xleft(ind1) xleft(ind2)],thsh);
    boundaryThickRight = interp1([yright(ind1Right) yright(ind2Right)],[xright(ind1Right) xright(ind2Right)],thshRight);
    
    par.([namestrc 'Umax']) = Umax;
    par.([namestrc 'UmaxRight']) = UmaxRight;
    
    %par.([namestrc 'Dist']) = dist; %Soeren
    %par.([namestrc 'Alignm']) = alignm; %Soeren
    par.([namestrc 'BoundaryThickness']) = max(xleft)-boundaryThick;
    par.([namestrc 'BoundaryThicknessRight']) = boundaryThickRight-min(xright);
    
    fprintf('\n\n Boundary thickness is %.3g and on the right side %.3g mm\n',[max(xleft)-boundaryThick boundaryThickRight-min(xright)])
    fprintf('\n\n Umax terminal velocity is %.6g and on the right side %.6g m/d\n',[Umax UmaxRight])
      
    
    hold on
    plot(x,y,'o-')
    plot(xleft,yleft,'b.','MarkerSize',20)
    plot(xright,yright,'b.','MarkerSize',20)
    plot(boundaryThick,thsh,'r.','MarkerSize',25)
    plot(boundaryThickRight,thshRight,'r.','MarkerSize',25)
    
    title(namestrc)
else
    figure
    % 1 diameter above
    y = flipud(-unique(par.y));
    x = flipud(par.avgMagnitude(:,par.centroPiv0(1)-dist*round(par.esd/par.wdw*par.mmpropix)));
    
    len = length(y);
    
    x(isnan(x)) = 0;
 
    ybot = y(1:len-par.centroPiv0(2));
    xbot = x(1:len-par.centroPiv0(2));
    
    findcutPos = max(find(ybot <= y(len-par.centroPiv0(2))-cutths));    

    ybottom = ybot(findcutPos:end);
    xbottom = xbot(findcutPos:end);
     
    if isempty(ybottom)
        disp('Matlab Error: Out of bounds, set to max bounds')
        ybottom = ybot;
        xbottom = xbot;
    end
    
    Umax = quantile(xbottom,0.85);
    thsh = Umax * 0.95;
    
    yvaltemp = ybottom(xbottom>thsh);
    ind1 = find(ybottom == max(yvaltemp));
    ind2 = ind1+1;
    
    boundaryThick = interp1([xbottom(ind1) xbottom(ind2)],[ybottom(ind1) ybottom(ind2)],thsh);
    
    par.([namestrc 'Umax']) = Umax;
    par.([namestrc 'UmaxRight']) = 'NaN';
    
    %par.([namestrc 'Dist']) = dist;
    %par.([namestrc 'Alignm']) = alignm;
    par.([namestrc 'BoundaryThickness']) = max(ybottom)-boundaryThick;
    par.([namestrc 'BoundaryThicknessRight']) = 'NaN';
    
    %else
    %    par.boundaryThick = (max(xleft)-boundaryThick);
    %    par.dist = [alignm num2str(dist)];
    %    par.align = alignm;
    %    par.Umax = Umax;
    %end
    hold on
    plot(xbottom,ybottom,'b.','MarkerSize',25)
    plot(x,y,'o-')
    plot(thsh,boundaryThick,'r.','MarkerSize',25)
    title(namestrc)
end

if isfield(par,'Slice_names')
    par.Slice_names{end+1} = {namestrc};
else
    par.('Slice_names') = {namestrc};
end
