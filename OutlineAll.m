%% cd C:\Users\admin\Desktop\piv_processing\pivclara\5\processed2
%
% Set Parameter
mscfolders = dir('St*');

for n = 1:length(mscfolders)

    cd(mscfolders(n).name)
    aggfolders = dir('./*');
    
    for agg_m = 3:length(aggfolders)
        cd(aggfolders(agg_m).name)
        
        repl = dir('./*');
        cd(strcat(repl(end).name,'/processed'));
        
        wdInitial = 32;
        overlap = 0.5;
        wdw = wdInitial*overlap;
        mmpropix = 1/140;
        sec = 0.04;

        threshold = 0.11;
        % Pictures to choose and cut around aggregate
        indi = [1:100];

        m = [1:100];
        indi(m)
        % Mask the Aggregate
        mkdir('results')
        [mask, centro, maxMin, area, perimeter, angle, xpl, ypl, xSm, ySm] = ... 
                maskAgg('../',indi(1),threshold,10);

        par.mmpropix = mmpropix;
        par.wdw = wdw;
        par.overlap = overlap;
        par.sec = sec;
        par.wdw = wdw;    
        par.centro = centro * mmpropix;
        par.area = area * mmpropix^2;
        par.angle = round(radtodeg(angle));
        par.perimeter = perimeter * mmpropix;
        par.esd = 2 * sqrt(area/pi) * mmpropix;

        par.lenX = (maxMin(2)-maxMin(1)) * mmpropix;
        par.lenY = (maxMin(4)-maxMin(3)) * mmpropix;
        fprintf('\n\n Aggregate area is %.6g mm^2\n',par.area)
        fprintf('\n Equivalent sphere diameter is %.6g mm\n',par.esd)
        fprintf('\n Aggregate perimeter is %.6g mm\n',par.perimeter)
        fprintf('\n Degrees of principle axis %g° %g°\n',par.angle)
        fprintf('\n Length along x-axis is %.6g mm\n',par.lenX)
        fprintf('\n Length along y-axis is %.6g mm\n\n',par.lenY)

        fact = 0.8;
        fact2 = 2;
        plot([(par.centro(1)+par.esd*fact)/mmpropix (par.centro(1)+par.esd*fact+1)/mmpropix], ...
            [(par.centro(2)+par.esd*fact)/mmpropix (par.centro(2)+par.esd*fact)/mmpropix],'w-','LineWidth',4)
        han = text((par.centro(1)+par.esd*fact+0.5)/mmpropix,(par.centro(2)+par.esd*fact)/mmpropix+26,'1 mm');
        set(han,'Color','white','HorizontalAlignment','center','FontSize',14)
        xlim([(par.centro(1)-par.esd*fact2)/mmpropix, (par.centro(1)+par.esd*fact2)/mmpropix])
        ylim([(par.centro(2)-par.esd*fact2)/mmpropix, (par.centro(2)+par.esd*fact2)/mmpropix])
        fig = gcf;
        fig.InvertHardcopy = 'off';
        
        aggfile = strcat('../',mscfolders(n).name,'-agg',aggfolders(agg_m).name,'-repl',repl(end).name);
        
        %print('./results/aggregateTracking','-dpng')
        cd ../../..
        %print(strcat(aggfile,'.png'),'-dpng')
        saveas(fig,strcat(aggfile,'.png'))
        pause(0.1)
        
        newmat = [xpl,ypl,xSm,ySm];
        save(strcat(aggfile,'.csv'),'newmat','-ASCII')
    end
    
    cd ..
end
