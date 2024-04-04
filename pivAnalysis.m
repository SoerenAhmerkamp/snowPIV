%% Set Parameter
clear all;

load Calibration.txt
%% cd C:\Users\admin\Desktop\piv_processing\pivclara\5\processed2
% Set Parameter

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
%mkdir('results')
[mask, centro, maxMin, area, perimeter, angle, xpl, ypl, xSm, ySm, Feret] = ... 
        maskAgg('./',indi(1),threshold/1,10);

par.feret = Feret*mmpropix;
par.mmpropix = mmpropix;
par.wdw = wdw;
par.overlap = overlap
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
fprintf('\n Ferret diameter is %.6g mm\n\n',par.feret)

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
%print('./results/aggregateTracking','-dpng')
%saveas(fig,'./results/aggregateTracking.png')

save D:\ClaraPIVData\DataExport_properties\dataMSC06_9_3 par

%% Particle tracking (based on cross correlation)
[posX, posY] = particleTracking('../',indi,maxMin,mask)
par.posX = posX;
par.posY = posY; 

hold on
plot(par.posX);
plot(par.posY);
legend('X-displacement','Y-displacement')
print('./results/aggregateDisplacement','-dpng')

%% Calculate the cutoff values (determined by the average of the picture)
cutYneg = centro(2) - max(-posY);
cutYpos = size(mask,1) - centro(2) - max(posY);
cutY = [cutYneg cutYpos];

cutXneg = centro(1) + min(posX);
cutXpos = size(mask,2) - centro(1) - max(posX);
cutX = [cutXneg cutXpos];

%% Test particle tracking
imgSave = 0;
cutAgg = round((maxMin(2)-maxMin(1))/2)+200;
cutAggy = round((maxMin(4)-maxMin(3))/2)+200;
% The following routine is just implemented to test the particle tracking
% How to adjust indi to test x- and y-displacement /Clara
for n = m
    Str = sprintf(['../' 'Frame_%07g.bmp'],indi(n));   
    img = double(imread(Str));
    %n
    centro = round(centro);
    xVal = centro(1)-cutAgg:centro(1)+cutAgg;
    yVal = centro(2)-cutAggy:centro(2)+cutAggy;
    imgTemp = img(yVal+posY(n),xVal+posX(n));
    imgSave = imgTemp + imgSave;
    subplot(1,2,1)
    imagesc(img(yVal+posY(n),xVal+posX(n))./max(img(:))), axis equal, axis tight
    caxis([0 0.2])
    subplot(1,2,2)
    imagesc(img(yVal,xVal)./max(img(:))), axis equal, axis tight
    pause(0.2)
end
figure;
imshow(uint8(imgSave./length(indi(m))))
% Less than 8 px should be fine
correction = 1;

print('./results/aggregate','-dpng')
% plot(posX,posY,'o')
%% Calculate Average flow fields
% 
close all
clear avgMagnitude avgUx avgUy xVal0 yVal0 avgNumb centroPiv
cutXPiv = round(cutX ./ wdw);
cutYPiv = round(cutY ./ wdw);
posXPiv = round(posX ./ wdw);
posYPiv = round(posY ./ wdw);

conv = [];

cutYPiv(2) = cutYPiv(2) - 2;
cutYPiv(1) = cutYPiv(1) - 1;

cutXPiv(1) = cutXPiv(1) - 1;
cutXPiv(2) = cutXPiv(2) - 2;

for n = m
    n
    Str = sprintf('Frame_%07g_piv.dat',indi(n)); % indi(n)
    disp(['Processing... ', Str])
    pivdata = readPIVView(Str);
    
    magnitude = medfilt2(pivdata.mag,[3, 3]);
    Ux = medfilt2(pivdata.dx,[3,3]);
    Uy = medfilt2(pivdata.dy,[3,3]);
    
    y = pivdata.iy;
    y = y - mod(min(y(:)),wdw);
    corr = flipud(pivdata.corr);
    corrBin = corr>10;
    magnitude = flipud(magnitude);
    Ux = flipud(Ux);
    Uy = flipud(Uy);
    x = pivdata.ix;
    
    [xTemp,yTemp] = find(x==round(centro(1)/wdw)*wdw & y == round(centro(2)/wdw)*wdw);
    centroPiv = [yTemp,xTemp];
    
    xVal = centroPiv(1)-cutXPiv(1)+posXPiv(n):centroPiv(1)+cutXPiv(2)+posXPiv(n);
    yVal = centroPiv(2)-cutYPiv(1)+posYPiv(n)+1:centroPiv(2)+cutYPiv(2)+posYPiv(n);
    magnitude = magnitude(yVal,xVal);
    Ux = Ux(yVal,xVal);
    Uy = Uy(yVal,xVal);
    x = x(yVal,xVal);
    y = y(yVal,xVal);
    corrBin = corrBin(yVal,xVal);
        
    if ~exist('avgMagnitude')
        xVal = centroPiv(1)-cutXPiv(1)+posXPiv(1):centroPiv(1)+cutXPiv(2)+posXPiv(1);
        yVal = centroPiv(2)-cutYPiv(1)+posYPiv(1)+1:centroPiv(2)+cutYPiv(2)+posYPiv(1);
        centroPiv0 = -[-cutXPiv(1)+posXPiv(n) -cutYPiv(1)+posYPiv(n)];
        xVal0 = xVal;
        yVal0 = yVal;
        avgMagnitude = zeros(size(magnitude));
        avgNumb = zeros(size(magnitude));
        avgUx = zeros(size(Ux));
        avgUy = zeros(size(Uy));
    end
    
    %sum(pivdata.snr(:))./(pivdata.I*pivdata.J)
    
    diff = avgMagnitude.^2 ./ n^.2 - magnitude.^2;
    conv(end+1) = sum(diff(:))./sum(avgMagnitude(:).^2);
   
    avgNumb = avgNumb + corrBin;
    magnitude(~corrBin) = 0;
    Ux(~corrBin) = 0;
    Uy(~corrBin) = 0;
    avgMagnitude = avgMagnitude + magnitude;
    avgUx = avgUx + Ux;
    avgUy = avgUy + Uy;
    
    subplot(1,2,1)
    hold on
    pcolor(x,-y,avgMagnitude./avgNumb), shading flat, axis equal, axis tight
    xlabel('x in pixel')
    ylabel('y in pixel')
    title('Averaged')
    plot(centro(1),-centro(2),'wx')
    
    subplot(1,2,2)
    hold on
    pcolor(x,-y,magnitude), shading flat, axis equal, axis tight
    caxis([0,12])
    xlabel('x in pixel')
    ylabel('y in pixel')
    title('Instantaneous')
    plot(centro(1),-centro(2),'wx')
    
    % Use % to suppress image saving
    %print(['./results/pivResultate' num2str(indi(n))],'-dpng')

    pause(0.5)
end

avgMagnitude = avgMagnitude ./ avgNumb;
avgUx = avgUx ./ avgNumb;
avgUy = avgUy ./ avgNumb;

par.centroPiv0 = centroPiv0;

par.posXm = posX(m);
par.posYm = posY(m);

fprintf('\n\n Maximum particle movement in x-direction: %.6g px\n',range(par.posXm))
fprintf('\n\n Maximum particle movement in y-direction: %.6g px\n\n',range(par.posYm))

print('./results/pivResultateAvg','-dpng');

%% Mask Aggregate
close all
xS = wdw:wdw:((size(mask,1)-mod(size(mask,1),wdw))-wdw);
yS = wdw:wdw:((size(mask,2)-mod(size(mask,2),wdw))-wdw);
maskT = mask(xS,yS);

maskT = maskT(yVal0,xVal0);

avgMagnitude(maskT==1) = NaN;
avgUx(maskT==1) = NaN;
avgUy(maskT==1) = NaN;

par.avgMagnitude = avgMagnitude * mmpropix / sec * 3600 / 1000 * 24;
par.avgUx = avgUx * mmpropix / sec * 3600 / 1000 * 24;
par.avgUy = avgUy * mmpropix / sec * 3600 / 1000 * 24;

par.x = x * mmpropix;
par.y = y * mmpropix;

figure
pcolor(par.x,-par.y,par.avgUy), shading flat, axis equal, axis tight
colormap('jet'), colorbar
xlabel('x in mm')
ylabel('y in mm')
title('y - Velocity in m / d')
print('./results/pivResultateMaskUy','-dpng')

figure
pcolor(par.x,-par.y,par.avgUx), shading flat, axis equal, axis tight
colormap('jet'), colorbar
xlabel('x in mm')
ylabel('y in mm')
title('x - Velocity in m / d')
print('./results/pivResultateMaskUx','-dpng')

figure
pcolor(par.x,-par.y,par.avgMagnitude), shading flat, axis equal, axis tight
colormap('jet'), colorbar
xlabel('x in mm')
ylabel('y in mm')
title('Velocity-Magnitude in m / d')
print('./results/pivResultateMaskMag','-dpng')
%%

xNorm = (par.x-par.centro(1))./par.esd;
yNorm = (par.y-par.centro(2))./par.esd;

figure
pcolor(xNorm,-yNorm,par.avgMagnitude), shading flat, axis equal, axis tight
colormap('jet'), colorbar
xlabel('x in mm')
ylabel('y in mm')
title('x - Velocity in m / d')


%% Boundary Layer Calculations and velocity field
% second number is distance to aggregate center
% third number is alignment (1=horizontal, 2=vertical)
% last number is cut-off value for DBL- and Umax-calculations (to be adjusted)

par.Slice_names = {};

[par] = profileCalc(par,0,1,'SliceEquator',3);
print(char(strcat('./results/',par.Slice_names{end})),'-dpng')

[par] = profileCalc(par,3,1,'SliceFront',4);
print(char(strcat('./results/',par.Slice_names{end})),'-dpng')

[par] = profileCalc(par,0,2,'SliceVertical',4.5);
print(char(strcat('./results/',par.Slice_names{end})),'-dpng')


%% Write into xls
save ./results/data.mat
piv2XLS(par,'./results/output.xlsx')


%% 
figure
startx = 4:0.04:9;
starty = 2*ones(size(startx));

%[xG,yG] = meshgrid(min(par.x(:)):0.001:max(par.x(:)),min(par.y(:)):0.1:max(par.y(:)));
%avgGUx = griddata(par.x,par.y,par.avgUx,xG,yG);
%avgGUy = griddata(par.x,par.y,par.avgUy,xG,yG);
%avgGmagnitude = griddata(par.x,par.y,par.avgMagnitude,xG,yG);
xI = 1:size(img,2)*par.mmpropix;
yI = 1:size(img,1)*par.mmpropix;

hold on
%imagesc(xI-1,yI,flipud(img)), colormap('gray')
pcolor(par.x,par.y,flipud(par.avgMagnitude)),shading flat,axis equal, axis tight

h = streamline(stream2(unique(par.x),unique(par.y),par.avgUx,par.avgUy,startx,starty));
%h = streamline(stream2(unique(xG),unique(yG),avgGUx,avgGUy,startx,starty));
set(h,'color','w')

%while 1
%    [startx,starty] = ginput(1);
%    h = streamline(mmstream2(unique(par.x),unique(par.y),par.avgUx,par.avgUy,startx,starty));
%    set(h,'color','w')
%end





