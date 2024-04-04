function [mask, centro, maxMin, area, perimeter, angle, xplot, yplot, xSm, ySm, Ferret] = maskAgg(fname,indi,ths,smPar)
%% [x, y, Ux, Uy, magnitude] = readPIVView(fname)
% reads Pivview data
%

if nargin > 1
    len = length(indi);
    for n = indi
        Str = sprintf([fname 'Frame_%07g.bmp'],n);
        if exist('img','var')
            img = img + double(imread(Str));
        else
            img = double(imread(Str));
        end
    end
    img = img / len;
else
    img = imread(fname);
end

yS = size(img,1);
xS = size(img,2);

[xv yv] = meshgrid(1:xS,1:yS);

xv = xv(:);
yv = yv(:);

imshow(uint8(img))


imgBin = im2bw(medfilt2(img,[3 3])./max(img(:)),ths);
CC = bwconncomp(imgBin);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
imgBin(:) = 0;
values = CC.PixelIdxList{idx};
imgBin(values) = 1;
[B,mask] = bwboundaries(imgBin,'noholes');
%plot(xv(values),yv(values),'.')
y = B{1}(:,1);
x = B{1}(:,2);

%pause
%[x, y] = getline;
[geom,~,cpmo] = polygeom(x,y);
centro = [geom(2) geom(3)];
area = geom(1);
perimeter = geom(4);
angle = [cpmo(2) cpmo(4)];


%mask = inpolygon(xv,yv,x,y);
maxMin = [min(x) max(x) min(y) max(y)];
xv = reshape(xv,[yS xS]);
yv = reshape(yv,[yS xS]);
mask = reshape(mask,[yS xS]);

% show results
x_cen = geom(2);
y_cen = geom(3);

I1 = cpmo(1);
angle1 = cpmo(2);
I2 = cpmo(3);
angle2 = cpmo(4);
 
% plot outline
xplot = x( [ 1:end 1] );
yplot = y( [ 1:end 1] );
rad = 250;
x1 = [ x_cen-rad*cos(angle1)  x_cen+rad*cos(angle1) ];
y1 = [ y_cen-rad*sin(angle1)  y_cen+rad*sin(angle1) ];
x2 = [ x_cen-rad*cos(angle2)  x_cen+rad*cos(angle2) ];
y2 = [ y_cen-rad*sin(angle2)  y_cen+rad*sin(angle2) ];
hold on
plot( xplot,yplot,'b', x_cen,y_cen,'ro', ...
      x1,y1,'g:', x2,y2,'r:'  )

[xi,yi] = polyxpoly(x1,y1,xplot,yplot);
plot(xi(1),yi(1),'r.','MarkerSize',20)
plot(xi(end),yi(end),'r.','MarkerSize',20)

%Ferret = norm([xi(1),yi(1)]-[xi(end),yi(end)]);

% Find which two boundary points are farthest from each other.
maxDistance = zeros(n,length(xplot));

for k = 1 : length(xplot)
    distances = sqrt( (xplot(k) - xplot) .^ 2 + (yplot(k) - yplot) .^ 2 );
    [maxDistance(k)] = max(distances);
end
%figure, plot(maxDistance,'o')
%maxDistance
Ferret = max(maxDistance);

xSm = smooth(xplot,smPar);
ySm = smooth(yplot,smPar);
plot(xSm,ySm,'w') 
end