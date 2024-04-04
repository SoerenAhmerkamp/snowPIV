function [posX, posY] = particleTracking(fname,indi,maxMin,mask)
  
dX = [];
dY = [];
maxMin = round(maxMin);

Str = sprintf([fname 'agg8_highRe_2017-06-09-023234-%04g.tif'],indi(1));   
img = double(imread(Str));
 
h = waitbar(0);
%imgRef = medfilt2(img,[6 6]);
imgRef = img;
imgRef(~mask) = 0;
for n = 1:length(indi)-1

    %imgRef = medfilt2(img,[4 4]);
    Str = sprintf([fname 'agg8_highRe_2017-06-09-023234-%04g.tif'],indi(n+1));
    img = double(imread(Str,'tif'));
    %binMat = zeros(size(imgRef));
    
    %binMat(maxMin(3)-20:maxMin(4)+20,maxMin(1)-20:maxMin(2)+20) = 1;
    %imgRef = imgRef.*binMat;
    %imgRef = medfilt2(imgRef,[4 4]);
    %img = medfilt2(img,[4 4]);
    
    autocorrImg = fftshift(ifft2(fft2(imgRef).*fft2(rot90(imgRef,2))));
    xcorrImg = fftshift(ifft2(fft2(img).*fft2(rot90(imgRef,2))));

    [yMR xMR] = find(autocorrImg == max(autocorrImg(:)));
    [yM xM] = find(xcorrImg == max(xcorrImg(:)));

    dX(end+1) = xM-xMR;
    dY(end+1) = yM-yMR;
    
    waitbar(n / length(indi),h,'Tracking Particle');

end

delete(h)
posX = [0 dX];
posY = [0 dY];

%plot(posX,posY);










