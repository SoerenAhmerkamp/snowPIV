% PIV Processing

data = load('proc.txt')
%load('SCmap')

xG = data(:,1);
yG = data(:,2);
Ux = data(:,4);
Uy = data(:,5);
Uz = data(:,6);
U = data(:,7);

N = length(xG)./size(unique(xG));
Nx = N(1);
Ny = length(yG)./N(1);

pxtocm = 6 / 2176;

x = xG(1:Ny);
y = yG(1:Ny:length(yG));

U = ((reshape(U,[Ny,Nx])));
U = U';
% Interp
xInt = linspace(min(x),max(x),1000);
yInt = linspace(min(y),max(y),1000);
[xIntG yIntG] = meshgrid(xInt,yInt);
zInt = interp2(x,y,U,xIntG,yIntG);


fvel = pxtocm*10;
%
hold on
colormap(SCmap)
imagesc(xInt*pxtocm,yInt*pxtocm,zInt*fvel)
quiver(xG(1:4:end)*pxtocm,yG(1:4:end)*pxtocm,Ux(1:4:end),Uy(1:4:end),1.2,'w')
axis tight,axis equal
caxis([0,30].*fvel)
colorbar, 

whitebg(1,'k')
box on
set(gca,'color','black')
set(gcf,'color','black')
c = colorbar;
c.Color = 'w';
% Resolution ~ 0.100 mm 

fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)

