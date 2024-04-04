function piv2XLS(par,filename)
%%
xlswrite(filename,{'PIV Settings'},'Parameter','A1')
xlswrite(filename,{'Calibration','mm/px',par.mmpropix}','Parameter','A2:A4')
xlswrite(filename,{'Recording','sec',par.sec}','Parameter','B2:B4')
xlswrite(filename,{'Window Size','px',par.wdw/par.overlap}','Parameter','C2:C4')
xlswrite(filename,{'Window Overlap','',par.overlap}','Parameter','D2:D4')
xlswrite(filename,{'Final Window Size','px',par.wdw}','Parameter','E2:E4')


xlswrite(filename,{'Aggregate Parameter'}','Parameter','F1')
xlswrite(filename,{'Area','mm^2',par.area}','Parameter','F2:F4')
xlswrite(filename,{'Perimeter','mm',par.perimeter}','Parameter','G2:G4')
xlswrite(filename,{'Equivalent Sphere Diameter','mm',par.esd}','Parameter','H2:H4')
xlswrite(filename,{'Diameter x-axis','mm',par.lenX}','Parameter','I2:I4')
xlswrite(filename,{'Diameter y-axis','mm',par.lenY}','Parameter','J2:J4')
xlswrite(filename,{'Angles of principle axis1','\circ',par.angle(1)}','Parameter','K2:K4')
xlswrite(filename,{'Angles of principle axis2','\circ',par.angle(2)}','Parameter','L2:L4')

xlswrite(filename,{'Aggregate Displacement Range x','px',range(par.posXm)}','Parameter','M2:M4')
xlswrite(filename,{'Aggregate Displacement Range y','px',range(par.posYm)}','Parameter','N2:N4')

%xlswrite(filename,{'Flow Field Characteristics (Slices)'},'Parameter','O1')

N=1;
% 
for n = 1:length(par.Slice_names) 
    n
    xlswrite(filename,{char(par.Slice_names{n})},...
        'Parameter',[char(78+N) num2str(1)])
    
    xlswrite(filename,{'Umax','m/d',par.(char(strcat(par.Slice_names{n},'Umax')))}', ...
    'Parameter',[char(78+N) num2str(2) ':' char(78+N) num2str(4)])
    N = N+1;
    
    xlswrite(filename,{'Umax Right','m/d',par.(char(strcat(par.Slice_names{n},'UmaxRight')))}', ...
    'Parameter',[char(78+N) num2str(2) ':' char(78+N) num2str(4)])
    N = N+1;

    xlswrite(filename,{'Boundary Layer Thickness','mm',par.(char(strcat(par.Slice_names{n},'BoundaryThickness')))}', ...
    'Parameter',[char(78+N) num2str(2) ':' char(78+N) num2str(4)])
    N = N+1;
    
       xlswrite(filename,{'Boundary Layer Thickness Right','mm',par.(char(strcat(par.Slice_names{n},'BoundaryThicknessRight')))}', ...
    'Parameter',[char(78+N) num2str(2) ':' char(78+N) num2str(4)])
    N = N+1;

end




