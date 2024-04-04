function [pivdata] = readPIVView(fname)
%% [x, y, Ux, Uy, magnitude] = readPIVView(fname)
% reads Pivview data
%
    file = fopen(fname,'r');

    n = 1;
    while n
        temp = fgets(file);
        if ~strcmp(temp(1),'#')
           n=0; 
        end
    end
    pivdata = {};
    
    temp = strrep(temp,'"','''');
    eval(temp)
    
    pivdata.title = TITLE;
    
    temp = fgets(file);
    temp = [temp(1:12) '{' temp(13:end) '}'];
    temp = strrep(temp,'"','''');
    eval(temp);

    temp = fgets(file);
    temp = temp(5:18);
    eval(temp);
    
    pivdata.I = I;
    pivdata.J = J;
    pivdata.variables = VARIABLES;
    pivdata.N = length(VARIABLES);
    
    mat = fscanf(file,'%g',[pivdata.N I*J]);

    for n = 1:length(VARIABLES);
       pivdata.(VARIABLES{n}) = reshape(mat(n,:),I,J)';
    end
    
    %x = reshape(mat(1,:),I,J)';
    %y = reshape(mat(2,:),I,J)';
    
    %Ux = flipud(reshape(mat(4,:),i,j)');
    %Uy = flipud(reshape(mat(5,:),i,j)');

    %pivdata.magnitude = flipud(sqrt(Ux.^2+Uy.^2));


    fclose(file);
end