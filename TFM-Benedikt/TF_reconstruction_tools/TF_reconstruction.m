% This program (regularized fourier transform traction force
% reconstruction) was produced at the University of Heidelberg, BIOMS
% group of Ulich Schwarz. It calculates traction from a gel displacement
% field.
%
% Benedikt Sabass 20-5-2007

function TF_reconstruction
    where_am_I = mfilename('fullpath');
    [my_directory,name] = fileparts(where_am_I);
    Hansen_sub = fullfile(my_directory,'Hansen');
    if isempty(strfind(path,my_directory))
        path(path, my_directory);
    end
    if isempty(strfind(path,Hansen_sub))
        path(path, Hansen_sub);
    end
    get_data;
end