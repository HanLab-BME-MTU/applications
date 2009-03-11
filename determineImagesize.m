function [data] = determineImagesize(data);
% determine image size of all fields in the experiment

for i=1:length(data)
    
    path  = data(i).source;
    od = cd;
    
    cd(path);
    cd('images283');
    loadfile = imread('cmask3001.jpg');
    imagesize = size(loadfile);
    
    data(i).imagesize = imagesize;
    cd(od);
end

end % of function