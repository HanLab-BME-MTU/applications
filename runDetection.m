function [exp] = runDetection(exp);
% runDetection automatically runs Henry's detection software - with iclean
% 1 on all the movies specified in the experiment structure
% SYNOPSIS [exp] = loadAndSaveDetection(exp);
%
% INPUT     exp    :    experiment structure, which has to contain a field
%                       .source
%                       which is the path to the data location
% OUTPUT               
% REMARKS   The function only performs the detection for a given movie if 
%           there are no detection data (folder called maxdata) already 
%           present in the specified directory; if you want a partial or
%           faulty detection to be replaced, you need to DELETE the folders
%           first!!!!
%
% Dinah Loerke, last modified Feb 2008

od = cd;

% loop over all entries in the structure to enter the image data necessary
% for the detection input
for i=1:length(exp)
    
    path = exp(i).source;
    % change to source directory
    cd(path);
    
    [oriImageName, oriImagePath] = uigetfile('.tif',['Select first original image for movie #',num2str(i)]); %HRJ'.tif'
    
    imData(i).oriImageName = oriImageName;
    imData(i).oriImagePath = oriImagePath;
    exp(i).firstImageName = oriImageName;
    
    cd(od);
end


for i = 1:length(exp)
    
    path = exp(i).source;
    % change to source directory
    cd(path);
    
    
    oriImageName = imData(i).oriImageName;
    oriImagePath = imData(i).oriImagePath;
    
    runvar = 0;
    % if detection results folder already exist, skip this entry 
    if  ~ (exist('maxdata283')==7) 
        runvar = 1;
    end
    
    if runvar==1
        main283AUTO(1,oriImageName,oriImagePath);
    end
    
    cd(od);
    
    
end % of for i

end % of function