
function [] = pcMD2tifs(MD,params,dirs)

% Assaf Zaritsky, Oct. 2017

if ~isfield(params,'sTime')
    params.sTime = 1;
end

for t = params.sTime : params.nTime
    imgFname = [dirs.images sprintf('%03d.tif',t)];        
    
    if exist(imgFname,'file') && ~params.always      
        continue;
    end
            
    I = MD.getChannel(1).loadImage(t);    
    
    imwrite(I,imgFname);    
end
end