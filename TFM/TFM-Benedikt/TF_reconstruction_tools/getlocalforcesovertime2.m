function [localtraction,averaged_localtraction,contraction] = getlocalforcesovertime2(TFM_results, refbrightness)
    [FileName,PathName] = uigetfile('*.*','Select image file to locate region');
    bild = imread(fullfile(PathName,FileName));
%      bild = im2double(bild);
    h = figure; imagesc(bild); colormap gray; hold;
%    plot(TFM_results(1).pos(:,1),TFM_results(1).pos(:,2),'.y'); 
    [horpos,vertpos] = meshgrid(1:size(bild,1),1:size(bild,2));
%    if nargin < 2
%        title('Select region for reference');
%        b = getline('closed');
%        plot(b(:,1),b(:,2),'r'); 
%        
%        insideref = inpolygon(vertpos, horpos,b(:,1),b(:,2));
%    end
    
    title('Select region of interest');
    a = getline('closed');
    plot(a(:,1),a(:,2),'g');
    
    hold off;
    dir_struct = dir(fullfile(PathName,'*.tif'));
    [sorted_names,sorted_index] = sortrows({dir_struct.name}');

   
    insidebild = inpolygon(vertpos, horpos,a(:,1),a(:,2));
  
        
    for frame = 1:size(TFM_results,2)
 %       imname = strcat(PathName, filesep, sorted_names(frame));
 %       bild = imread(imname{1});
        
        inside = inpolygon(TFM_results(frame).pos(:,1),TFM_results(frame).pos(:,2),a(:,1),a(:,2));
        localtraction(frame) = mean(TFM_results(frame).traction_magnitude(inside,:));
        localtraction_vec(frame).data = TFM_results(frame).traction_magnitude(inside,:);
%        brightness(frame) = sum(sum(bild.*uint16(insidebild')))./nnz(insidebild);
%        if nargin < 2
%            refbrightness(frame) = sum(sum(bild.*uint16(insideref')))./nnz(insideref);
%        end
    end
 
    maxaverage = floor(size(TFM_results,2)/3)*3;
    %{
    for frame = 1:3:maxaverage-3
        averaged_localtraction(frame) = mean(localtraction(frame:frame+3));
    end
%}
    averaged_localtraction = [];
%    normbrightness = brightness-refbrightness;
   
    contraction = [];
    for frame = 1:size(TFM_results,2)
        contraction = vertcat(contraction,localtraction_vec(frame).data);
    end
    figure;
    
    cla;hold on;
    plot(1:size(TFM_results,2),localtraction,'r-*');
    plot(1:3:length(averaged_localtraction),averaged_localtraction(1:3:end),'b-*'); hold off;
    title('Mean Traction in selected Region over time');
 
%    figure;
%    cla;hold on;
%    plot(1:size(TFM_results,2),normbrightness,'g-*'); hold off;
%    title('Mean Brightness in selected Region over time');
%    
%    disp('Forces extracted');
end    