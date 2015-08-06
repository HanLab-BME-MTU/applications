function [ h ] = protrusionMapMine( imgData,pixSizeNm,frameInSec,umPerMinFlag,CLim,blackOutHi,blackOutLo)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here 
% if nargin<3 
%     makeMovie =0 ;
% end 
%setFigureMoviePlots(200,200,'off')

h = setAxis;



Colormap='jet';
x = colormap;
if blackOutLo == 1
% colormap set to go from lowest to highest
x(1,:) = [1,1,1]; % set the original data to white
end 
if blackOutHi ==1
x(end,:) = [0,0,0]; % set the outlier data to black
end 
colormap(x)


colorbar


% convert imageData to nm per sec
                imgData = imgData.*pixSizeNm; % currently 216 nm per pixel % convert to um
                if umPerMinFlag ==1 
                    imgData = imgData./1000; 
                end 
                                             
                imgData= imgData./frameInSec;  % 5 sec per frame ( convert to minutes)
                
                if umPerMinFlag ==1 
                    imgData = imgData.*60; 
                end 
                
imagesc(imgData,CLim);



axis xy
% white out NaN values
%                              alphamask =true(size(imgData));
%                              alphamask(isnan(imgData))=false;
%                              set(h,'AlphaData',alphamask,'AlphaDataMapping','none');
ylabel('Window Number');

if umPerMinFlag ==1
forLabel = 'min'; 
else 

forLabel = 'sec'; 
end 

xlabel('Time (s)' )

% if umPerMinFlag ==1 
%     converted = converted./60;  
% end 
%end 
 
 

%set(gca,'FontSize',10,'FontName','Arial');
colorbar
[ny,nx] = size(imgData); 
axis([1,nx-1,1,ny]); 
%  title(name,sfont{:});
frameNum = get(gca,'XTickLabel');
%frameInMin = frameInSec*60; 

converted = str2double(frameNum).*frameInSec; % SET THE TIME INTERVAL
num2str(converted);
set(gca,'XTickLabel',converted);
%else 
% if makeMovie~=1
% saveas(gcf,[saveDir filesep filename '.eps'],'psc2');
% saveas(gcf,[saveDir filesep filename '.fig']);
% saveas(gcf,[saveDir filesep filename '.png']);
%end % else add scatters
end

