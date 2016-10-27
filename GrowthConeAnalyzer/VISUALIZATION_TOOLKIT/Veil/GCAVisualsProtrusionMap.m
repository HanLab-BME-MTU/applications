function [ h ] = GCAVisualsProtrusionMap(imgData, varargin)
%GCAVisualsProtrusionMap 

%%Input check
ip = inputParser;

ip.CaseSensitive = false;


ip.addRequired('imgData'); % what windowing method to use - currently assumes you will use one of the
% windowing methods associated with the MD.
% defaultOutDir =  [movieData.outputDirectory_ filesep 'VisualizationOverlays' filesep 'WholeNeurite' filesep 'VeilWindows'];
% 
% ip.addParameter('OutputDirectory_', defaultOutDir);

% ip.addParameter('pixSizeNm',0.218); 

ip.addParameter('frameInSec',5); 

ip.addParameter('umPerMin',false); % default assumed in nm/sec  

ip.addParameter('CLim',[-100,100]); 

ip.addParameter('blackOutHi',false); 

ip.addParameter('blackOutLo',false); 

ip.addParameter('colorMap','jet'); 

ip.addParameter('FontSize',20); 

ip.addParameter('ShowColorBar',true); 

ip.addParameter('alphaMaskNaN',false); 

ip.addParameter('title',[]); 
ip.addParameter('visible','off'); 

ip.parse(imgData,varargin{:});

%% 
% pixSizeNm = ip.Results.pixSizeNm; 
 frameInSec = ip.Results.frameInSec; 

%setFigureMoviePlots(200,200,'off')

h = setAxis(ip.Results.visible,0.95,ip.Results.FontSize);
x = ip.Results.colorMap;
if ip.Results.blackOutLo
    % colormap set to go from lowest to highest
    x(1,:) = [1,1,1]; % set the original data to white
end
if ip.Results.blackOutHi
    x(end,:) = [0,0,0]; % set the outlier data to black
end
colormap(x)
% convert imageData to nm per sec
%                 imgData = imgData.*pixSizeNm; % currently 216 nm per pixel % convert to um
%                 if ip.Results.umPerMin 
%                     imgData = imgData./1000; 
%                 end 
%                                              
%                % imgData= imgData./frameInSec;  % 5 sec per frame ( convert to minutes)
%                 
%                 if ip.Results.umPerMin ==1 
%                     imgData = imgData.*60; 
%                 end 
                
h = imagesc(imgData,ip.Results.CLim);

axis xy
% white out NaN values
%                              alphamask =true(size(imgData));
%                              alphamask(isnan(imgData))=false;
%                              set(h,'AlphaData',alphamask,'AlphaDataMapping','none');
ylabel('Window Number');
% 
if ip.Results.umPerMin
    forLabel = 'min';
else
    
    forLabel = 'sec';
end

xlabel(['Time ' '( ' forLabel ' )'] )

% if umPerMinFlag ==1 
%     converted = converted./60;  
% end 
%end 
 
 

%set(gca,'FontSize',10,'FontName','Arial');
if ip.Results.ShowColorBar
    colorbar
end
[ny,nx] = size(imgData); 
axis([1,nx-1,1,ny]); 
%  title(name,sfont{:});
frameNum = get(gca,'XTickLabel');
%frameInMin = frameInSec*60; 

converted = str2double(frameNum).*frameInSec; % SET THE TIME INTERVAL
num2str(converted);
set(gca,'XTickLabel',converted);

if ip.Results.alphaMaskNaN
    
 % mask out NaNs 
alphaMask = false(size(imgData));
alphaMask(isnan(imgData)) = true;
set(gcf,'AlphaMap', [1,0.8]);
% first number is will set the transparancy of
% of the true values in the alphaMask while the second number will set the
% transparancy of the false number of the alpha mask
set(h, 'alphaData', alphaMask,'alphaDataMapping','direct');
end 

if ~isempty(ip.Results.title)
    title(ip.Results.title); 
end 

%else 
% if makeMovie~=1
% saveas(gcf,[saveDir filesep filename '.eps'],'psc2');
% saveas(gcf,[saveDir filesep filename '.fig']);
% saveas(gcf,[saveDir filesep filename '.png']);
%end % else add scatters
end

