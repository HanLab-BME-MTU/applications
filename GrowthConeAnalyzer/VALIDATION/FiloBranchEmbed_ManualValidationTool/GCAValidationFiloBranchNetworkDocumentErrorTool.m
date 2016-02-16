function [ output_args ] = GCAValidationFiloBranchNetworkDocumentErrorToolFinal(reconstructDir,varargin)
% GCAValidationFiloBranchNetworkDocumentErrorTool
%
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('reconstructDir');
ip.addParameter('type','filoBranch');

ip.parse(reconstructDir,varargin{:});

type = ip.Results.type;

x = dir(reconstructDir);
x = x(3:end);
list = arrayfun(@(i) x(i).name,1:length(x),'uniformoutput',0);

% will eventually have a filter here that will check if done and only load
% those not finished if user desires.


idxInclude  = listSelectGUI(list,[],'move');


% listSelectGUI to select a file to work on
toTest = list(idxInclude);

for iProj = 1: numel(toTest)
    % make documentation directory
%     overlayDirC = [reconstructDir filesep toTest{iProj} filesep  'ValidationOverlays' filesep ];
%     
%     if ~isdir(overlayDirC)
%         mkdir(overlayDirC)
%     end
    
    switch type
        case 'filoBranch'
            listOfImages =  searchFiles('.png',[],[reconstructDir filesep toTest{iProj} filesep ...
                'OverlayFiloBranch' filesep 'Reconstruct_Movie'],0,'all',1);
              
            overlayDirC = [reconstructDir filesep toTest{iProj} filesep 'OverlayFiloBranch' filesep ... 
                'ValidationOverlays' filesep ];
            
            
            
        case 'embed'
            listOfImages = searchFiles('.png',[],[reconstructDir filesep toTest{iProj} filesep ...
                'OverlayEmbed'],0,'all',1);
              
            overlayDirC = [reconstructDir filesep toTest{iProj} filesep 'OverlayEmbed' filesep ... 
                'ValidationOverlays' filesep ];
            
       
            
            
    end
    
    if ~isdir(overlayDirC)
        mkdir(overlayDirC);
    end
    
    
    imgRaw = imread(listOfImages{1});
    imgOverlay = imread(listOfImages{end});
    imgLarge = [imgOverlay imgRaw];
    [nyLarge,nxLarge,~] = size(imgLarge);
    
    % move to next non-documented frame?
    
    % choose a file manually?
    
    %
    % load the chosen project directory
    
    errorName{1} = 'False Pos'; 
    errorName{2} = 'False Neg'; 
    errorName{3} = 'Mis'; 
    
    colors{1} = 'b'; 
    colors{2} = 'r'; 
    colors{3} = 'g'; 
    
    for iError = 1:3
        % The first and last frames in the Troubleshoot Reconstruction Directory
        % will be what we load.
        [ny,nx,~] = size(imgRaw);
        
        
        setFigure(nxLarge,nyLarge,'on');
        imshow(imgLarge,[]);
        hold on
        coords =  gcaValidationDocument(errorName{iError},[ny,nx],colors{iError});
        name = strrep(errorName{iError},' ' ,'_'); 
        save([overlayDirC filesep name '_Coords.mat'],'coords');
        saveas(gcf, [overlayDirC filesep name '_Overlay.fig']);
        saveas(gcf,[overlayDirC filesep name '_Overlay.png']); 
        saveas(gcf,[overlayDirC filesep name '_Overlay.eps'],'psc2'); 
        close gcf
        clear coords
    end
  
    
end % iProj



