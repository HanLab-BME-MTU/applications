 function [ output_args ] = GCAValidationFiloBranchNetworkDocumentErrorTool(reconstructDir,varargin)
% GCAValidationFiloBranchNetworkDocumentErrorTool
%
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('reconstructDir');
defaultType{1} = 'VeilFilo';
defaultType{2} = 'Branch';
defaultType{3} = 'Embed';
ip.addParameter('type',defaultType); % embed, veilFilo, branch


ip.parse(reconstructDir,varargin{:});
%%
typesAll = ip.Results.type;

x = dir(reconstructDir);
x = x(3:end);
list = arrayfun(@(i) x(i).name,1:length(x),'uniformoutput',0);

% will eventually have a filter here that will check if done and only load
% those not finished if user desires.


idxInclude  = listSelectGUI(list,[],'move');


% listSelectGUI to select a file to work on
toTest = list(idxInclude);

     
colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'g';

for iType = 1:numel(typesAll)
    
    typeC = typesAll{iType};
    
    
    
    
    for iProj = 1: numel(toTest)
        
        errorName{1} = ['False Pos ' typeC ];
        errorName{2} = ['False Neg ' typeC ] ;
        
        if strcmpi(typeC, 'branch')
            errorName{3} = ['Misconnection ' typeC];
        end
        
        
        listOfImages =  searchFiles('.png',[],[reconstructDir filesep toTest{iProj} filesep ...
            'Overlay' typeC ],0,'all',1);
        
        overlayDirC = [reconstructDir filesep toTest{iProj} filesep 'Overlay' typeC filesep ...
            'ValidationOverlays' filesep ];
        
        
        if ~isdir(overlayDirC)
            mkdir(overlayDirC);
        end
        
        
        imgRaw = imread(listOfImages{1});
        imgOverlay = imread(listOfImages{end});
        imgLarge = [imgOverlay imgRaw];
        [nyLarge,nxLarge,~] = size(imgLarge);
        
          
        for iError = 1:numel(errorName)
            % The first and last frames in the Troubleshoot Reconstruction Directory
            % will be what we load.
            [ny,nx,~] = size(imgRaw);
            
            
            setFigure(nxLarge,nyLarge,'on');
            imshow(imgLarge,[]);
            hold on
            [coords,filoCount] =  gcaValidationDocument(errorName{iError},[ny,nx],colors{iError});
            name = strrep(errorName{iError},' ' ,'_');
            save([overlayDirC filesep name '_Coords.mat'],'coords','filoCount');
            saveas(gcf, [overlayDirC filesep name '_Overlay.fig']);
            saveas(gcf,[overlayDirC filesep name '_Overlay.png']);
            saveas(gcf,[overlayDirC filesep name '_Overlay.eps'],'psc2');
            close gcf
            clear coords
        end
         clear errorName
        
    end % iProj
    
end % iType

