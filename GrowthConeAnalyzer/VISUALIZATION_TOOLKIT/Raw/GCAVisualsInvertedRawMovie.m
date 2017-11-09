function [ success,saveDir] = GCAVisualsInvertedRawMovie(MD,varargin)
%GCAVisualsInvertedRawMovie
% small function to



%% Check input 
ip = inputParser; 
ip.CaseSensitive = false;

ip.addParameter('ScaleBar',true); 
ip.addParameter('Timer',true); 
ip.addParameter('TreatmentFrame',[]); 
ip.addParameter('TreatmentName','CK666'); 
defaultOut = [MD.outputDirectory_  filesep  'VisualizationOverlays' filesep 'Raw'];
ip.addParameter('OutputDirectory',defaultOut); 
ip.addParameter('Overwrite',false); 
ip.addParameter('NeuriteID',true); 
ip.addParameter('NeuriteIDText', []); % if isempty and true find the neurite ID otherwise use input
ip.addParameter('startFrame',1); 
ip.addParameter('endFrame',[]); 
 

ip.parse(varargin{:});

%% 
% for now just keep success to 1
success = 1;

run = 0; 

% Perform only if the directory does not exist OR if user wants to
% overwrite
% check if the directory exists : if yes skip
if ~isdir(ip.Results.OutputDirectory)
    mkdir(ip.Results.OutputDirectory)
    run = 1;
end
  
if ip.Results.Overwrite 
    run =1; 
end 


% check to make sure you have the files you need
names = MD.getImageFileNames{1};
if run == 1;
    if ~isempty(ip.Results.endFrame)
        endFrame = ip.Results.endFrame;
    else
        endFrame = numel(names);
    end
    for iFrame = ip.Results.startFrame:endFrame
        % load image
        
        img = double(imread([MD.getChannelPaths{1} filesep names{iFrame}]));
        [ny,nx] = size(img);
        setFigure(nx,ny,'off');
        imshow(-img,[]);
        hold on
        
        if ip.Results.ScaleBar == true
            pixSizeMic = MD.pixelSize_/1000; 
            width  = 10/pixSizeMic; 
            plotScaleBar(width,2,'Location','SouthWest','Color',[0 0 0]);
        end
        
        if ip.Results.Timer == true
            text(nx-35,ny-10,[num2str(iFrame*MD.timeInterval_ - MD.timeInterval_),' s'] ,'Color', 'k',...
                'FontSize',10,'FontWeight','Bold');
        end
        
        if ~isempty(ip.Results.TreatmentFrame)
            if iFrame<=ip.Results.TreatmentFrame
                text(10,25,'Pre-Treatment','FontSize',10,'Color','k');
            else
                text(10,25,['Post-Treatment: ' ip.Results.TreatmentName],'FontSize',10,'Color','k');
            end
            
        end
        
        
        if ip.Results.NeuriteID
            if ~isempty(ip.Results.NeuriteIDText)
                neuriteID = ip.Results.NeuriteIDText;
            else
                neuriteID = gcaGetNeuriteID(MD.outputDirectory_);
            end
            
            text(10,10,neuriteID , 'FontSize', 10, 'Color', 'k');
        else
            
        end
        
        
        saveas(gcf,[ip.Results.OutputDirectory filesep num2str(iFrame,'%03d') '.png']);
        saveas(gcf,[ip.Results.OutputDirectory filesep num2str(iFrame,'%03d') '.fig']); 
        
        saveas(gcf,[ip.Results.OutputDirectory filesep num2str(iFrame,'%03d') '.eps'], 'psc2'); 
        close gcf
    end % iFrame 
    
    %% use ffmpeg to make movie (on windows having some trouble with odd sizes) 
%     cd(['C:\Users\Maria\Desktop\' ...
%        'ffmpeg-20141210-git-ae81680-win64-static\ffmpeg-20141210-git-ae81680-win64-static\bin']);
%     
%     name = projList{iProj,2};
%       ffmpeg -i pathToFrames/frame_%0Nd.tif -r 15 -b 20000k filename.mp4
%     execute = ['ffmpeg -r 5 -i ' saveDir filesep '%03d.png' ...
%     ' -crf 22 -pix_fmt yuv420p -b 20000k ' saveDir filesep name '.mp4'];
%     system(execute)
%     cd(saveDir)
else
    display(['Image Directory Already Found for ' MD.outputDirectory_ ' Skipping']);
end




