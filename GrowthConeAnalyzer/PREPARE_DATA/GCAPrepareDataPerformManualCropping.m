

function stk2tiffDirsWithSeriesCropMultChannel(varargin)
% stktiffDirs splits STKs in input directory into folders with TIFF files.
% 
% Synopsis:    stk2tiffDirs(path)
%              stk2tiffDirs(path,'Crop','on')
%
% Input:
%      path : optional - path to the directory containing STK files. Can be
%      a string or a cell array of strings. If not input, the user will be
%      asked to select a folder.
%
%     saveDir: optional - path to the directory 
%
%      Optional parameter/value pairs
%           Crop ('on'/'off') :  if 'on', a window will open asking the
%           user to select the region to crop.
%
% Francois Aguet, 09/01/2010
% Maria Bagonis 12/29/2013 : added scrolling option for cropping, added
% linescan tool. 

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('saveDir',[],@(x) ischar(x) || isempty(x) || iscell(x)); 
ip.addParamValue('Crop', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
%ip.addParamValue('Line','off',@(x) strcmpi(x,'on') |strcmpi(x,'off')); 
ip.addParamValue('BioSensors','off', @(x) ischar(x) || strcmpi(x,'off'));
% ip.addParamValue('ChannelOrder', 'interleaved', @(x) any(strcmpi(x, {'interleaved', 'consecutive'})));

% ip.addParamValue('ChannelNames', cell(1));
ip.parse(varargin{:});
stkpath = ip.Results.path;
saveDir = ip.Results.saveDir; 
crop = ip.Results.Crop;
%nc = ip.Results.Channels;
bioOn = ip.Results.BioSensors; 
%cdir = ip.Results.ChannelNames;
% if nc>1
%     if isempty(cdir{1})
%         cdir = arrayfun(@(i) ['c' num2str(i) filesep], 1:nc, 'UniformOutput', false);
%     else
%         cdir = cellfun(@(i) [i filesep], cdir, 'UniformOutput', false);
%     end
% end

if isempty(stkpath)
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end
end

% Recursive call if input is cell array
if iscell(stkpath), 
    cellfun(@(x)stk2tiffDirs(x,'Crop',crop),stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = [dir([stkpath '*.tif']) dir([stkpath '*.tiff']) dir([stkpath '*.stk'])];

N = length(stkList);
if N==0
    fprintf('No TIFF files found in input directory.\n');
end

if isempty(saveDir)
    saveDir = upDirectory(stkpath,1); % make the directory above where the stks are stored
end



% if biosensor option on
% make folders with channels
%if bioOn == 1
tifNames = vertcat({stkList.name})';
movieNamesAll= cellfun(@(x) upDirectory(strrep(x,'_',filesep),1),tifNames,'uniformoutput',0);

movieNames = unique(movieNamesAll);
dataSet.movieNames = movieNames; 

for iMovie = 1:length(movieNames)
    cMovie = movieNames{iMovie};
    
    
    idx = cellfun(@(x) strcmpi(x,cMovie),movieNamesAll);
    
    cTifs = tifNames(idx);
    
    % find the mCherry channel
    [~,cNames] = cellfun(@(x) upDirectory(strrep(x,'_',filesep),1),cTifs,'uniformoutput',0);
    %     idx = strcmpi(cNames,'mCherry');
    %     toCrop = cTifs(idx);
    
    fprintf('Converting: %s\n', cMovie);
   if strcmpi(ip.Results.BioSensors,'on');  

    
    NChannels = 3;
   else 
       NChannels = 1; 
   end 
    % [~,stkname] = fileparts(stkList(k).name);
    
    % get stack info
    %[ny,nx,NFrames] = size(stack);
    stackForCrop = stackRead([stkpath filesep cTifs{end}]); % for now set up to use mCherry which will always be the end channel
    
    [ny,nx,nFrames] = size(stackForCrop);
    dataSet.nFrames = nFrames; 
    % initiate matrix and get other channels in order/ likely faster way to do this
    % by just sorting cTims
    stackAllOrig = zeros(ny,nx,NChannels,nFrames);
    if strcmpi(ip.Results.BioSensors,'on');  
        
    stackAllOrig(:,:,3,:) = stackForCrop;
    else 
        stackAllOrig(:,:,1,:) = stackForCrop; 
    end 
    
    if strcmpi(ip.Results.BioSensors,'on');  
    % get the other two channels
    for iCh = 1:2 
        stackC = stackRead([stkpath filesep cTifs{iCh}]);
        
        stackAllOrig(:,:,iCh,:) = stackC;
        
        
    end
    
    end 
    
    
    
    
    %
    if strcmpi(ip.Results.Crop, 'on')
        
        imseriesshow(stackForCrop);
        set(gcf,'Name','Crop the Image; Double Click Retangle When Finished','NumberTitle','off')
        %  hMsg = msgbox({'Use the expandable rectangle to crop the image'; 'Double Click when finished'},'help');
        
        %stackOrig= stack;
        % h = figure;
        % imagesc(stack(:,:,1)); colormap(gray(256)); axis image;
        
        hr = imrect();
        
        pos = round(wait(hr));
        % cropRoi{k} = pos;
        stackCrop = stackForCrop(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
        %
        close all
        %end
        
        % reshow the new cropped series
        
        reply = 'Yes';
        while strcmpi(reply,'Yes');
            imseriesshow(stackCrop)
            hMsg=  msgbox({'Click through your movie to make sure cropping is correct';'When Finished Click OK'},'help');
            uiwait(hMsg);
            reply = questdlg('Re-do the cropping?');
            if strcmpi(reply,'Yes');
                close gcf
                imseriesshow(stackForCrop);
                hr = imrect();
                pos = round(wait(hr));
                
                stackCrop = stackForCrop(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
                %
                
                
            end % if strcmpi
            
        end
        % finish cropping all channels
        x = pos(1);
        x2 = pos(1)+pos(3);
        y = pos(2);
        y2= pos(2)+pos(4);
        deltX = abs(x2-x+1);
        deltY = abs(y2-y+1);
        stackAllCrop = zeros(deltY,deltX,NChannels,nFrames); % NOTE make channel number fliz
        % for each channel put into the multiChannelStack
        % stackAll = y,x,channel,frame
        % Crop All Channels
        for iChannel = 1:NChannels
            stackAllCrop(:,:,iChannel,:) = stackAllOrig(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),iChannel,:);
        end
       
        
        close gcf
   
        
        
        
        
    end % iCrop
    
    % make directories
    example = strrep(cTifs{1},'_',filesep);
    newDir = [saveDir filesep upDirectory(example,1)];
    if ~isdir(newDir) 
        mkdir(newDir); 
    end 
    if strcmpi(ip.Results.Crop,'on')
    % record cropping coords
        cropRoi{iMovie,1} = pos;
        save([newDir filesep 'cropRegion.mat'],'pos');
    end 
    
    for iChannel = 1:NChannels
        cChanDir = [newDir filesep 'Channels' filesep 'C' num2str(iChannel) '_' cNames{iChannel}];
        % make movieData Object 
        % make channel object 
        channel(iChannel) = Channel(cChanDir);
      %  channel(iChannel).fluorophore_= cNames{iChannel};
        if strcmpi(cNames{iChannel},'mCherry');  
          %  channel(iChannel).emissionWavelength_=name2wavelength(cNames{iChannel})*1e9;
              % Set some channel properties
            % channel.fluorophore_='yfp';
             channel(iChannel).emissionWavelength_=name2wavelength('mCherry')*1e9;
            % channel.imageType_='TIRF';
            %% NOTE need to check if these are the same for Ludo's Experiments 
            NA = 1.4; 
            M = 60; 
            cameraPix= 6.45*1e-6;% in m
%             sigma = getGaussianPSFsigma(NA,M,2*cameraPix,name2wavelength('mCherry')); 
%             channel(iChannel).psfSigma_ = sigma; % used for the fitting. 
         end 
        % channel.imageType_='TIRF';
        
        
        if ~isdir(cChanDir)
            
            mkdir(cChanDir);
        end
        
        stackAllCrop = uint16(stackAllCrop);
        
        for iFrame = 1:nFrames
            imwrite(stackAllCrop(:,:,iChannel,iFrame),[cChanDir filesep 'C' num2str(iChannel) '_' num2str(iFrame,'%03d') '.tif'],'tif');
        end
    end
    imseriesshow(stackForCrop); 
  %% document Characteristic flags : this will be a personal addendum that may 
  % want to add to the neurite profile for future grouping .
  neighborhood = questdlg('Your neurite is...', ...
       'Define Neurite Neighborhood', 'Touching Another Cell?', 'Nearby Another Cell?', 'In Isolation?', 'In Isolation?' ); 
  bifNum = inputdlg('How many bifurcations does your neurite have?');
  failFlag = questdlg('Maria do you think this seg will fail?','My Check', 'yes','no','no') ; 
   neuriteProfile.neighborhood = strrep(neighborhood,'?','');      % 'touch', 'inRegion', 'isolated
   neuriteProfile.bifurcations =  bifNum ;   % 0,1,2... numeric 
   neuriteProfile.failFlag = failFlag; 
%   neuriteProfile.mayFail =    % numeric 0 or 1 - this is mroe a flag for me if I see characteristics that 
  % I know will be handled poorly by the reconstruction 
 
  
  close gcf
    
    % from here set up MovieData and automatically run biosensors and
    % analysis.
    frames = [1 round(nFrames/2) nFrames];
    %% setting up the montage stack
    if iMovie == 1
        mStack = zeros(ny,nx,1,NChannels*N/NChannels);% 
    end
    for iChannel = 1:NChannels
        for iCol = 1:length(frames)
            img = stackAllOrig(:,:,iChannel,frames(iCol));
            
            img = (img-min(img(:)))/(max(img(:))-min(img(:)));
            if strcmpi(ip.Results.BioSensors,'on')
                invert = 1; 
            else 
                invert = -1; 
            end 
            mStack(:,:,iChannel,iMovie*3-3 + iCol) = invert.*img;
        end
    end
    
   % finish making movieData object 
      saveFolder = [ newDir  filesep 'ANALYSIS'];
      if ~isdir(saveFolder) 
          mkdir(saveFolder)           
      end 
      % save neurite neighborhood info
        save([saveFolder filesep 'neuriteProfile.mat'],'neuriteProfile');
      % Constructor needs an array of channels and an output directory (for analysis)
      MD = MovieData(channel,saveFolder);
      % Set the path where to store the MovieData object.
      MD.setPath(saveFolder);
      MD.setFilename('movieData.mat');
          % Run sanityCheck on MovieData. 
        % Check image size and number of frames are consistent. 
        % Save the movie if successfull
        MD.sanityCheck;

        % Set some additional movie properties
      %  MD.numAperture_=1.4;
        MD.pixelSize_=215;
        MD.timeInterval_=5;
        MD.camBitdepth_=16;
       % MD.binning = 2; 
       
        %MD.notes_='Created for test purposes'; 

        % Save the movie
        MD.save;

clear MD
      
      
      
      
     
    
end % for iMovie




close all
%% save Montage Info 
montageDir = [stkpath filesep 'Montages'];
if ~isdir(montageDir) 
mkdir(montageDir)
numDataSetRun =0; 
else 
    load([montageDir filesep 'dataSet.mat']) 
    numDataSetRun = length(dataSet);
    
end 

% maybe cluster the mStacks according to dataset 
dataSet(numDataSetRun+1).mStack = mStack; 
dataSet(numDataSetRun+1).cropRegion = cropRoi; 
save([montageDir filesep 'dataSet.mat'],'dataSet'); 


% make the montage 
GCAPrepareDataMakeCropMontage(dataSet,montageDir); 

%N = N/3; 

% make montages no longer than 3 movies so readable
% NMontages = ceil(N/3);
% test = mod(N,3);
% 
% % loop through making the montages
% for iMontage = 1:NMontages
%     if (iMontage  == NMontages-1 && test >0) ||( NMontages==1 && test>0)
%         final = iMontage*3+test-3;
%     else
%         final  = iMontage*3;
%     end
%     
%     start  = iMontage*3-2;
%     
%    
%     
%    
%     
%     mStackPlot = mStack; 
%     % NOTE just for now to see overlays 
%     % set the third channel to zero 
%    % mStackPlot(:,:,3,:) = 0; 
%     
%     
%     montage(mStackPlot(:,:,:,2*iMontage*3-2:final*3),'DisplayRange',[],'Size',[N,length(frames)]);
%     
%     hold on
%     moviesC = start:final;
%     % plot the cropped regions
%     
%     for iMovie = 1:length(moviesC)
%         pos = cropRoi{moviesC(iMovie),1};
%         x = pos(1);
%         x2 = pos(1)+pos(3);
%         y = pos(2);
%         y2= pos(2)+pos(4);
%         add = [0,nx,nx*2]; % change if make length of frames variable
%         addNy = ny*iMovie-ny;
%         if strcmpi(ip.Results.BioSensors,'on'); 
%             c = 'y'; 
%             cScale = [1,1,0]; 
%         else 
%             c = 'k'; 
%             cScale = [0,0,0]; 
%         end 
%         
%         for iCol = 1:length(frames)
%             % if crop plot the region
%             line([x+add(iCol),x+add(iCol)],[y+addNy +y2+addNy],'color','r','Linewidth',2);
%             line([x+add(iCol),x2+add(iCol)],[y+addNy,y+addNy],'color','r','Linewidth',2);
%             line([x+add(iCol),x2+add(iCol)],[y2+addNy,y2+addNy],'color','r','Linewidth',2);
%             line([x2+add(iCol),x2+add(iCol)],[y+addNy,y2+addNy],'color','r','Linewidth',2);
%             name = [strrep(movieNames{iMovie}, filesep,' '), 'Frame: ' num2str(frames(iCol))];  
%             text(10+add(iCol),20+addNy,name ,'color', c,'FontName','Arial','FontSize',12,'FontWeight','Bold');
%         end
%         
%     end
% pixSizeMic = 0.215; 
%      pixels = round(5/pixSizeMic); %  
%      plotScaleBar(pixels,pixels/5,'Label','5 um','Color',cScale ,'FontSize',5);
%      if ~isdir([saveDir filesep 'Montages'])
%          mkdir([saveDir filesep 'Montages']);
%      end
%     
%     saveas(gcf,[ saveDir filesep 'Montages' filesep  num2str(iMontage) '.fig']); 
%     saveas(gcf,[saveDir filesep 'Montages' filesep num2str(iMontage) '.eps'],'psc2'); 
%   
% %     % save the mStack used to make the montage and the crop region so can
% %     % remake if needed (ie with different overlays etc) 
% %     montageInfo(iMontage).mStack = mStackPlot(:,:,:,iMontage*3-2:final*3); 
% %     montageInfo(iMontage).cropRoi = cropRoi{moviesC}; 
%    
%     
%     close gcf
% end
%  save([saveDir filesep 'Montage' filesep 'montageInfo' num2str(iMontage) '.mat','montageInfo']); 
