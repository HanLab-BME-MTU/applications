function plotSelectedWindowsOverlay(MD,winCluster,iChan,frameRange,movieOn)
%      plotSelectedWindowsOverlay(MD,winCluster,iChan,frameRange,movieOn)
  
iProt = MD.getProcessIndex('ProtrusionProcess');

if isempty(iProt)
    error('This movie does not have a valid protrusion process! Please run protrusion calculation first!')
end


nFrames = MD.nFrames_;
nChan = numel(MD.channels_);

if nChan == 1
    iChan = 1;
end


if isempty(frameRange)
    frameRange = 1:nFrames;
end

prot    = MD.processes_{iProt}.loadChannelOutput;
winIdx  = MD.getProcessIndex('WindowingProcess');

%% ----- Figure Making ------ %%

ny = MD.imSize_(1);
nx = MD.imSize_(2);


% Configure figure

h = figure('Visible', 'on', 'Position', [50 50 nx ny]);

iptsetpref('ImshowBorder','tight');

set(h, 'InvertHardcopy', 'off');

set(h, 'PaperUnits', 'Points');

set(h, 'PaperSize', [nx ny]);

set(h, 'PaperPosition', [0 0 nx ny]); % very important

set(h,'DefaultLineLineSmoothing','on');

set(h,'DefaultPatchLineSmoothing','on');

axes('Position',[0 0 1 1]) 
hold on

nCluster  = numel(winCluster);
colorCode = {'b','r'} ;

if movieOn
    
    
    framesPath = [MD.outputDirectory_ filesep 'movieOverlay'];
    if ~isdir(framesPath)
        mkdir(framesPath)
    end

    
    fmt = ['%0' num2str(ceil(log10(numel(frameRange)))) 'd'];
    ext = '.png';
    zoom = 1;
    
end


for iFrame = frameRange

   iWindow = MD.processes_{winIdx}.loadChannelOutput(iFrame,iChan) ;
    K = imadjust(MD.channels_(iChan).loadImage(iFrame));
   %imagesc(1:nx,1:ny,MD.channels_(iChan).loadImage(iFrame))
   imagesc(1:nx,1:ny,K)
   %imshow(MD.channels_(iChan).loadImage(iFrame))
   %imagesc(1:nx,1:ny,MD.processes_{8}.loadChannelOutput(2,iFrame))
   axis([1 nx 1 ny])
   colormap(1-gray)

   hold on
   cellfun(@(x,y) plotWindows(iWindow(x),{y,'FaceAlpha',1},0,'bandMin',1,'bandMax',1),winCluster,colorCode,'Unif',0)
   hold off
   
   if movieOn
       print(gcf, '-dpng', '-loose', ['-r' num2str(zoom*72)], [framesPath filesep 'frame' num2str(iFrame, fmt) ext]);
       cla(gca)
   else
       pause(0.01)
   end
   
end