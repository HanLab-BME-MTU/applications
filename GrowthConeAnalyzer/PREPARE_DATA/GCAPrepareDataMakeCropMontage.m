function [ output_args ] = makeMontages(dataSet,saveDir)
%
%% save Montage Info 



cropRoi = dataSet.cropRegion; 
mStack = dataSet.mStack; 
% calc the number of total movies
N = numel(cropRoi); 

nFrames = 3;
% make montages no longer than 3 movies so readable
NMontages = ceil(N/3);
test = mod(N,3);
[ny,nx,~,~] = size(mStack); 
% loop through making the montages
for iMontage = 1:NMontages
    % if i montage isn't the last 
    if (iMontage <  NMontages ) ||( iMontage==NMontages && test==0)
       final  = iMontage*3;
    else
         final = iMontage*3+test-3;
        
    end
    
    start  = iMontage*3-2;% movie number start
    
   
    
   
    
    mStackPlot = mStack; 
    % NOTE just for now to see overlays 
    % set the third channel to zero 
   % mStackPlot(:,:,3,:) = 0; 
       frames = [1 round(nFrames/2) nFrames];
    NMovieToPlot = length(start:final); 
    montage(mStackPlot(:,:,:,((iMontage-1)*9)+1:final*3),'DisplayRange',[],'Size',[NMovieToPlot,length(frames)]);
    
    hold on
    moviesC = start:final;
    % plot the cropped regions
    movieNames = dataSet.movieNames(start:final); 
    
    for iMovie = 1:numel(moviesC)
        pos = cropRoi{moviesC(iMovie),1};
        x = pos(1);
        x2 = pos(1)+pos(3);
        y = pos(2);
        y2= pos(2)+pos(4);
        add = [0,nx,nx*2]; % change if make length of frames variable
        addNy = ny*iMovie-ny;
%         if strcmpi(ip.Results.BioSensors,'on'); 
%             c = 'y'; 
%             cScale = [1,1,0]; 
%         else 
            c = 'k'; 
            cScale = [0,0,0]; 
%         end 
        
        for iCol = 1:length(frames)
            % if crop plot the region
            line([x+add(iCol),x+add(iCol)],[y+addNy +y2+addNy],'color','r','Linewidth',2);
            line([x+add(iCol),x2+add(iCol)],[y+addNy,y+addNy],'color','r','Linewidth',2);
            line([x+add(iCol),x2+add(iCol)],[y2+addNy,y2+addNy],'color','r','Linewidth',2);
            line([x2+add(iCol),x2+add(iCol)],[y+addNy,y2+addNy],'color','r','Linewidth',2);
            name = [strrep(movieNames{iMovie}, filesep,' '), 'Frame: ' num2str(frames(iCol))];  
            text(10+add(iCol),20+addNy,name ,'color', c,'FontName','Arial','FontSize',5,'FontWeight','Bold');
        end
        
    end
pixSizeMic = 0.215; 
     pixels = round(5/pixSizeMic); %  
     plotScaleBar(pixels,pixels/5,'Label','5 um','Color',cScale ,'FontSize',5);
     if ~isdir([saveDir filesep 'Montages'])
         mkdir([saveDir filesep 'Montages']);
     end
    
    saveas(gcf,[ saveDir filesep 'Montages' filesep  num2str(iMontage) '.fig']); 
    saveas(gcf,[saveDir filesep 'Montages' filesep num2str(iMontage) '.eps'],'psc2'); 
    saveas(gcf,[saveDir filesep 'Montages' filesep num2str(iMontage) '.tif']); 
%     % save the mStack used to make the montage and the crop region so can
%     % remake if needed (ie with different overlays etc) 
%     montageInfo(iMontage).mStack = mStackPlot(:,:,:,iMontage*3-2:final*3); 
%     montageInfo(iMontage).cropRoi = cropRoi{moviesC}; 
   
    
    close gcf
end

end

