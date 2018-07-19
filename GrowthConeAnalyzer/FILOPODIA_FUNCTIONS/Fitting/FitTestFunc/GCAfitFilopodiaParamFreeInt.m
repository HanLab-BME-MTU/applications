function [ filoInfo] = GCAfitFilopodia( filoInfo,img,varargin)


%% INPUT:
%
%  filoInfo: (REQUIRED) : Rx1 structure arr
%       where R is the number of filopodia like ridges detected
%       .filoInfo has many subfields that store information regarding the
%       filopodias coordinates,orientation,grouping(for branching) etc.
%       notes that each independent ridge is considered a "filopodia"
%       and if branch structures do exist these "filopodia" substructures
%       are marked as related by their fields .groupCount (the scalar
%       identifying the group label, .conIdx (the filoInfo ID of each filopodia
%       ridge branch structrue directly attached to each filopodia ridge)
%       .conXYCoords (the linear indexing of the attachment sites)
%  img: (REQUIRED)
%
%
%  InternalFiloOn: (PARAM) :  scalar
%      1: will fit external only
%      2: will fit internal only
%      3. will fit both
%
% ('NumPixForFitBack' -> scalar) Optional Default = 10 Pixels
%         This parameter dictates the number of pixels back along the filopodia
%         that will be used in the fit relative to the end of the high confidence filopodia tip
%         estimated via the steerable filter response thresholding.
%         Note the signal is often quite noisy along the filopodia - ie there can be
%         multiple possible viable sigmoidal fits. This is especially the
%         case when fitting a lifeAct reporter signal- so one just wants to fit the
%         signal in a local area around the putative tip of the filopodia
%         If the number of pixels in the filopodia is less than this value
%         the entire filopodia length from the thresholded steerable filter response
%         will be used for the respective filopodia tip localization fit.
%
% NOTE TO SELF: this was actually dictated first in the original walkFiloForandBack function
%         ('NumberPixForFitFor' -> scalar) Optional Default = 10 Pixels
%         This parameter dictates the number of pixels forward (relative to
%         the steerable filter threshold set in the filopodia
%         reconstruction) used in the fitting.
%
%
%
% ('ValuesForFit' -> character 'Intensity','Response','Both') Optional
% Default = 'Intensity'
%         The values to use the fitting:
%         Intensity: The values corresponding to the intensity will be fit
%         to a sigmoid
%         Response: The values corresponding to the NMS response from the
%         steerable filter will be fit to a sigmoid.
%         Both: Both of the above operations will be performed- this is
%         mainly for comparison in the early stages of development
%
%('Averaging' -> character ('Gaussian Weighted' , 'Perpendicular'
% For Fitting: Response or
%

%
%('SavePlots' -> logical) Optional Default = 1
%
%%Input check
ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('TSOverlays',true,@(x) islogical(x));

ip.addParameter('InternalFiloOn',3,@(x) isscalar(x));
%ip.addParameter('NumPixForFitBack',10,@(x) isscalar(x));
ip.addParameter('fitLengthInPix',10,@(x) isscalar(x)); 
ip.addParameter('ValuesForFit','Intensity',@(x) ischar(x)); % maybe remove
ip.addParameter('PSFSigma',0.43,@(x) isnumeric(x)) ; %% NOTE CHANGE THIS TO BE READ IN FROM MD.
ip.addParameter('OutputDirectory',pwd,@ischar);
ip.addParameter('filterByBackEst',true); 

ip.parse(varargin{:});
p = ip.Results;
%% Initialize 

toAddCell{1} = 'Ext_';
toAddCell{2} = 'Int_';

% for now int will just be the flag for internal versus external
switch p.InternalFiloOn
    case 1
        typeEnd =1;
        typeStart =1;
        
        if p.TSOverlays
            outDir1 =  [p.OutputDirectory filesep toAddCell{1}];
            if ~isdir(outDir1)
                mkdir(outDir1);
            end
        end
        
    case 2
        typeEnd = 2;
        typeStart = 2;
        if p.TSOverlays
            outDir2 =  [p.OutputDirectory filesep toAddCell{2}];
            if ~isdir(outDir2)
                mkdir(outDir2);
            end
        end
        
        
    case 3
        typeStart  = 1;
        typeEnd = 2;
        if p.TSOverlays
            outDir1 =  [p.OutputDirectory filesep toAddCell{1}];
            if ~isdir(outDir1)
                mkdir(outDir1);
            end
            outDir2 =  [p.OutputDirectory filesep toAddCell{2}];
            if ~isdir(outDir2)
                mkdir(outDir2);
            end
        end
end

%% not sure why I have this here ... think this might be old remanants
% shouldn't need a restart any more I don't think possible delete 20140517
if isfield(filoInfo,'endpointCoordFitPix');
    idx2fill = arrayfun(@(x) find(isempty(x.endpointCoordFitPix)),filoInfo);
else
    idx2fill = 1:numel(filoInfo);
end
numFilo2Fit = length(idx2fill);
%% Filter Image 
% currently use the built in matlab function
H = fspecial('gaussian',3,p.PSFSigma);
imgFilt = imfilter(img,H); % for weighted averaging

%% Start
for iType = typeStart:typeEnd
    
    toAdd = toAddCell{iType};
    
    for ifilo = 1:numFilo2Fit
        
        idxCurrent = idx2fill(ifilo);
        
        % convert to distance along filo in pixels: think about a prettier way to do this.
        xyFilo = filoInfo(idxCurrent).([toAdd 'coordsXY']);
        
        % sometimes particularly in int filo structure have []
        
        if ~isempty(xyFilo);            
            % sometimes if the filopodia is at the border I will pad with NaN - old
            % prompt would skip these.
            toRemove = isnan(xyFilo);
            toRemove = toRemove(:,1);
            xyFilo = xyFilo(~toRemove,:);
        end
        
        % first test to make sure that don't have the walk forward error noted
        % 20141029 (few filopodia have this - will check why to see if we can
        % avoid however the best solution is likely simply to filter it).
        if isempty(intersect(filoInfo(idxCurrent).Ext_pixIndicesFor,filoInfo(idxCurrent).Ext_pixIndicesBack));
         
            if ~isempty(xyFilo); 
                                
                xFilo = xyFilo(:,1);
                yFilo = xyFilo(:,2);
                deltX = diff(xFilo);
                deltY = diff(yFilo);
                
                deltFilo =  arrayfun(@(i) sqrt(deltX(i)^2+deltY(i)^2),1:length(deltX));
                distFilo = zeros(length(xyFilo),1);
                
                % think I can use cumsum here...
                toAddNum = 0;
                for i =2:length(deltFilo)+1
                    distFilo(i) = toAddNum + deltFilo(i-1);
                    toAddNum = distFilo(i);
                end
                
                % save the distance along the filo in a new field %% FIXED
                % 20151123
                filoInfo(idxCurrent).([toAdd 'distFilo']) = distFilo; % the distance (in pixels) along the filopodia- will be helpful for calculating branch distances could test if this is a field. .
                
                
                if strcmpi(p.ValuesForFit,'Intensity') % might remove this option entirely
                    maskIndices = filoInfo(idxCurrent).([toAdd 'maskIndices']);
                    % quick fix 20141023
                    
                    if length(distFilo) ~= length(maskIndices(:,1))
                        display(['Maria you need to make sure to truncate internal maskIndices too that span body for Filo' num2str(ifilo) '!'])
                        maskIndices = maskIndices(1:length(distFilo),:);   
                    end
                                        
                    % First column is of the mask indices is the primary
                    % values, rows from veil to tip 
                    yData = imgFilt(maskIndices(:,1)); % 
                    % save the yData for later testing intensities before
                    filoInfo(idxCurrent).([toAdd 'weightedAvg']) = yData;
                    
                    % Reformat mask indices for plotting. 
                    % I saved these originally with the first column the center line and
                    % the perpendicular pixels outward... might want to save in this
                    % form to avoid future confusion (might have decided not to change
                    % just because of other input that would make it a pain...check 20140517)
                    maskIndices = maskIndices(:,[5,3,1,2,4]);
                    maskIndices = maskIndices'; % switch for plotting
                    
                else % to remove
                    
                    yData = filoInfo(idxCurrent).([toAdd valuesForFit]);
                    
                end
                
                yData = yData(~isnan(yData)); % sometimes I had to pad with NaNs
                
         if  (~isempty(yData) && length(yData) > 4)
                
                %% Quick and Dirty : Estimate the number of sigmoidals
                % get the values where the sigmoidal is maximally
                % increasing or decreasing
                
                
                [slopeMaxNeg,slopeMaxPos,valuesNeg,valuesPos] = gcaFindPotentialSigmoidals(yData,'makePlot',false,'outPath',[p.OutputDirectory filesep toAdd],...
                    'forTitle', ['Filopodia_' num2str(idxCurrent,'%03d') toAdd '_der']);
                
                %% Quick and Dirty : Estimate the mean intensity background
                forSearch =filoInfo(idxCurrent).([toAdd 'pixIndicesFor']);
                nSearch = length(forSearch(~isnan(forSearch))); % make sure don't include nans (sometimes had to truncate the 
                % search forward if the filopodia are located at the
                % border. 
                % Check the size of the mask indices 
               

                % get the pixel indices from the forward search mask 
                % (minus the central position which is still potentially 
                % filopodia signal 
                indicesBack = maskIndices([1:2,4:5],end-nSearch:end); 
                
                % take out nan
                indicesBack = indicesBack(~isnan(indicesBack));
                
                % take the average background signal intensity from the
                % filtered image
                backEst = imgFilt(indicesBack(:)) ; %
                backEstMean = nanmean(backEst); 
                
               
                % sometimes there isn't a well defined local negative slope
                % use the end point from the thresholded steerable filter
                % NMS response
                if isempty(slopeMaxNeg)
                    slopeMaxNeg = length(filoInfo(idxCurrent).([toAdd 'pixIndicesBack']));
                end
               
                            
                %% Additional External Quick and Dirty: (Noise Filter) find the first slopeMaxNeg before descent into the mean background estimation
                
                % Potentially treat embedded versus external filopodia differently
                % embedded search for signal closer to the veil and cut off
                % after first significant increase in the signal.
                switch toAdd
                    case 'Ext_'
                        if ip.Results.filterByBackEst
                           % slopeMaxNeg(yData(slopeMaxNeg)<backEstMean) = [];
                            
                             noiseDip = find(yData<backEstMean,1,'first');
                             if ~isempty(noiseDip)
                                slopeMaxNeg = slopeMaxNeg(slopeMaxNeg < noiseDip); 
                             end 
                        end
                        
                       
                        
                        if isempty(slopeMaxNeg)
                            slopeMaxNeg = length(filoInfo(idxCurrent).([toAdd 'pixIndicesBack']));
                        end
                        
                      
                                               
                        %%  perform fitting around the first sigmoidal closest to the background
                        
                        if length(yData)-slopeMaxNeg(end) <ip.Results.fitLengthInPix ;
                            endFit = length(yData);
                        else
                            endFit  = slopeMaxNeg(end)+ip.Results.fitLengthInPix;
                        end
                        
                        if slopeMaxNeg(end)<=ip.Results.fitLengthInPix;
                            startFit = 1;
                        else
                            startFit = slopeMaxNeg(end)-ip.Results.fitLengthInPix;
                        end
                        
                        if length(slopeMaxNeg)>1
                            if slopeMaxNeg(end-1) >= startFit;
                                startFit= slopeMaxNeg(end-1) +2; % start after the slope fit
                            end
                        end
                        
                        
                        
                        % but check to make sure not including increasing portions
                        % (should check if increasing portions signficiant if want to be rigourous..);
                        
                        if ~isempty(slopeMaxPos) % you have portions of the linescan with positive slope
                            
                            % make sure that do not have an increasing portion
                            % between tentative sigmoidal mean and end point of fit
%                             if (sum(endFit>slopeMaxPos)>0 && sum(slopeMaxPos > slopeMaxNeg(end))>0);
%                                 idx =  find(endFit>slopeMaxPos);
%                                 endFit = slopeMaxPos(idx(1)) -2;
%                             end

                            tentMean = slopeMaxNeg(end); 
                            % get the indices between 
                            between1 = tentMean:endFit; 
                            between2 = startFit:tentMean; 
                            test1 = intersect(between1,slopeMaxPos); 
                            test2 = intersect(between2,slopeMaxPos); 
                          
                            if ~isempty(test1)
                                endFit = between1(between1==test1(1))-2; 
                            end 
                            
                            if ~isempty(test2)
                                startFit = between2(between2==test2(1))+2; 
                            end 
                            
                            % make sure that do not have an increasing portion
                            % start of fit tentative sigmoidal mean
%                             if (sum(startFit<slopeMaxPos)>0 && sum(slopeMaxPos < slopeMaxNeg(end))>0);
%                                 idx = find(startFit<slopeMaxPos);
%                                 startFit = slopeMaxPos(idx(1))+2;mp
%                             end
                              


                        end
                        
                     
                        
                        
                    case 'Int_'
                        
                        % first test if there is any increasing portions of the curve if so it is a flag you are in noise.
                        
                        if ~isempty(slopeMaxPos)
                            testEnd = slopeMaxPos(slopeMaxPos ==min(slopeMaxPos));
                        else
                            testEnd = length(yData);
                        end
                        
                        % find the last negative slope before increase
                        slopeMaxNeg(slopeMaxNeg>testEnd) = [];
                        if isempty(slopeMaxNeg)
                            slopeMaxNeg = length(filoInfo(idxCurrent).([toAdd 'pixIndicesBack']));
                        end
                        
                        if length(yData)-slopeMaxNeg(end) <ip.Results.fitLengthInPix ;
                            endFit = length(yData);
                        else
                            endFit  = slopeMaxNeg(end)+ip.Results.fitLengthInPix;
                        end
                        
                        if slopeMaxNeg(end)<=ip.Results.fitLengthInPix;
                            startFit = 1;
                        else
                            startFit = slopeMaxNeg(end)-ip.Results.fitLengthInPix;
                        end
                        
                        if length(slopeMaxNeg)>1
                            if slopeMaxNeg(end-1) >= startFit;
                                startFit= slopeMaxNeg(end-1) +2; % start after the slope fit
                            end
                        end
                        
                        if ~isempty(slopeMaxPos) % you have portions of the linescan with positive slope
                            
                             tentMean = slopeMaxNeg(end); 
                            % get the indices between 
                            between1 = tentMean:endFit; 
                            between2 = startFit:tentMean; 
                            test1 = intersect(between1,slopeMaxPos); 
                            test2 = intersect(between2,slopeMaxPos); 
                          
                            if ~isempty(test1)
                                endFit = between1(between1==test1(1))-2; 
                            end 
                            
                            if ~isempty(test2)
                                startFit = between2(between2==test2(1))+2; 
                            end 
                            
                            
                            
                            
                            % make sure that do not have an increasing portion
                            % between tentative sigmoidal mean and end point of fit
%                             if (sum(endFit>slopeMaxPos)>0 && sum(slopeMaxPos > slopeMaxNeg(end))>0);
%                                 idx =  find(endFit>slopeMaxPos);
%                                 endFit = slopeMaxPos(idx(1)) -2;
%                             end
%                             
%                             % make sure that do not have an increasing portion
%                             % start of fit tentative sigmoidal mean
%                             if (sum(startFit<slopeMaxPos)>0 && sum(slopeMaxPos < slopeMaxNeg(end))>0);
%                                 idx = find(startFit<slopeMaxPos);
%                                 startFit = slopeMaxPos(idx(1))+2;
%                             end
                        end
                        if endFit < startFit 
                            display(['Check filo' num2str(ifilo)]); 
                        end 
                        
                        
                end % toAdd
                
                % Define the data for the fitting.
                yDataFit = yData(startFit:endFit);
                distFiloFit = distFilo(startFit:endFit);
                
                % final check 
%                 if isempty(yDataFit)
%                     
%                 end 
         else 
             yDataFit = []; 
             distFiloFit = []; 
         end 
%% Perform the Fit      

                if  (~isempty(yDataFit) && length(yDataFit) > 4) % make a cut off that length yData needs to be at least four pixels.           
                      
                    % CHECK these starts (20140517)
                    starts = [max(yDataFit)-min(yDataFit),length( filoInfo(idxCurrent).([toAdd 'pixIndicesBack'])),5,...
                        min(yDataFit)];
                    upb = [max(yData); max(distFiloFit);20 ; max(yData)];
                    lwb = [0;min(distFiloFit); 1; min(yData)];
     
                    [params, resnorm,resid,exitFlag] = lsqcurvefit(@testfun,starts,distFiloFit,yDataFit,lwb,upb);
                    
                    % Add Parameters from fit to structure 
                    filoInfo(idxCurrent).([toAdd 'params'])  = params;
                    
                    %find an approximate single value index of the filo tip position 
                    % (as measured along filo: 0  = the base of the filo)
                    % so can get a quick xy coord for plotting.
                    delt = abs(distFilo-params(2)); % Whoops 20160128 need this to be distFilo now (notDistFiloFit) (pixIndices are not truncated)
                 
                    pixIdx = find(delt==min(delt));
                    
                    pixIndices = filoInfo(idxCurrent).([toAdd 'pixIndices']);
                    if pixIdx == 0 % indicates the fit was bad set fit parameters to NaN
                       
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) = [NaN NaN];
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
                        filoInfo(idxCurrent).([toAdd 'length']) = NaN;
                        filoInfo(idxCurrent).([toAdd 'std']) = NaN;     
                    else
                        
                        % Record the image coordinates of the fitted filopodia endpoint
                        % NOTE:  should someday try to translate 
                        % the fitted coordinate system to get the subpixel
                        % coordinates in the image for plotting.
                        % However, at this point we just need the accurate
                        % sub-pixel length metrics, and the eye can't really 
                        % tell much of a difference between a pixel vs
                        % subpixel coordinate. 
                        
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = pixIndices(pixIdx);% for graphing - 
                        [ny,nx] = size(img);
                        [y,x] = ind2sub([ny,nx],filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']));
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) = [x,y]; % for graphing an approximate xy coord.
                                        
                        % Think about removing 'length' and 'std' field as a bit redundant
                        % with above - however keep for now as most of the 
                        % downstream functions use this field. 
                        filoInfo(idxCurrent).([toAdd 'length']) = params(2); % 
                        
                        filoInfo(idxCurrent).([toAdd 'std']) = params(3); % std includes that of the psf of the microscope.
                        
                    end
                    
                    filoInfo(idxCurrent).([toAdd 'exitFlag']) = exitFlag; % always include exit flag
                    filoInfo(idxCurrent).([toAdd 'resnorm']) = resnorm;
                    %filoInfo(idxCurrent).([toAdd 'resnorm']) = resnorm/length(distFiloFit); % divide by the number of points
                    filoInfo(idxCurrent).([toAdd 'resid']) = resid; 
                   
               
                    
                    yFit = params(1)*(1-0.5*(1+erf((distFilo-params(2))/(params(3)*sqrt(2)))))+params(4);
                    
                   %% Optional Overlays
                   
                    if p.TSOverlays 
                        
                        % Initiate
                        fsFigure(0.75,'visible','off'); 
                        
                        %% First Subplot: View Local Detection of Ridge Signal 
                        if sum(isnan(maskIndices(:)))==0; % if no-mask NaNs (ie not at the border); 
                            subplot(3,2,(1:2));
                            imgFilo = img(maskIndices);
                            imagesc(-imgFilo,[-max(img(:)),-min(img(:))]);
                            
                            % mark pixels for potential background
                            % estimation 
                            
                            backMask = zeros(size(imgFilo)); 
                            backMask([1:2,4:5],end-nSearch:end) = 1;
                            idx = find(backMask); 
                            [y,x] = ind2sub(size(backMask),idx); 
                            hold on 
                            scatter(x,y,100, [0.5020, 0.4510, 0.6745],'x'); % color in purple for now 

                            ylabel('Filopodia Width (Pixels)'); 
                            axis([0.5,size(maskIndices,2),0.5,5.5])
                            title(['Filopodia Detection ' num2str(idxCurrent) ' : Purple Crosses Mark Local Background'],'FontSize',10,'FontName','Arial');  
                        end
                        
                        %% Second Subplot: Intensity Values Along Automated Linescan Etc
                        subplot(3,2,(3:4)); 
                        hData = scatter(distFilo,yData,50,'k');
                        hold on 
                        
                        % mark the estimated local background value
                        hBack = line([0 distFilo(end)], [backEstMean,backEstMean],'color','k','Linestyle','--');
                        
                        % mark the data used for the fit 
                        hDataFit =  scatter(distFiloFit,yDataFit,'k','filled'); % color in the data specifically used for the fitting.
                        
                        % plot the fit
                        hFit =  plot(distFilo,yFit,'r'); 
                        
                        hMean = line([params(2),params(2)],[min(yData),max(yData)],'Color','k','Linewidth',2); 
                        
                        % plot sites of potential decay and positive slopes
                        % for troubleshooting
                        hDec = scatter(distFilo(slopeMaxNeg),yData(slopeMaxNeg),50,[ 0.0039  ,  0.4264 ,   0.3848],'filled'); % cyan 
                        hInc = scatter(distFilo(slopeMaxPos),yData(slopeMaxPos),50,[ 0.6471 ,        0 ,   0.1490  ] ,'filled');
                        
                         legend([hData,hDataFit,hFit,hDec,hBack,hMean,hInc],'Raw Data', 'Data For Fit','Fit','Potential Signal Decay',... 
                             'Background Estimate','Tip Position','Potential Signal Rise', 'Location','Northeast','FontSize',6,'FontName','Arial'); 
                         legend('boxoff');
                        xlabel('Distance Along Filopodia (Pixels)','FontSize',10,'FontName','Arial');
                        ylabel('Fluorescence Intensity (AU)','FontSize',10,'FontName','Arial'); 
                    if exitFlag>1 
                        title({['Length (Mean Gaussian Survival) = ' num2str(params(2).*0.216,3) ' um'] , ['Intensity Filament (Amplitude) = ' num2str(params(1,1),3)] })                       
                    end
                        
                        %% Third Subplot: Filopodia Overlay  
                        subplot(3,2,5);
                        imshow(-img,[]) ;
                        hold on
                        filoInfoC = filoInfo(idxCurrent);
                        text(5,5,['Filopodia ' num2str(idxCurrent)],'FontSize',6,'FontName','Arial');
 
                        pixelsF = zeros(size(img));
                        pixF= filoInfoC.([toAdd 'pixIndicesFor']);
                        pixF = pixF(~isnan(pixF));
                        pixelsF(pixF)=1;
                        % pixelsB = zeros(size(img));
                        % pixelsB(filoInfoC.([toAdd 'pixIndicesBack']))=1;
                        [yb,xb] = ind2sub(size(img), filoInfoC.([toAdd 'pixIndicesBack'])); 
                        scatter(xb,yb,10,[ 0.0039  ,  0.4264 ,   0.3848],'filled'); 
%                         spy(pixelsF,'g',10);
%                         spy(pixelsB,'r');
                        [yf,xf] = ind2sub(size(img),filoInfoC.([toAdd 'pixIndicesFor'])); 
                        scatter(xf,yf,10,[0.5020, 0.4510, 0.6745],'filled');  % purple                     
                        %% 4th Subplot the residuals compared to the amplitude
                        subplot(3,2,6);
                        boxplot(resid);
                        hold on
                        thresh = prctile(resid,95); 
                        if params(1) > thresh; 
                            add = '*'; 
                        else 
                            add = []; 
                        end 
                        yLimHi = max(img(:));
                        yLimLo = min(resid)-2;
                        axis([0.5 1.5 yLimLo yLimHi]);
                        
                        scatter(1,params(1),'r','filled');
                        text(1,params(1),['Intensity Filament ' add],'FontSize',6,'FontName','Arial'); 
                        ylabel({'Residuals' ; 'From Local Fit'},'FontSize',6,'FontName','Arial'); 
                        
                                            
                        filename{1} = ['Filopodia_' num2str(idxCurrent,'%03d') toAdd '.fig' ];
                        filename{2} = ['Filopodia_' num2str(idxCurrent,'%03d') toAdd '.png'];
                        
                        for i = 1:2
                            saveas(gcf, [ p.OutputDirectory filesep toAdd filesep filename{i}]);
                        end
                        saveas(gcf,[p.OutputDirectory filesep toAdd filesep 'Filopodia_' num2str(idxCurrent,'%03d') toAdd '.eps'],'psc2'); 
                    end
                else
                    filoInfo(idxCurrent).([toAdd 'params']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) =[ NaN,NaN];
                    filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'exitFlag']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'resnorm']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'resid']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'length']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'std']) = NaN;
                    % new field added 05/17/2014
                    filoInfo(idxCurrent).([toAdd 'distFilo']) = NaN; %
                                   
                end
            else
                filoInfo(idxCurrent).([toAdd 'params']) = NaN;
                filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) =[ NaN,NaN];
                filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
                filoInfo(idxCurrent).([toAdd 'exitFlag']) = NaN;
                filoInfo(idxCurrent).([toAdd 'resnorm']) = NaN;
                filoInfo(idxCurrent).([toAdd 'resid']) = NaN; 
                filoInfo(idxCurrent).([toAdd 'length']) = NaN;
                filoInfo(idxCurrent).([toAdd 'std']) = NaN;
                % new field added 05/17/2014
                filoInfo(idxCurrent).([toAdd 'distFilo']) = NaN; %
                
                
                
                
            end % if isnan
        else
            filoInfo(idxCurrent).([toAdd 'params']) = NaN;
            filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) =[ NaN,NaN];
            filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
            filoInfo(idxCurrent).([toAdd 'exitFlag']) = -20; % make a new error exit flag ...
            filoInfo(idxCurrent).([toAdd 'resnorm']) = NaN;
            filoInfo(idxCurrent).([toAdd 'length']) = NaN;
            filoInfo(idxCurrent).([toAdd 'std']) = NaN;
            % new field added 05/17/2014
            filoInfo(idxCurrent).([toAdd 'distFilo']) = NaN; %
        end % is empty
        close all
    end % ifilo 
end % end iType

end
function F= testfun(p,xData)
%F = p(1).*xData + p(2);
F = p(1)*(1-0.5*(1+erf((xData-p(2))/(p(3)*sqrt(2)))))+p(4);
end



