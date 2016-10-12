function [ filoInfo] = GCAfitFilopodia( filoInfo,img,varargin)


%% INPUT:
%
%  filoInfo: (REQUIRED) : Rx1 structure array
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
ip.addParameter('NumPixForFitBack',10,@(x) isscalar(x));
ip.addParameter('ValuesForFit','Intensity',@(x) ischar(x)); % maybe remove 
ip.addParameter('PSFSigma',0.43,@(x) isnumeric(x)) ; %% NOTE CHANGE THIS TO BE READ IN FROM MD. 
ip.addParameter('OutputDirectory',pwd,@ischar); 

ip.parse(varargin{:});
p = ip.Results;
%%

toAddCell{1} = 'Ext_';
toAddCell{2} = 'Int_';

% for now int will just be the flag for internal versus external
switch p.InternalFiloOn
    case 1
        typeEnd =1;
        typeStart =1;
    case 2
        typeEnd = 2;
        typeStart = 2;
    case 3
        typeStart  = 1;
        typeEnd = 2;
end




%% not sure why I have this here ... think this might be old remanants
% shouldn't need a restart any more I don't think possible delete 20140517
if isfield(filoInfo,'endpointCoordFitPix');
    idx2fill = arrayfun(@(x) find(isempty(x.endpointCoordFitPix)),filoInfo);
else
    idx2fill = 1:numel(filoInfo);
end

numFilo2Fit = length(idx2fill);
%%
% currently use the built in matlab funtction but potentially use other
H = fspecial('gaussian',3,p.PSFSigma);
imgFilt = imfilter(img,H); % for weighted averaging


for iType = typeStart:typeEnd
    toAdd = toAddCell{iType};
    for ifilo = 1:numFilo2Fit
        
        
        idxCurrent = idx2fill(ifilo);
        
        numBack = length(filoInfo(idxCurrent).([toAdd 'pixIndicesBack']));
        numFor = length(filoInfo(idxCurrent).([toAdd 'pixIndicesFor']));
        
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
            
            
            if ~isempty(xyFilo); % NOTE needed to change as sometimes have NaN padded on here
                % an example of this is Filo 44 FRAME 08 CDC42KD
                
                
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
                
                % save the distance along the filo in a new field
                filoInfo(idxCurrent).distFilo = distFilo; % the distance (in pixels) along the filopodia- will be helpful for calculating branch distances could test if this is a field. .
                
                
                if strcmpi(p.ValuesForFit,'Intensity')
                    maskIndices = filoInfo(idxCurrent).([toAdd 'maskIndices']);
                    % quick fix 20141023
                    
                    if length(distFilo) ~= length(maskIndices(:,1))
                        display(['Maria you need to make sure to truncate internal maskIndices too that span body for Filo' num2str(ifilo) '!'])
                        maskIndices = maskIndices(1:length(distFilo),:);
                        
                    end
                    yData = imgFilt(maskIndices(:,1)); % just take the center line check to see the order of the data row 1 = base I think wherease end = tip
                    if p.TSOverlays  == true
                        if sum(isnan(maskIndices(:)))==0;
                            figure('visible','off')
                            subplot(3,2,(1:2));
                            % I saved these originally with the first column the center line and
                            % the perpendicular pixels outward... might want to save in this
                            % form to avoid future confusion (might have decided not to change
                            % just because of other input that would make it a pain...check 20140517)
                            maskIndices = maskIndices(:,[5,3,1,2,4]);
                            % maskIndices = maskIndices(:,[5,2,3,1,4]);
                            maskIndices = maskIndices';
                            imgFilo = img(maskIndices);
                            imagesc(imgFilo,[min(img(:)),max(img(:))]);
                        end
                    end
                    
                    
                    
                else
                    
                    yData = filoInfo(idxCurrent).([toAdd valuesForFit]);
                    
                end
                
                yData = yData(~isnan(yData)); % sometimes I had to pad with NaNs
                numForFit = p.NumPixForFitBack;
                % plot(distFilo,yData);
                if numBack > p.NumPixForFitBack
                    %xData = 1:length(filoInfo(idxCurrent).pixIndices(
                    
                    
                    yDataFit = yData(end-(numFor+numForFit)+1:end); % truncate to just fit the end pixels and the forward projections
                    distFiloFit = distFilo(end-(numFor+numForFit)+1:end);
                    % possible change
                    measEnd = numForFit+1;
                    pixIndices = filoInfo(idxCurrent).([toAdd 'pixIndices'])(end-numFor-numForFit+1:end);
                    %
                else
                    yDataFit = yData; %use the whole thing
                    % possible change
                    measEnd = numBack+1;
                    pixIndices = filoInfo(idxCurrent).([toAdd 'pixIndices']);
                    %
                    distFiloFit = distFilo;
                    
                end % if numBack
                
                % I might want to change how I save this..
                filoInfo(idxCurrent).([toAdd 'pixIndicesUsedForFit']) = pixIndices;
                
                
                
                %          x1 = -numBack:-1;
                %          x2 = 0:numFor-1;
                %          xData = [x1,x2];
                % yData = filoInfo(idxCurrent).response;
                
                %yData = yData(~isnan(yData));
                % measEnd = numBack+1;
                %xData  = 1:length(yData); % ok we'll see if this fucks up
                
                
                
                % plot(distFilo,yData)
                
                
                % yData = yData/max(yData); % normalize;
                % if (sum(~isnan(filoInfo(idxCurrent).([toAdd p.ValuesForFit]))) == sum(~isnan(filoInfo(idxCurrent).([toAdd 'pixIndices']))) ...
                if      (~isempty(yData) && length(yData) > 4) % make a cut off that length yData needs to be at least four pixels..
                    % I am getting some errors in the fitting function (and they aren't
                    % relavent with the very short filo) Note added 03-02-2013
                    % i am not sure why they wouldn't be but this would be bad (i'll add it given my previous check below 02-22
                    %     if (length(xData) == length(yData) && numBack>3) % hmmm not sure why I made this criteria that the data had to be equal?? maybe sometimes the response was nan? now should always be equal
                    if p.TSOverlays ==1
                        
                        %     figure('visible','on');
                        subplot(3,2,3:4)
                        
                        linescan = scatter(distFilo,yData);
                        xlabel('Distance Along Filopodia Pixels');
                        hold on
                    end
                    % line([measEnd,measEnd],[max(yData),min(yData)],'color','g');%
                    
                    % calculate starts
                    % 4 params corresponding to
                    % p1 = intensity without background
                    % p2 = mean filo length
                    % p3 = combined std of filo and PSF
                    % p4 = background intensity
                    
                    % OLD STARTS: PRE 2014_01_31
                    % CHECK these starts (20140517)
                    starts = [max(yDataFit)-min(yDataFit),length( filoInfo(idxCurrent).([toAdd 'pixIndicesBack'])),5,...
                        min(yDataFit)];
                    upb = [max(yData); max(distFiloFit);20 ; max(yData)];
                    lwb = [0;min(distFiloFit); 1; min(yData)];
                    
                    % calculate lb and ub
                    
                    
                    [params, resnorm,~,exitFlag] = lsqcurvefit(@testfun,starts,distFiloFit,yDataFit,lwb,upb);
            %%        
                    % test if u = min(distFiloFit)
                    if (abs(params(2)-distFiloFit(1))<1 && params(2)> 1)
                        distFiloFitOld = distFiloFit;
                        yDataFitOld = yDataFit;
                        distFiloFit = distFilo;
                        
                        yDataFit = yData;
                        
                        starts = [max(yDataFit)-min(yDataFit),length( filoInfo(idxCurrent).([toAdd 'pixIndicesBack'])),5,...
                            min(yDataFit)];
                        
                        upb = [max(yData); max(distFiloFit);20 ; max(yData)];
                        lwb = [0;min(distFiloFit); 1; min(yData)];
                        
                        
                        % try again using whole values
                        [params,resnorm,~,exitFlag] = lsqcurvefit(@testfun,starts,distFiloFit,yDataFit,lwb,upb);
                        % reset the pix indices to use the whole thing
                        pixIndices = filoInfo(idxCurrent).([toAdd 'pixIndices']);
                        filoInfo(idxCurrent).([toAdd 'pixIndicesUsedForFit']) = pixIndices;
                        
                        
                    end % abs params(2)
             %%       
                    if p.TSOverlays
                        line([params(2),params(2)],[max(yData),min(yData)]);
                        %text(params(2),(max(yDataFit)-min(yDataFit)./2)+min(yDataFit),['mean = ' num2str(params(2),3)]);
                    end
                    
                    filoInfo(idxCurrent).([toAdd 'params'])  = params;
                    
                    % need to change now that do the distance calc...
                    %  pixIdx = round(params(2)); % here rounding the pix indices not very precise calculation should fix
                    % guess this is still decent other option is to have an idx
                    delt = abs(distFiloFit-params(2)); % find the pixIdx (as measured along filo: 0  = the base of the filo so can get a quick xy coord for plotting.
                    % changed above from distFilo to distFiloFit = 20141017 think this
                    % is why the indexing was off
                    
                    pixIdx = find(delt==min(delt));
                    
                    
                    %   pixIndices = filoInfo(idxCurrent).pixIndices;
                    if pixIdx == 0 % indicates the fit was bad set fit parameters to NaN
                        %filoInfo(idxCurrent).endpointCoordFit = NaN;
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) = [NaN NaN];
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
                        filoInfo(idxCurrent).([toAdd 'length']) = NaN;
                        filoInfo(idxCurrent).([toAdd 'std']) = NaN;
                        
                        
                    else
                        
                        % pixIndices = flip(pixIndices); % NOTE 20141017 try flipping pixIndices
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = pixIndices(pixIdx);% for graphing
                        [ny,nx] = size(img);
                        [y,x] = ind2sub([ny,nx],filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']));
                        filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) = [x,y]; % for graphing an approximate xy coord.
                        
                        % actually don't think i need this now that I have converted to
                        % distances the mean should still represent the shifted
                        % sigmoidal.
                        filoInfo(idxCurrent).([toAdd 'length']) = params(2); % maybe fix this later
                        
                        
                        %             if numBack > numForFit
                        %                 numPixNotInFit =  numBack - numForFit;
                        %             filoInfo(idxCurrent).([toAdd 'length']) = params(2) + numPixNotInFit; % use the fitted param value for the length calc
                        %             else
                        %                 filoInfo(idxCurrent).([toAdd 'length']) = params(2);  % use whole thing
                        %             end % if numBack
                        
                        
                        
                        % filoInfo(idxCurrent).length = params(2) ; % ultimately would like a subpixel accuracy here... however would have to do some fitting
                        % i think as these can be curved therefore not a direct
                        % conversion between 'pixels' and xy coordinate absolute
                        % lengths;
                        filoInfo(idxCurrent).([toAdd 'std']) = params(3); % std includes that of the psf.
                        
                    end
                    
                    filoInfo(idxCurrent).([toAdd 'exitFlag']) = exitFlag; % always include exit flag
                    
                    filoInfo(idxCurrent).([toAdd 'resnorm']) = resnorm/length(distFiloFit); % divide by the number of points
                    
                    if exitFlag>1 && p.TSOverlays == 1
                        title({['Length = ' num2str(params(2).*0.216,3) ' um'] , ['p1 = ' num2str(params(1,1))] , ['ExitFlag = ' num2str(exitFlag)]})
                        
                    end
                    %yFit = params(1)*distFilo+params(2);
                    
                    yFit = params(1)*(1-0.5*(1+erf((distFilo-params(2))/(params(3)*sqrt(2)))))+params(4);
                    
                    
                    hold on
                    if p.TSOverlays ==1
                        plot(distFilo,yFit,'r')
                        scatter(distFiloFit,yDataFit,'b','filled'); % color in the data specifically used for the fitting.
                        if exist('distFiloFitOld','var')
                            scatter(distFiloFitOld,yDataFitOld,'g');
                            clear distFiloFitOld yDataFitOld
                            
                        end
                        %             title({['Filo' num2str(idxCurrent)];'Green Line: Endpoint Skel'; 'Blue Line: Endpoint Fit'});
                        filename{1} = ['Filopodia_' num2str(idxCurrent,'%03d') toAdd '.fig' ];
                        filename{2} = ['Filopodia_' num2str(idxCurrent,'%03d') toAdd '.png'];
                        
                        %             if (exitFlag >= 1 && ~isnan(filoInfo(idxCurrent).endpointCoordFitXY(1,1)))
                        %                 currentSaveDir = goodFitSaveDir;
                        %                 maskGoodFit(filoInfo(idxCurrent).([toAdd 'pixIndices'])) =1; % add to good fit mask
                        %
                        %             else
                        %                 currentSaveDir = badFitSaveDir;
                        %                 maskBadFit(filoInfo(idxCurrent).([toAdd 'pixIndices']))=1;
                        %             end
                        
                        subplot(3,2,5:6);
                        imshow(-img,[]) ;
                        hold on
                        filoInfoC = filoInfo(idxCurrent);
                        %GCAVisualsMakeOverlaysFilopodia(filoInfoC,[ny,nx],1,iType,[],0); Ex
                        pixelsF = zeros(size(img));
                        pixF= filoInfoC.([toAdd 'pixIndicesFor']);
                        pixF = pixF(~isnan(pixF));
                        pixelsF(pixF)=1;
                        pixelsB = zeros(size(img));
                        pixelsB(filoInfoC.([toAdd 'pixIndicesBack']))=1;
                        spy(pixelsF,'g',10);
                        spy(pixelsB,'r');
                        
                        for i = 1:2
                            saveas(gcf, [ p.OutputDirectory filesep filename{i}]);
                        end
                        
                        %             test = zeros(imgSize(1),imgSize(2));
                        %             idx = filoInfo(idxCurrent).([toAdd 'pixIndices']);
                        %             test(idx) = 1;
                        %
                        %             h = setFigure(imgSize(1),imgSize(2));
                        %             imshow(totalMask,[]);
                        %             hold on
                        %              spy(test,'r');
                        %              % only plot endpointCoord if actually have an appropriate exit
                        %              % flag
                        %              if (filoInfo(idxCurrent).([toAdd 'exitFlag']) >=1 && ~isnan(filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY'])(1,1)))
                        %                  scatter(filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY'])(:,1),filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY'])(:,2),'y','filled');
                        %                  testpix = zeros(imgSize(1),imgSize(2));
                        %                  testpix(pixIndices) = 1;
                        %                  spy(testpix,'y');
                        %
                        %
                        %              end
                        %
                        %              saveas(h,[currentSaveDir filesep 'Filopodia_' num2str(idxCurrent) ,'mask.eps'],'psc2');
                        
                        
                        
                    end
                else
                    filoInfo(idxCurrent).([toAdd 'params']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) =[ NaN,NaN];
                    filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'exitFlag']) = NaN;
                    filoInfo(idxCurrent).([toAdd 'resnorm']) = NaN;
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
    
    % collect pixels for all bad fits
    
end % end iType


end
function F= testfun(p,xData)
%F = p(1).*xData + p(2);
F = p(1)*(1-0.5*(1+erf((xData-p(2))/(p(3)*sqrt(2)))))+p(4);
end



