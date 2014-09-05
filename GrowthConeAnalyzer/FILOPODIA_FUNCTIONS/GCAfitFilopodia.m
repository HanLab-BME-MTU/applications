function [ filoInfo] = GCAfitFilopodia( filoInfo,img,paramsIn)

% Perform fits on the filopodia reconstruction in formation 
% currently designed to take in/read out information per frame. 
% was fitLinescansNewDistFinal until 20140529
p = paramsIn; 

% for now int will just be the flag for internal versus external 
if p.InternalFiloOn ==1; 
    toAdd = 'Int_'; 
else toAdd = 'Ext_'; 
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
H = fspecial('gaussian',3,p.sigma); 
imgFilt = imfilter(img,H); % for weighted averaging 

for ifilo = 1:numFilo2Fit
   idxCurrent = idx2fill(ifilo);
   
   numBack = length(filoInfo(idxCurrent).([toAdd 'pixIndicesBack']));
   numFor = length(filoInfo(idxCurrent).([toAdd 'pixIndicesFor']));
   
    % convert to distance along filo in pixels: think about a prettier way to do this.  
    xyFilo = filoInfo(idxCurrent).([toAdd 'coordsXY']);
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
   
 
if strcmpi(paramsIn.ValuesForFit,'Intensity')
   maskIndices = filoInfo(idxCurrent).Ext_maskIndices; 
   yData = imgFilt(maskIndices(:,1)); % just take the center line check to see the order of the data row 1 = base I think wherease end = tip 
   if p.SavePlots  == 1 
       if sum(isnan(maskIndices(:)))==0; 
           
       subplot(2,2,(1:2)); 
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
        if p.SavePlots ==1
            
       %     figure('visible','on');
        subplot(2,2,3:4)     
            
            linescan = scatter(distFilo,yData);
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
        starts = [max(yDataFit)-min(yDataFit), filoInfo(idxCurrent).([toAdd 'endpointCoord']),5,...
             min(yDataFit)];
         upb = [(max(yDataFit)-min(yDataFit)), length(distFiloFit),20,max(yDataFit)];
         lwb = [upb(1)/2, 0, 1, min(yDataFit)];
        
        % calculate lb and ub
        
         
        [params, resnorm,~,exitFlag] = lsqcurvefit(@testfun,starts,distFiloFit,yDataFit,lwb,upb);
        if p.SavePlots
        line([params(2),params(2)],[max(yData),min(yData)]);
        %text(params(2),(max(yDataFit)-min(yDataFit)./2)+min(yDataFit),['mean = ' num2str(params(2),3)]); 
        end 
        
        filoInfo(idxCurrent).([toAdd 'params'])  = params;
        
        % need to change now that do the distance calc... 
      %  pixIdx = round(params(2)); % here rounding the pix indices not very precise calculation should fix
        % guess this is still decent other option is to have an idx 
         delt = abs(distFilo-params(2)); % find the pixIdx (as measured along filo: 0  = the base of the filo so can get a quick xy coord for plotting. 
         pixIdx = find(delt==min(delt)); 
        
        
     %   pixIndices = filoInfo(idxCurrent).pixIndices;
        if pixIdx == 0 % indicates the fit was bad set fit parameters to NaN
            %filoInfo(idxCurrent).endpointCoordFit = NaN;
            filoInfo(idxCurrent).([toAdd 'endpointCoordFitXY']) = [NaN NaN]; 
            filoInfo(idxCurrent).([toAdd 'endpointCoordFitPix']) = NaN; 
            filoInfo(idxCurrent).([toAdd 'length']) = NaN; 
            filoInfo(idxCurrent).([toAdd 'std']) = NaN; 
            
            
        else
            
            
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
        
        filoInfo(idxCurrent).([toAdd 'resnorm']) = resnorm/distFilo(end); % divide by the number of points
        
        if exitFlag>1 && p.SavePlots == 1
            title(['Length = ' num2str(params(2).*0.216,3) ' um']) 
        
        end 
        
        
        yFit = params(1)*(1-0.5*(1+erf((distFilo-params(2))/(params(3)*sqrt(2)))))+params(4);
        hold on
        if p.SavePlots ==1
            plot(distFilo,yFit,'r')
            scatter(distFiloFit,yDataFit,'b','filled'); % color in the data specifically used for the fitting. 
%             title({['Filo' num2str(idxCurrent)];'Green Line: Endpoint Skel'; 'Blue Line: Endpoint Fit'});
         filename = ['Filopodia_' num2str(idxCurrent,'%03d') '.fig' ];
%             
            
%             if (exitFlag >= 1 && ~isnan(filoInfo(idxCurrent).endpointCoordFitXY(1,1))) 
%                 currentSaveDir = goodFitSaveDir; 
%                 maskGoodFit(filoInfo(idxCurrent).([toAdd 'pixIndices'])) =1; % add to good fit mask 
%                 
%             else 
%                 currentSaveDir = badFitSaveDir; 
%                 maskBadFit(filoInfo(idxCurrent).([toAdd 'pixIndices']))=1; 
%             end 
%                 
            saveas(linescan, [ p.OutputDirectory filesep filename]);
            close all
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
             
            close all
            
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
    
  
    
   

    
end
 
% collect pixels for all bad fits


 

end
function F= testfun(p,xData)
F = p(1)*(1-0.5*(1+erf((xData-p(2))/(p(3)*sqrt(2)))))+p(4);
end



