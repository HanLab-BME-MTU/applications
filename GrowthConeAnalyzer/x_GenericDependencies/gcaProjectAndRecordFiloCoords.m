function [filoInfo] = gcaProjectAndRecordFiloCoords(verticesEP,edgePathCoord,maxTh,maxRes,img,varargin)
%gcaProjectAndRecordFiloCoords.m
%
% Inputs the ordered coordinates cooresponding to a detected filopodia
% Assumes that the end of the ridge response is likely a mid point of an
% intensity sigmoidal decay at the tip of the filopodia along the direction
% of the filopodia length.
%
% 
%was walkfFiloForAndBack.m until 20150528
%% INPUT:

% verticesEP :  1x2 double
%   of the yx coords indicating the endpoint of the ridge
%
% edgePathCoord : rx2 double
%   where r indicates the number of coordinates along the ridge path
%   and column(1) is the y coord and column(2) is the x coord.
%   The edgePath should be ordered from endpoint r(1) to the intersection
%   with the veilStemMask. (Excluding the endpoint of the ridge)
%
% maxTh: rxc double
%   where r (rows) is the size of the original image in y and c (cols) is
%   the size of the image in x.
%   Ridge orientations as measured via the steearble filter.
%   See output of multiscalesteerableDetector.m
%
% maxRes: rxc double
%   where r (rows) is the size of the original image in y and c (cols) is
%   the size of the image in x.
%   Ridge response as measured via the steearble filter.
%   See output of multiscalesteerableDetector.m
%
% img: rxc double
%   where r (rows) is the size of the original image in y and c (cols) is
%   the size of the image in x.
%   image intensity values
%
% filoInfo : (OPTIONAL)  rx1 structure
%   where r is the number of filopodia recorded in the structure array
%   Only required filoInfo needs to be appended.
%   Default = [];
%
% veilStemMaskC : (OPTIONAL) rxc logical
%   Required only if the embeddedFlag is set to true
%   Default = [];
%
% embeddedFlag: (PARAM) logical
%   If true, indicates coordinate correspond to an embedded fit, therefore it will check
%   to make sure that the ridge detection does not span
%   beyond the veilStemMask.
%   Default = false
%
% numPixSearchForward: (PARAM) positive scalar
%   Number of pixels to project forward/record from the tip of the ridge response
%   in the direction of the filopodia (ridge)
%   Default: 10 pixels (For Filopodia outside of Veil);
%   
% propBasedOnResp : (PARAM) false ; % option to propogate forward based on the
% of the response at each pixel (might either be better or get
% noisy)
%% Check Input
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;

% REQUIRED
ip.addRequired('verticesEP');
ip.addRequired('edgePathCoord');
ip.addRequired('maxTh');
ip.addRequired('maxRes');
ip.addRequired('img');

% OPTIONAL
ip.addOptional('filoInfo',[]); 
ip.addOptional('veilStemMaskC',[],@(x) islogical(x)); 

% PARAM
ip.addParam('embeddedFlag',false,@(x) islogical(x));
ip.addParam('numPixSearchForward',10,@(x) isscalar(x));
ip.addParam('propBasedOnResp',false,@(x) islogical(x)); 

%% INITIATE
[ny nx] = size(img);
%% GET AVG RESPONSES AND INTENSITIES FOR FILO COORDS AND WALK FORWARD BASED ON STEERABLE FILTER DIRECTION
for ifilo = 1:numFilo
    %  filoCount = numFiloPrev +ifilo;
    filoCount = ifilo;
    filoCoords= edgePathCoord{ifilo};
        
        pixIdxBack = sub2ind(size(maxTh),filoCoords(:,1),filoCoords(:,2));
        
        pixIdxBack=  flipdim(pixIdxBack,1);
        filoCoords = flipdim(filoCoords,1);
        %%%%%%%%%%%%%% AVERAGE RESPONSE BACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avgResBack = nan(length(pixIdxBack),1);
        avgIntBack = nan(length(pixIdxBack),1);
        
        [yMax,xMax] = size(maxTh);
      
        indicesForMaskBack = zeros(length(pixIdxBack),5);
        for i= 1:length(pixIdxBack)
            coordsPer = zeros(5,2);
            dx = round(cos(maxTh(pixIdxBack(i))));
            dy = round(sin(maxTh(pixIdxBack(i))));
            
            coordsPer(1,1) = filoCoords(i,1); %
            coordsPer(1,2) = filoCoords(i,2);
            
            
            coordsPer(2,1) = coordsPer(1,1) + dy;
            coordsPer(2,2) = coordsPer(1,2) + dx;
            
            coordsPer(3,1) = coordsPer(1,1) -dy;
            coordsPer(3,2) = coordsPer(1,2) - dx;
            
            coordsPer(4,1) = coordsPer(2,1) +dy;
            coordsPer(4,2) = coordsPer(2,2) +dx;
            coordsPer(5,1) = coordsPer(3,1) - dy;
            coordsPer(5,2) = coordsPer(3,2) - dx;
            
            
            idxOutOfBounds = find(coordsPer(:,1) > yMax | coordsPer(:,2) > xMax |...
                coordsPer(:,1) <= 0 | coordsPer(:,2) <= 0);
            if ~isempty(idxOutOfBounds)
                coordsPer(idxOutOfBounds,:) = nan;
            end
           
            indicesForAvgBack = sub2ind(size(img),coordsPer(:,1),coordsPer(:,2));
            indicesForMaskBack(i,:) = indicesForAvgBack';
            indicesForAvgBack = indicesForAvgBack(~isnan(indicesForAvgBack));
                
            avgResBack(i) = nanmean(maxRes(indicesForAvgBack));
            avgIntBack(i) = nanmean(img(indicesForAvgBack));         
        end
      
        %%% WALK FORWARD %%%
        % project the filopodia forward based on the direction of the response
        % from the steerable filter
        
        n =  ip.Results.numPixSearchForward ; % how many pixels you will project forward
        
        % initialize (from endpoint)
        coordsFor = nan(n+1,2);
        dirFeat = nan(n+1,1);
        
        coordsFor(1,1) = verticesEP(ifilo,1); % y
        coordsFor(1,2) = verticesEP(ifilo,2); % x
        
        % NOTE direction stored in  is perpendicular to filo length
        dirFeat(1) = maxTh(coordsFor(1,1),coordsFor(1,2));
        % find which direction oriented
        
        % Define the direction forward from tip
        filoDeltX = coordsFor(1,2)-filoCoords(end,2);
        filoDeltY = coordsFor(1,1)-filoCoords(end,1);
        filoMag = sqrt(filoDeltX^2 + filoDeltY^2);
        
        dX1 = round(cos(dirFeat(1)-pi/2));
        dY1 = round(sin(dirFeat(1)-pi/2));
        
        dX2 = round(cos(dirFeat(1)+pi/2));
        dY2 = round(sin(dirFeat(1)+pi/2));
        
        cosV1Filo = (filoDeltX.*dX1 + filoDeltY.*dY1)/filoMag;
        cosV2Filo = (filoDeltX.*dX2 + filoDeltY.*dY2)/filoMag;
        
        test1 = (1-cosV1Filo);
        %test2 =(1-cosV2Filo);
        
        if abs(test1) < 1
            dX = dX1;
            dY = dY1;
            add = -pi/2;
        else
            dX = dX2;
            dY = dY2;
            add = pi/2;
        end
        
        coordsFor(2,2) = coordsFor(1,2) + dX;
        coordsFor(2,1) = coordsFor(1,1) + dY;
        
        
        
        % use the direction at tip for all projected values imshow(img,[]);
        
        for i = 2:n
            if (round(coordsFor(i,1))<=0 || round(coordsFor(i,2))<=0 || round(coordsFor(i,1))>ny || round(coordsFor(i,2))>nx)%
                % break if hit end of image
                
                idx = ~isnan(coordsFor(:,1)); % shorten the coords
                coordsFor = coordsFor(idx,:);
                coordsFor = coordsFor(1:end-1,:); % delete the problem coord
                
                
                break
                
            else
                if ip.Results.propBasedOnResp == 1
                    % recalculate dir based on new coord
                    
                    dirFeat(i) = maxTh(round(coordsFor(i-1,1)),round(coordsFor(i-1,2)));
                    dX = round(cos(dirFeat(i)+add));
                    dY = round(sin(dirFeat(i)+add));
                    coordsFor(i+1,2) = coordsFor(i,2) + dX;
                    coordsFor(i+1,1) = coordsFor(i,1) + dY;
                else  % project coordinates only based on the direction of the
                    % response at the tip of the filo
                    coordsFor(i+1,2) = coordsFor(i,2) + dX;
                    coordsFor(i+1,1) = coordsFor(i,1) + dY;
                end
            end
        end
        
        idxOutOfBounds = find(coordsFor(:,1) > yMax | coordsFor(:,2) > xMax |...
            coordsFor(:,1) <= 0 | coordsFor(:,2) <= 0);
        if ~isempty(idxOutOfBounds)
            coordsFor(idxOutOfBounds,:) = nan ;
        end
        
        
        pixIdxFor =  sub2ind(size(maxTh),round(coordsFor(:,1)),round(coordsFor(:,2)));
        
      
        % get average response 2 pixels on each sideof coordinate
        % for now just do direction feature of tip %%
        % SHOULD CHANGE TO INCORPORATE INFORMATION OF SCALES!!!
        avgResFor = nan(n+1,1);
        avgIntFor = nan(n+1,1);
        % stupid way to do this think about changing.
        coordsPer = zeros(5,2);
        
        indicesForMaskFor = zeros(length(pixIdxFor),5);
        
        for i= 1:length(pixIdxFor)
            
            dx = round(cos(dirFeat(1)));
            dy = round(sin(dirFeat(1)));
            
            coordsPer(1,1) = coordsFor(i,1);% maybe add a third d here: therefore save as
            coordsPer(1,2) = coordsFor(i,2); % first d is the center coord (xy) - maybe just make  coordsPer(i,1,1) etc
            % where i is the ID along the filo - though the easiest way to do it for
            % visualization is just dilation....also could just get the dilataion
            % along each point and weight center high surrounding lower - in this
            % manner front and behind will be taken in the signal ... could likely
            
            
            
            coordsPer(2,1) = coordsPer(1,1) + dy; % first set of flanking coords
            coordsPer(2,2) = coordsPer(1,2) + dx;
            
            coordsPer(3,1) = coordsPer(1,1) -dy;
            coordsPer(3,2) = coordsPer(1,2) - dx;
            
            coordsPer(4,1) = coordsPer(2,1) +dy; % second set of flanking coords
            coordsPer(4,2) = coordsPer(2,2) +dx;
            coordsPer(5,1) = coordsPer(3,1) - dy;
            coordsPer(5,2) = coordsPer(3,2) - dx;
            
            idxOutOfBounds = find(coordsPer(:,1) > yMax | coordsPer(:,2) > xMax |...
                coordsPer(:,1) <= 0 | coordsPer(:,2) <= 0);
            if ~isempty(idxOutOfBounds)
                coordsPer(idxOutOfBounds,:) = nan ;
            end
                 
            indicesForAvg = sub2ind(size(img),coordsPer(:,1),coordsPer(:,2));
            indicesForMaskFor(i,:)  = indicesForAvg';
            indicesForAvg = indicesForAvg(~isnan(indicesForAvg));
            
            
            avgResFor(i) = nanmean(maxRes(indicesForAvg));
            avgIntFor(i) = nanmean(img(indicesForAvg));
        end
        
        
        
        %%% SAVE RESP AND INTENSITY INFO %%%
        
        avgResFilo = [avgResBack;avgResFor] ;
        avgIntFilo = [avgIntBack;avgIntFor];
        %% Make sure to filter out pixels that extend beyond the body Mask (however need to make sure pixIdxFor is in same order
        if embeddedFlag == true % you need to filter to make sure the pixels don't exceed the body estimation
            pixelsForBody = find(veilStemMaskC==1);
            pixIdxFor = intersect(pixIdxFor,pixelsForBody,'stable');
        end
        pixIndices = [pixIdxBack; pixIdxFor];
    
    
    if embeddedFlag == true;
        toAdd = 'Int_';
    else toAdd = 'Ext_';
    end
    
    % think I want to save everything then filter to fit.
    filoInfo(filoCount).([toAdd 'response' ]) = avgResFilo;
    filoInfo(filoCount).([toAdd 'intensities' ]) = avgIntFilo;
    filoInfo(filoCount).([toAdd 'maskIndices']) = [indicesForMaskBack; indicesForMaskFor];
    %filoInfo(filoCount).([toAdd 'maskPostFit']) =   ;
    filoInfo(filoCount).([toAdd 'dirAtTip' ]) = dirFeat(:);
    filoInfo(filoCount).([toAdd 'pixIndices' ]) =  pixIndices; % both forward and back
    [y,x] = ind2sub(size(img),pixIndices);
    filoInfo(filoCount).([toAdd 'coordsXY' ]) =[x,y]; % all potential filoCoords;
    
    filoInfo(filoCount).([toAdd 'endpointCoord' ]) = verticesEP(ifilo,:); % based on response
    filoInfo(filoCount).([toAdd 'pixIndicesFor' ]) = pixIdxFor; % for plotting
    filoInfo(filoCount).([toAdd 'pixIndicesBack' ]) = pixIdxBack; % for plotting
    filoInfo(filoCount).([toAdd 'vectFilo']) = [x(end,1)-x(1,1),y(end,1)-y(1,1)];
    
    
end %ifilo

end


