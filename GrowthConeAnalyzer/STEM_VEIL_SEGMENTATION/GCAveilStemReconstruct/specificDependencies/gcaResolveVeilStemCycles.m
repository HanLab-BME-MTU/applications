function [resolvedVeilStemMask,cycleFlag,TSFigs] = gcaResolveVeilStemCycles(backbone2Dil,veilStemNodeMask,backboneInfo,img,varargin)
% gcaResolveVeilStemCycles : this function tests for and resolved veil stem
% cycles 


%% INPUT 
% veilStemMask: (REQUIRED)  an rxc logical array (binary mask) 
%               where r is the height (ny)
%               of the input img and c is the width (nx)
%               marking the current binary veil/stem mask for which to test
%               for cycles. 
%
% backbone2Dil : (REQUIRED)  an rxc logical array (binary mask) 
%
%
% veilStemNodeMask: (REQUIRED) an rxc logical array (binary mask) 
%
% TSOverlays : (PARAMS) 
% 
%% INPUT PARSER
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;
%REQUIRED
ip.addRequired('backbone2Dil');
ip.addRequired('veilStemNodeMask');
ip.addRequired('backboneInfo'); 
ip.addRequired('img'); 
%ip.addRequired('BBScaleC'); 


%ip.addOptional('img',[]); 

% PARAMETERS
ip.addParameter('TSOverlays',true,@(x) islogical(x));

ip.parse(backbone2Dil,veilStemNodeMask, backboneInfo,varargin{:});
p = ip.Results;

%% Initiate 
if p.TSOverlays == true 
    iFig = 1; 
end 
[ny,nx] = size(backbone2Dil); 
TSFigs = []; 
%% TEST FOR CYCLES AND CORRECT
        dilBBMask = imdilate(backbone2Dil,strel('disk',4));
        %[~,~,~,scaleMapFine] = gcaMultiscaleSteerableDetector(img,4,'sigmaArray',[1:0.5:6]);
        %dilBBMask =  gcaImdilateWithScale(backbone2Dil,scaleMapFine,[1:0.5:6]); 
        
        
        
         
        
        fullMask = dilBBMask | veilStemNodeMask;
        
        if p.TSOverlays == true
            
            
            TSFigs(iFig).h = setFigure(nx,ny,'off');
            TSFigs(iFig).name = 'Ridge Radius Estimation';
            if ~(isempty(ip.Results.img)) ;
                
                imshow(-img,[]);
                
                hold on
            end 
                
                spy(backbone2Dil,'k');
%                 idx =  find(backbone2Dil); 
%                   values = scaleMap(idx); 
%                   cmap = colormap('jet',128); 
%                 [y,x] = ind2sub([ny,nx],idx); 
%                 scatter(x(:),y(:),values'filled'
                
                roiYX = bwboundaries(dilBBMask);
               
                roiYX2 = bwboundaries(fullMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','b'),roiYX2);
                 cellfun(@(x) plot(x(:,2),x(:,1),'color','r'),roiYX);
                text(5,5,'NMS from 1st Scale Integration to Dilate','Color','k'); 
                text(5,15,'Radius Estimation','Color','r'); 
                text(5,25,'Final Mask','Color','b'); 
                
                
            
            iFig = iFig+1; 
        end
        
        if p.TSOverlays == true 
           TSFigs(iFig).h = setFigure(nx,ny,'off');
           TSFigs(iFig).name = 'Scale Integration Original';
           imagesc(backboneInfo.scaleMapLarge); 
           colorbar
           iFig = iFig+1;   
        end 
        
%         if p.TSOverlays == true 
%             TSFigs(iFig).h = setFigure(nx,ny,'off'); 
%             TSFigs(iFig).name = 'Scale Integration Fine';
%             imagesc(scaleMapFine); 
%             colorbar
%             iFig = iFig+1; 
%             
%             
%         end 
%         
%         if p.TSOverlays == true 
%            TSFigs(iFig).h = setFigure(nx,ny,'off'); 
%            TSFigs(iFig).name = 'Scale Integration Fine Plus Overlay'; 
%            imagesc(scaleMapFine); 
%            hold on 
%            cellfun(@(x) plot(x(:,2),x(:,1),'color','w','Linewidth',2),roiYX2); % fullMask
%            cellfun(@(x) plot(x(:,2),x(:,1),'color','w','Linewidth',2),roiYX); % dilated region
%            roiYX3 =  bwboundaries(backbone2Dil);
%            cellfun(@(x) plot(x(:,2),x(:,1),'color','w','Linewidth',2),roiYX3); 
%            colorbar 
%         end 
%         
        
        % take largest cc and fill holes
        fullMask = logical(getLargestCC(fullMask));
        prefill = fullMask; % added 20140819
        fullMask = imfill(fullMask,'holes');
        cycleFlag = 0;
        if ~isequal(prefill,fullMask);% you have cycle.
            cycleFlag = 1;
            display('you have a cycle');
            % deconstruct the body as a graph
            % label the new body mask
            labelsBody = bwlabel(veilStemNodeMask);
           
            % Problem 20141009 becomes that dilation of 4 is can cross the
            % small body pieces resulting in a merging of the two paths
            % quick solution is use a smaller dilation for the label making
            % In the end we want to get a better estimate for the redilation
            % of the paths anywa
            %dilBBForLabels =  imdilate(backbone2Dil,strel('disk',2));
            
            
            
            %% (FIX)  Small test to mask sure dilation does not merge edge paths
            % the idea here is want to dilate for the intensity
            % integration reponse metrics but don't want this to be -
            % think about reworks in this coding for the final release as
            % it is a bit rough.
            dilBBMask = imdilate(backbone2Dil,strel('disk',4));
            CCPreDil = bwconncomp(backbone2Dil);
            CCEdges = bwconncomp(dilBBMask);
            stopFlagLowerDil = CCPreDil.NumObjects  > CCEdges.NumObjects;
            countDilDec = 1;
            while stopFlagLowerDil >0
                dilBBMask = imdilate(backbone2Dil,strel('disk',4-countDilDec));
                CCEdges = bwconncomp(dilBBMask);
                stopFlagLowerDil = CCPreDil.NumObjects  > CCEdges.NumObjects;
                countDilDec = 1 + countDilDec;
            end % while
            % END TEST 1
            %% 2nd Test for problems in the case dilation was too large.
            labelsC = bwlabel(dilBBMask);
            % for each label get the pixels that overlap body labels%
            CCEdges = bwconncomp(labelsC);
            edges = cellfun(@(x) unique(labelsBody(x)),CCEdges.PixelIdxList,'uniformoutput',0);
            % dilate each piece and give
            %conForLabel = imdilate(dilBBMask,strel('disk',3));
            edges =  cellfun(@(x) x(x~=0)',edges,'uniformoutput',0);
            
            % test for problems in the dilation
            numVertices = cellfun(@(x) length(x) ,edges);
            stopFlagLowerDil = sum(numVertices>2);% initiate stop flag
            % Also check to make sure that CC before dilation remains
            % consistent
            
            countDilDec = 1; % initiate count
            %if sum(numVertices > 2) ~=0 % test for problem cases where the dilation of 4 was too large
            % and spanned the body the small node...
            while stopFlagLowerDil >0
                dilBBMask = imdilate(backbone2Dil,strel('disk',4-countDilDec));
                labelsC = bwlabel(dilBBMask);
                CCEdges = bwconncomp(labelsC);
                edges = cellfun(@(x) unique(labelsBody(x)),CCEdges.PixelIdxList,'uniformoutput',0);
                % dilate each piece and give
                %conForLabel = imdilate(dilBBMask,strel('disk',3));
                edges =  cellfun(@(x) x(x~=0)',edges,'uniformoutput',0);
                
                
                numVertices = cellfun(@(x) length(x),edges);
                stopFlagLowerDil = sum(numVertices>2);
                countDilDec = 1 + countDilDec;
            end % while    
            
%% Calculate Scores For the Edges : RESPONSE SCORE 
            
            % most probable paths have high response steerable filter
            % response values.
            responseMap = backboneInfo.maxNMSLarge; % currently do NOT save the full
            % response in the backboneInfo need to check if this is more
            % helpful.
            
            resScore = cellfun(@(x) responseMap(x),CCEdges.PixelIdxList,'uniformoutput',0);
            resScore = cellfun(@(x) mean(x(x~=0)),resScore); % take out zer values from NMS.
            
            resScore = resScore./max(resScore);
            
            
            % Make Trouble Shooting Plots for the Steerable Filter Response
            % Score (Mean of Path)
            if p.TSOverlays == true
               
                TSFigs(iFig).h =  setFigure(nx,ny,'off');
                TSFigs(iFig).name = 'Response Score'; 
                imagesc(responseMap)
                hold on
                text(20,20,'Trouble Shoot Response Score');
                colorbar
                % plot the outline of the paths considered
                roiYX = bwboundaries(dilBBMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
                roiYXPieces = bwboundaries(veilStemNodeMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYXPieces);
                % document response score for each path
                centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
                arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(resScore(i),3),'color','w'),1:length(centers));
                
                iFig = iFig+1; 
            end % p.plots == 1
            
 %% Calculate Scores For the Edges : SCALE SCORE       
            
            % most probable paths have larger scale ridges (at least when the
            % scale estimate is working correctly - currenlty there seems to
            % think the max response is defaulting to max scale tested due
            % to either potential bug or something have to actually work
            % out.
            scaleMap = backboneInfo.scaleMapLarge ;
            scaleScore =  cellfun(@(x) mean(scaleMap(x)),CCEdges.PixelIdxList);
            scaleScore = scaleScore./max(scaleScore); % make between 0 and 1
            
            % Make Trouble Shooting Plots for the Steerable Filter Response
            % Score (Mean of Path)
            if p.TSOverlays == true
%                 TBScalePath =  [p.OutputDirectory filesep 'TroubleShootScaleScores'];
%                 if ~isdir(TBScalePath)
%                     mkdir(TBScalePath)
%                 end
                TSFigs(iFig).h = setFigure(nx,ny,'off');
                TSFigs(iFig).name = 'ScaleScore'; 
                imagesc(scaleMap)
                hold on
                text(20,20,'Trouble Shoot Scale Score');
                colorbar
                % plot the outline of the paths considered
                roiYX = bwboundaries(dilBBMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYXPieces);
                % document response score for each path
                centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
                arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(scaleScore(i),3),'color','w'),1:length(centers));
%                 saveas(gcf,[TBScalePath filesep num2str(iFrame,'%03d') '.tif']);
%                 saveas(gcf,[TBScalePath filesep num2str(iFrame,'%03d') '.fig']);
%                 close gcf
            end % p.plots == 1
            
%% Calculate Scores For the Edges : INTENSITY SCORE 
            intScore = cellfun(@(x) img(x),CCEdges.PixelIdxList,'uniformoutput',0);
            intScore = cellfun(@(x) mean(x(x~=0)),intScore); % take out zer values from NMS.
            
            intScore = intScore./max(intScore);
            
            % Make Trouble Shooting Plots for the Steerable Filter Response
            % Score (Mean of Path)
            if p.TSOverlays == 1
               % TBIntPath =  [p.OutputDirectory filesep 'TroubleShootIntensity'];
                %if ~isdir(TBIntPath)
                 %   mkdir(TBIntPath)
                %end
                TSFigs(iFig).h = setFigure(nx,ny,'off');
                TSFigs(iFig).name = 'IntensityScore'; 
                
                imshow(-img,[])
                hold on
                text(20,20,'Trouble Shoot Scale Score');
                colorbar
                % plot the outline of the paths considered
                roiYX = bwboundaries(dilBBMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','y'),roiYX);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','r'),roiYXPieces);
                % document response score for each path
                centers  =  regionprops(CCEdges,'centroid'); % returns a struct with field centroid
                arrayfun(@(i) text(centers(i).Centroid(1),centers(i).Centroid(2),num2str(intScore(i),3),'color','w'),1:length(centers));
               % saveas(gcf,[TBIntPath filesep num2str(iFrame,'%03d') '.tif']);
               % saveas(gcf,[TBIntPath filesep num2str(iFrame,'%03d') '.fig']);
                iFig = iFig +1; 
            end % p.plots == 1
                
            finalScore = scaleScore + resScore +intScore; % for now just use a scale Score
%% TEST FOR PARALLEL EDGES
            edgesStr  = cellfun(@(x) num2str(x),edges,'uniformoutput',0); % put into string to use unique
            % find all parallel edges - stupid way for sure but lets do it
            % for now.
            [test,ic,iR]  = unique(edgesStr);
            if length(edgesStr) ~= length(test); % test for repeats.
                display('Parallel Connections To Same Node Found: Fixing');
                
                uiR = unique(iR); % get the indexes of the potential repeats
                idxDiscard = cell(length(uiR),1);
                for iPotRepeat = 1:length(uiR)
                    
                    repeatTest = sum(iR == uiR(iPotRepeat)); % problem if there are two repeats Fixed 20141129.
                    if repeatTest > 1; % parallel edge
                        % get the score for each edge and choose the max % or
                        % could possibly just
                        idxTest  = find(iR == uiR(iPotRepeat)); % indices of th
                        % getScores or repeats
                        scoresRepeats = finalScore(iR==uiR(iPotRepeat));
                        maxScoreInGroup = max(scoresRepeats);
                        idxDiscard{iPotRepeat} = idxTest(scoresRepeats~=maxScoreInGroup);
                    end  % repeat test
                end %
                idxDiscard = vertcat(idxDiscard{:}); % changed from horzcat... 20141207
                % if check
                check = 1;
                if check == 1
                    imshow(img,[]);
                    hold on
                    edgeMask = zeros(size(img));
                    edgeMask(vertcat(CCEdges.PixelIdxList{:})) = 1;
                    spy(edgeMask,'b');
                    edgeMask(vertcat(CCEdges.PixelIdxList{idxDiscard}))= 0 ;
                    spy(edgeMask,'r');
                end
                
                % discard that edge
                edges(idxDiscard) = [];
                dilBBMask(vertcat(CCEdges.PixelIdxList{idxDiscard}))=0;
                CCEdges.PixelIdxList(idxDiscard) = [];
                CCEdges.NumObjects = CCEdges.NumObjects - length(idxDiscard);
                finalScore(idxDiscard) = [] ;
                % CCEdges.NumObjects = CCEdges.
                %end % repeat test
                % end % iPotRepeat
                % for iRepeat = 1:length(uiR)
                % end
                %if length(iRepeat
                
                % seem to have a problem with putting parallel edges
                % set up weights
                
                
            end
            prefill = (veilStemNodeMask | dilBBMask);
            % try again to fix the cycle
            
            
            % retest for cycles.
            
            % set up as a min span tree (need to make a sparse array)
            
            % check for parallel edges by finding the edges repeats in
            % the edge map
            % get all repeats. choose the max score  path here.
            %Try to test again
            
            fullMask = imfill(prefill,'holes');
            diffMask= fullMask-prefill;
            numBodyNodes = cellfun(@(x) length(x),edges); % as I don't treat
            % the surrounding frame pixels as a body "node" sometimes these
            % an path can connect to only 1 body piece this will cause the
            % following steps to crash - therefore select for only those
            % paths that connect two well-defined body pieces)
            edges = edges(numBodyNodes ==2) ;
            finalScore= finalScore(numBodyNodes==2);
            CCEdges.PixelIdxList = CCEdges.PixelIdxList(numBodyNodes==2);
            CCEdges.NumObjects = CCEdges.NumObjects - sum(numBodyNodes==2) ;
            
            
            if (~isequal(fullMask,prefill) && sum(diffMask(:))>2 && ~isempty(edges)); % make the size
                % slightly larger as cycle test currently based on simple
                % fill criterion therefore not that stable. see if can make
                % more stable.
                % if ~isempty(edges) % again simply check if it is a viable cycle-
                % again this is a weakness in the way I implemented
                % this..
                cycleFlag = 2; % non-parallel cycles
                % perform minspanning tree..
                display('cylces NOT due to parallel paths: performing minspantree');
                % put the information into a sparse mtrix.
                vect1 = cellfun(@(x) x(1),edges);
                vect2 = cellfun(@(x) x(2),edges);
                % make so lower score more favorable but do not have a zero weight
                % as a sparse array as will remove that edge completely
                finalScore = max(finalScore) - finalScore + 0.01;
                UG = sparse([vect1 vect2], [vect2 vect1], [finalScore finalScore] ); % need to make undirected.
                
                gFinal =  graphminspantree(UG);
                weightsFinal = full(gFinal(gFinal~=0));
                % delete edges given min span tree output.
                weightDelete = setdiff(finalScore,weightsFinal);
                idxDiscard   = arrayfun(@(i) find(finalScore == weightDelete(i)),1:length(weightDelete));
                
                
                if p.TSOverlays == true
                    % create the directory if doesn't exist
                    %outPathTBTree = [p.OutputDirectory filesep 'MinSpan Tree Weights'];
%                     if ~isdir(outPathTBTree)
%                         mkdir(outPathTBTree)
%                     end
                    
                    TSFigs(iFig).h = setFigure(nx,ny,'off'); 
                    TSFigs(iFig).name = 'MinSpanTreeWeights';
                    imshow(-img,[]);
                    hold on
                    
                    spy(dilBBMask,'b');
                    hold on
                end
                % delete the edge from the dilated backbone mask
                dilBBMask(vertcat(CCEdges.PixelIdxList{idxDiscard}))= 0 ; %
                if p.TSOverlays == 1
                    
                    spy(dilBBMask,'r');
                    cellfun(@(x) plot(x(:,2),x(:,1),'color','b'),roiYXPieces);
                    %saveas(gcf, [outPathTBTree filesep num2str(iFrame,'%03d') '.tif']);
                    %saveas(gcf,[outPathTBTree filesep num2str(iFrame,'%03d') '.fig']);
                    iFig = iFig +1; 
                end
                
                
                % discard the appropriate edges from the list (ccs of the dilBBMask mask )
                CCEdges.PixelIdxList(idxDiscard) = [];
                CCEdges.NumObjects = CCEdges.NumObjects - length(idxDiscard);
                
                fullMaskPreFill = veilStemNodeMask | dilBBMask;
                % test for cycles one last time. - if still have cycles can flag
                fullMask = imfill(fullMaskPreFill,'holes');
                if ~isequal(fullMaskPreFill,fullMask)
                    display('You still have a hole that the algorithm cannot currently resolve: this frame will be flagged as a low confidence neurite body reconstruct');
                    cycleFlag= 3;
                end
                
                
            end % second check for cycles (before minSpanTree)
            
            
            
        end % all tests for cycles      
        resolvedVeilStemMask = fullMask;
end

