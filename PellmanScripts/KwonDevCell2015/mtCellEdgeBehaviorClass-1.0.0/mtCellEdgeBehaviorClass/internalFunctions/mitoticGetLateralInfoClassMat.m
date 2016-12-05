function [OrientVsDispVsDispMTVect] = mitoticGetLateralInfoClassMat(subRoiEdgeDir,varargin)
% INPUT: mitoticGetLateralInfoClassMat: measures the orientation/disp/dispMTvect for a single subroi project  

%  subRoiEdgeDir: the cortical subRoi File 
%        minDispVect: a scalar specifying the min disp vect to use for the
%        orient calc 
% 
%From plusTipTracker output : loads :
%                              /roi_X/meta/projData.mat
%
%From mitoticWrapperThresholding output: loads:
%                              /roi_X/masks  (These are perFrame masks of the cell edge)
%                              loads in mitoticGetClassificationMat
%
% INPUT:
% subRoiEdgeDir: (REQUIRED) : character array 
%                            subRoiDirectory where the projData.mat is stored                     
% 
% minDispVect: (PARAM) : scalar
%                        Default : 3 (in pixels)
%                        minimum displacement along the microtuble starting
%                        from the trajectory endpoint used for the orientation
%                        calculation. The first point along the MT trajectory
%                        to be > than this value is used for the calculation
%                        (Note in the future one could instead choose to smooth the
%                        trajectory and try to get an instantaneous
%                        orientation value at a set distance along the MT
%                        from this interpolation- The only caution would be
%                        how to choose the most appropriate smoothing parameter:
%                        We indeed thought about implementing this, but decided
%                        due to time constraints and the fact these
%                        approximate measurements
%                        were stable to small perturbation in the
%                        minDispVector used that changing the measurement approach wasn't
%                        necessary as the biological conclusion would be
%                        likely be maintained.
%
% % TSOverlays (PARAM)    : logical
%                         Default : true 
%                         Flag to make the troubleshooting (TS) overlays 
%                         that plot the classification measurements per track 
%                         extracted: note this will
%                         incease the time it take to run the results but 
%                         can be very helpful if you want to reset the
%                         classification cut-offs to target a specific biological 
%                         behavior.  
%% Check Input 

ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('minDispVect',3,@(x) isnumeric(x));

ip.addParameter('TSOverlays',true,@(x) islogical(x)); 

ip.parse(varargin{:});


%% 


listOfMasks = mitoticSearchFiles('.tif',[],[subRoiEdgeDir filesep 'masks'],0);

s= load([subRoiEdgeDir filesep 'meta' filesep 'projData.mat']);
projData = s.projData;

dispByFrame = projData.dispByFrame;
velsByFrame = projData.velsByFrame;
xMatInByFrame = projData.xMatIn_ByFrame;
yMatInByFrame = projData.yMatIn_ByFrame;
%         trckIdxInMultMasks = projData.trckIdxInMultMasks;
%         trckIdxInMultMasks_byFrame =  projData.trckIdxInMultMasks_ByFrame ;
xMatOutByFrame = projData.xMatOut_ByFrame;
yMatOutByFrame = projData.yMatOut_ByFrame;
% put the tracks back together 


for iMask = 1:length(listOfMasks)-1
    % get all displacement values for each mask (the displacement in the
    % mask saved before)
    dispValues = dispByFrame{iMask};
    velValues = velsByFrame{iMask};
    %
%                idxDiscard_perFrame{iMask} = find(dispValues<distCutoffDiscard); % for plotting
    %
%                 if ~isempty(idxDiscard_perFrame{iMask})
%                     nDiscard{iMask}= size(idxDiscard_perFrame{iMask},1);
%                     xCoordsDiscard{iMask} =  xMatInByFrame{iMask}(dispValues<distCutoffDiscard,:);
%                     yCoordsDiscard{iMask} = yMatInByFrame{iMask}(dispValues<distCutoffDiscard,:);
%     
%                 else
%                     nDiscard{iMask} = 0;
%                     xCoordsDiscard{iMask} = [];
%                     yCoordsDiscard{iMask} = [];
%                 end % isempty idxDiscard
    
    % put the track back together 
   
    
    xMatIn_iFrame = xMatInByFrame{iMask}; % don't filter yet you can do some filtering by displacement later
    yMatIn_iFrame= yMatInByFrame{iMask};
    xMatOut_iFrame = xMatOutByFrame{iMask};
    yMatOut_iFrame = yMatOutByFrame{iMask};
    
    
    if ~isempty(xMatIn_iFrame) % if tracks existed in current mask
        
        
        %%change to FULL TRACKS : 2013_09_21 Note before was using only the
        % part of the track IN the region.
        % this will give NaNs if the minimum displacement vector excedes that
        % minimum cut-off.
        % therefore put the tracks back together here and get the orientation from the
        % full track
        %
        
        % put the track back together: silly but need to do this because never
        % saved full track in data struct didn't think I would need it.
        % now very important now to make sure uniform filtering of minimum displacements for
        % discarding track information
        xMatIn_iFrame_ones = swapMaskValues(xMatIn_iFrame,nan,1);
        xMatOut_iFrame_ones = swapMaskValues(xMatOut_iFrame,nan,1);
        xMatFull_iFrame_ones = xMatIn_iFrame_ones.*xMatOut_iFrame_ones;
        xMatFull_iFrame = swapMaskValues(xMatFull_iFrame_ones,1,nan);
        
        yMatIn_iFrame_ones = swapMaskValues(yMatIn_iFrame,nan,1);
        yMatOut_iFrame_ones = swapMaskValues(yMatOut_iFrame,nan,1);
        yMatFull_iFrame_ones = yMatIn_iFrame_ones.*yMatOut_iFrame_ones;
        yMatFull_iFrame = swapMaskValues(yMatFull_iFrame_ones,1,nan);
        
        
        
        %     dispValuesFilt = dispValues(dispValues>0.3);
        %     xMatIn_iFrame= xMatIn_iFrame(dispValues>0.3,:);
        %     yMatIn_iFrame = yMatIn_iFrame(dispValues>0.3,:);
        %
        % Get the orientation of the MT trajectories relative to the cell 
        % edge 
        [incidence,dispMTVect] =  mitoticGetMTTrackOrientRelativeToCellEdge(projData,xMatFull_iFrame,...
            yMatFull_iFrame,iMask,'minDispVect', ip.Results.minDispVect,'TSOverlays',ip.Results.TSOverlays);
       
        incidenceByFrame{iMask} = rad2deg(incidence);
        dispMTVectByFrame{iMask} = dispMTVect;
     %   dispValuesFiltByFrame{iMask}  = dispValues;
    else
        incidenceByFrame{iMask} =[];
        dispMTVectByFrame{iMask} = [];
    end % isempty(xMatIn_iFrame)
    
    %% Curvature calc
%               if ~isempty(xMatIn_iFrame) 
%                    figure; 
%               % just do it stpupid way for now get curvature. 
%              for iTrack = 1:length(xMatIn_iFrame(:,1))
%                x = xMatIn_iFrame(iTrack,:); 
%                x = x(~isnan(x)); 
%              
%                y = yMatIn_iFrame(iTrack,:);
%                  y = y(~isnan(y));
% %                  x= 1:10; 
% %                  R  = 20; 
% %                  y = -sqrt(R^2-x.^2); 
%                s = spline(x,y);
%                
%                s1 = s; 
%                s1.order= s.order-1; 
%                s1.coefs = bsxfun(@times,s.coefs(:,1:end-1),s1.order:-1:1); 
%                
%                
%                s2= s1; 
%                s2.order = s1.order-1; 
%                s2.coefs= bsxfun(@times,s1.coefs(:,1:end-1),s2.order:-1:1); 
%               
%                xInt = linspace(x(1),x(end)); % try to end for now 
%                A1 = ppval(s1,x); 
%                A2 = ppval(s2,x); 
%                
%                kappa= A2./(1+A1.^2).^(3/2);
%                kappaPerTrack(iTrack) =  mean(kappa); 
%                
%              
%                % find the circle that goes through point x,y and has radius
%                % x. 
%                
%                % sanity check 
%                rPerTrack(iTrack) = 1/mean(kappa); 
%                
%                imshow(maskC,[]) ; 
%                hold on 
%                
%                plot(x,y,'y','Linewidth',2); 
%                
%                hold on 
%                
%                %plot(x,r,'b','Linewidth',2); 
%                
%              end 
%               
%               kappaByFrame{iMask} = kappa; 
%               clear kappa
%               if ~isdir([subRoiEdgeDir filesep 'sanityCheckKappa']) 
%                 mkdir([subRoiEdgeDir filesep 'sanityCheckKappa']); 
%               end 
%               saveas(gcf,[subRoiEdgeDir filesep 'sanityCheckKappa' filesep num2str(iMask,'%03d') '.tif']); 
%               
%               else 
%                   kappaByFrame{iMask} = []; 
%               end 
%%
clear xMatIn_iFrame yMatIn_iFrame xMatOut_iFrame yMatOut_iFrame
end % iMask 
%    % these should now be no longer filtered by min displacement. 
     dispValuesAll = vertcat(dispByFrame{:});
     orientValuesAll =vertcat( incidenceByFrame{:});
     dispMTVectAll = vertcat(dispMTVectByFrame{:}); 
     % 
OrientVsDispVsDispMTVect = [orientValuesAll  dispValuesAll dispMTVectAll];

end

