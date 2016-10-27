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
   
    
    xMatIn_iFrame = xMatInByFrame{iMask}; % don't filter yet you can do some filtering by displacement later
    yMatIn_iFrame= yMatInByFrame{iMask};
    xMatOut_iFrame = xMatOutByFrame{iMask};
    yMatOut_iFrame = yMatOutByFrame{iMask};
    
    
    if ~isempty(xMatIn_iFrame) % if tracks existed in current mask
        
        
        %%change to FULL TRACKS : Note before was using only the
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

