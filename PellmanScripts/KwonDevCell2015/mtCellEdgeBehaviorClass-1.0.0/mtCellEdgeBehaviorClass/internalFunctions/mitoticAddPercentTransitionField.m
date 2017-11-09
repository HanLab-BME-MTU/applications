
function [corticalData,projData] = mitoticAddPercentTransitionField(corticalData,projData,varargin) 
% mitoticAddPercentTransitionField : Adds a field to corticalData defining the lateral vs
% end on MTs for the set of classification parameters entered and
% calculates the percent lateral transitions for the given project.  
%
% INPUT:
%  % corticalData: (REQUIRED) : structure containing the cortical data
%                               information including the lateral
%                               classification information
% 
%    projData :  (REQUIRED) : output of plusTipSubRoiExtractTracksMITOTICPACKAGE_multMasks.m
%                             contains the information corresponding to the
%                             dwell lifetime and the classification 
%                             (terminal,shrinkage,or pause) of each
%                             subTrack in the region
%    
% % Parameters for classification cut-offs : 
%   Note: The classification cut-offs for orient, dispInRegion, and dispDiscard. 
%   are defined here. You can try multiple cut-offs and it will save this 
%   in the same data structure that is called in downstream functions. However, if
%   you don't run a combination here first and try to call that combination in subsequent steps
%   (such as visualization) the workflow will fail. 
% 
% minDispVect: (PARAM) : scalar
%                        Default : 3 pixels
%                        the min MT Track Length of interest 
%                        (currently in pixels)
%                        see Supp Figure 1 A of Kwon et al. Dev Cell 2015
%
% orientCutOff: (PARAM) : scalar
%                         Default : 60 (in degrees) 
%                         Local Direction MT Growth relative to the local
%                         cell edge normal.  
%                         lateral classifications as described 
%                         in Figure S7A of Kwon et al. Dev Cell 2015
%                         
%
% dispInRegion: (PARAM) : scalar
%                         Default : 0.7 ( in um) 
%                         MT Trajectories with Cort Dist < this value will 
%                         automatically be classified as end-on. 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%
% dispDiscard : (PARAM) : scalar 
%                         Default : 0.3 (im um)
%                         MT trajectories wth MT Cort Dist < this value discarded.                         
%                         see Figure S7A of Kwon et al. Dev Cell 2015
% 
% orientDiscard : (PARAM) : scalar 
%                           Default : 150 (in Degrees)
%                           MT trajectories with with Local Direction MT
%                           Growth relative to the local cell edge normal
%                           > than this value will be discarded (likely
%                           error)
%                           
% OUTPUT: 
%  Adds to corticalData.mat field
%                                   .params(minDispVect_orientCutOff_dispInRegion_dispDiscard)
%                                          .idxLatLogical 
%                                                a rx1 logical where r is
%                                                the number of MT tracks
%                                                extracted from the
%                                                subRegion AFTER filtering 
%                                                by disDiscard. 
%                                                true indicates 
%                                                the MT trajectory has been
%                                                characterized by lateral 
%                                                behavior using the
%                                                classification conditions
%                                                input above.
%                                                
%                                           .stats 
%                                                .percentLateral
%                                                    a scalar holding the 
%                                                    the value for
%                                                    percent lateral
%                                                    classifications
%                                                    for the given subRegion
% 
%% Check input
ip = inputParser;

ip.CaseSensitive = false;

% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));
ip.addParameter('orientDiscard',150, @(x) isnumeric(x)); 

ip.parse(varargin{:});
%% Set up 
% REMOVE some extraneous fields from previous calcs 
if isfield(corticalData,'stats')
    corticalData = rmfield(corticalData,'stats');
end

if isfield(corticalData,'idxLatLogical')
    corticalData = rmfield(corticalData,'idxLatLogical');
end
%% Classify and record percent lateral transitions
classMat = corticalData.(['OrientVsDispVsDispMTVect_' num2str(ip.Results.minDispVect)]);
dispValues = classMat(:,2);
% filter by discard disp
classMat = classMat(dispValues>ip.Results.dispDiscard,:);
% filter by discard orient
classMat = classMat(classMat(:,1)<ip.Results.orientDiscard,:); 
% calculate cut-offs
idxLatLogical  = (classMat(:,1) > ip.Results.orientCutOff & classMat(:,2) > ip.Results.dispInRegion);
numLat =  sum(idxLatLogical);
percentLateral = (numLat/length(classMat(:,1)))*100;
myStrdisp= num2str(ip.Results.dispInRegion);
myStrdisp = strrep(myStrdisp,'.','pt');
myStrdiscard = num2str(ip.Results.dispDiscard);
myStrdiscard = strrep(myStrdiscard,'.','pt');
myStrdiscardOrient = num2str(ip.Results.orientDiscard); 

% give it a complete identifier and store in projData for now to keep consistent
projData.stats.(['Percent_Transitions' num2str(ip.Results.minDispVect) '_' (num2str(ip.Results.orientCutOff)) '_' myStrdisp '_' myStrdiscard '_' myStrdiscardOrient]) = percentLateral;

%     corticalData.stats.(['Percent_Transitions_' num2str(ip.Results.minDispVect) '_' (num2str(ip.Results.orientCutOff)) '_' myStrdisp '_' myStrdiscard]) = percentLateral;
corticalData.(['params' num2str(ip.Results.minDispVect) '_' num2str(ip.Results.orientCutOff) '_' myStrdisp '_' myStrdiscard '_' myStrdiscardOrient]).idxLatLogical = idxLatLogical ; % for now save in logical form in case want to plot etc
corticalData.(['params' num2str(ip.Results.minDispVect) '_' num2str(ip.Results.orientCutOff) '_' myStrdisp '_' myStrdiscard '_' myStrdiscardOrient]).stats.percentLateral = percentLateral;

end


