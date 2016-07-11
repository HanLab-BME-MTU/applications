function [ output_args ] = mitoticMeasureLateralClassRelatedStats(groupList,varargin)
% mitoticMeasureLateralClassRelatedStats: wrapper function, 
% calls mitoticAddPercentTransitionField and/or mitoticAddEndOnInfoWithClassifications
% 
%                                                
% % INPUT:
% groupList: (REQUIRED) : Rx3 cell array
%                        holding directories of the subRegional MT tracks
%                        you would like to include is stored.
% 
% % Parameters for classification cut-offs
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
%
% 
% 
% Parameters: Type of Calculations
%
% runPercentLat : (PARAM) : logical 
%                            Default : true 
%                            Runs : mitoticAddPercentTransitionField
%                                   Adds to corticalData.mat 
%                                   field
%                                   .params(minDispVect_orientCutOff_dispInRegion_dispDiscard_orientDiscard)
%                                          .idxLatLogical 
%                                                a rx1 logical where r is
%                                                the number of MT tracks
%                                                extracted from the
%                                                subRegion AFTER filtering 
%                                                by dispDiscard and discardOrient. 
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
%                                   
% 
% runEndOn : (PARAM) : logical 
%                           Default : true 
%                           Runs :  mitoticAddEndOnInfoWithClassifications
%                                   Adds to corticalData.mat 
%                                           .params(minDispVect_orientCutOff_dispInRegion_dispDiscard_dispOrient)
%                                               .trackInfoEndOnST
%                                                             .dwells 
%                                                                   a rx1 vector where r is
%                                                                   the number of MT subTracks
%                                                                   extracted
%                                                                   in the
%                                                                   given
%                                                                   subRegion
%                                                                   that are classified as
%                                                                   end-on given the
%                                                                   classification conditions
%                                                                   and terminate in either 
%                                                                   a MT terminal or an MT shrinkage
%                                                                   event
%                                               .trackInfoEndOnP
%                                                             .dwells  
%                                                                   a rx1 vector where r is
%                                                                   the number of MT subtracks 
%                                                                   extracted in the given subRegion 
%                                                                   that are classified as
%                                                                   end-on given the
%                                                                   classification conditions
%                                                                   and
%                                                                   terminate
%                                                                   in a
%                                                                   pause.
%                                               Note that mitoticAddPercentTransitionField
%                                               needs to be run at some point for a given 
%                                               set of classification
%                                               parameters otherwise this
%                                               function will error. (as it
%                                               needs the idxLatLogical
%                                               field). 
%                                               
%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;


% Parameters for the classification cut-offs 
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x)); 
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x)); 
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x)); 
ip.addParameter('orientDiscard',150,@(x) isnumeric(x)); 

% type of calcs 
ip.addParameter('runPercentLat',true,@(x) islogical(x)); 
ip.addParameter('runEndOn',true,@(x) islogical(x)); 

ip.parse(varargin{:});

%%
if size(groupList,2) == 2
    projList =groupList(:,2); 
else 
    projList = groupList; 
end 
    

for iProj = 1:numel(projList)
    if ~exist([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
        display(['No Growth Tracks in subRoi for' projList{iProj}])
        
    else
        s1 = load([projList{iProj} filesep 'meta' filesep 'projData.mat']);
        projData = s1.projData;
        
        
        s2 = load([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
        corticalData = s2.corticalData;
        
        if ip.Results.runPercentLat
            
        [corticalData,projData] = mitoticAddPercentTransitionField(corticalData,projData,'minDispVect',ip.Results.minDispVect,... 
            'orientCutOff',ip.Results.orientCutOff,'dispInRegion',ip.Results.dispInRegion, ...
            'dispDiscard',ip.Results.dispDiscard); 
        end 
        
        if ip.Results.runEndOn
            [corticalData,projData] = mitoticAddEndOnInfoWithClassifications(corticalData,projData,'minDispVect',ip.Results.minDispVect,... 
            'orientCutOff',ip.Results.orientCutOff,'dispInRegion',ip.Results.dispInRegion, ...
            'dispDiscard',ip.Results.dispDiscard); 
        end
         
        save([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat'],'corticalData'); 
        save([projList{iProj} filesep 'meta' filesep 'projData.mat'],'projData'); 
        
        
    end
end 

