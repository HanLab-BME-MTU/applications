function [filoPercentEmbed] = GCAAnalysisExtract_percentOfActinBundlesVeilEmbedded(analInfo,filoFilterSet)
%% GCAAnalysisExtract_percentEmbed
% Calculates the percentage of total actin bundles detected with significant
% veil embedded detections. 
% Note we do attempt to detect completely embedded actin bundles in this 
% algorithm- 
% A significant exposed actin bundle detection 
% is required for the actin bundle to be considered. 
%
%
%INPUT:
%
%analInfo :  REQUIRED:  large rx1 structure (currently takes some time to load)
%                       where r is the number of frames
%                       that includes all the information one would want about the
%                       segmentation of the movie (including the filopodia information)
%                       (eventually the analInfo needs to be
%                       saved per frame I think to make faster...too big of file)
%
%filoFilterSet: PARAM:  rx2 cell array where r (row) is the number of
%                       frames and each cell is a rBundlex2 logical vector. 
%                       where rBundle is the number of Actin Bundles 
%                       detected in the frame. 
%                       Column 1 of rBundle places a 1 for each actin bundle 
%                       to be included in the analysis based on filtering 
%                       criteria
%                       Column 2 of rBundle places a 1 for each embedded
%                       detection of the actin bundle that is considered
%                       significant. 
%                       Output of GCACreateFilopodiaFilterSetWithEmbed.m 
%                       
%     
% OUTPUT:
% filoPercentEmbed:  rx1 cell array where r (row) is the number of frames
%                    and each cell holds the percentEmbedded calculation
%                    for the given frame - (output is a cell to maintain
%                    consistency with other measurements that produce
%                    distributions per frame)
%
%
%%
 
frames = length(analInfo);
filoPercentEmbed = cell(frames,1);

for iFrame = 1:length(analInfo)-1
       
    filoInfo = analInfo(iFrame).filoInfo;
    
      if ~isempty(filoInfo);
          % filter the filopodia information 
          
          filterFrameC= filoFilterSet{iFrame};
         
          % Sum the embedded filo post filtering ..but only keep those embedded with 
          % a viable external component (ie the same criteria used for total length). 
         
          
          NEmbed = sum(filterFrameC(:,2) & filterFrameC(:,1));
          
          % NTotal 
          NTot = sum(filterFrameC(:,1)); 
          
          pe = NEmbed/NTot; 
          
             filoPercentEmbed{iFrame} = pe;
    else
        filoPercentEmbed{iFrame} = [];
      end % isempty
    clear pe
          
          
end  % iFrame 
      
end 
          
        
          
    
    
 
