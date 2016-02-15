function [filoLengthToVeil] = GCAAnalysisExtract_filoLength(analInfo,filoFilterSet,filoRegion,includeAllInt)
%% GCAAnalysisExtract_filoLengthToVeil
% Collects Filopodia Lengths to the Veil for an entire movie for a
% Filtered Set of Filopodia
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
%filoFilterSet: PARAM:  rx1 cell array where r (row) is the number of
%                       frames and each cell is a r1x1 logical vector
%                       that places a 1 for each filopodia to be included
%                       in the analysis.
%
%umPerPixel     PARAM:  um per pixel % add as a parameter.
%
%filoRegion     PARAM: filoRegion: character 'Int_', 'Ext_', or 'Tot'
%                      specifying the part of the actin bundle "filopodia"
%                      you want to extract (DEFAULT: 'Ext_')
%
% OUTPUT:
% filoLengthToVeil:  rx1 cell array where r (row) is the number of frames
%                    and each cell holds the distribution of filo lengths
%                    for the given frame
%
%

%%

if nargin<4 
   includeAllInt = true;  
end


% 

if length(analInfo) == 1
    nFrames = 1; 
else 
    nFrames = length(analInfo)-1; 
end 
filoLengthToVeil = cell(nFrames,1);

for iFrame = 1:nFrames
    
    
    filoInfo = analInfo(iFrame).filoInfo;
    
    if ~isempty(filoInfo);
        % currently 2 columns of the filter - one for the external and one for the 
        % internal based on fitting criteria. 
        filterFrameAll= filoFilterSet{iFrame};
        
        if strcmpi(filoRegion,'Int_'); 
            filterFrameC = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) == 1); 
        else 
            filterFrameC = filterFrameAll(:,1); 
        end 
        
        filoInfoFilt  = filoInfo(filterFrameC);
        
        % collect lengths
        if ~strcmpi(filoRegion,'Tot');
            lengths =  vertcat(filoInfoFilt(:).([(filoRegion) 'length'])).*.216; % add as a parameter
            if (strcmpi(filoRegion,'Int_') && includeAllInt);
                lengths(isnan(lengths)) = 0 ;        
                
                
            end 
        else
            % for now filter out all the internal filo that do not pass the
            % fitting criteria 
            
             filterInt = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) ==0 ); 
             filterInt = filterInt(filterFrameC);  
             
            lengthsInt =  vertcat(filoInfoFilt(:).Int_length).*.216;
            lengthsExt = vertcat(filoInfoFilt(:).Ext_length).*.216;
            % convert NaN lengths of internal to zero
            lengthsInt(isnan(lengthsInt))=0;
            lengthsInt(filterInt) = 0; 
            
            lengths = lengthsInt + lengthsExt;
        end
        
        filoLengthToVeil{iFrame} = lengths;
    else
        filoLengthToVeil{iFrame} = [];
    end
    clear lengths
end
% if ~isempty(outDir)
% save([outDir filesep 'param_filoLengthToVeil.mat'],'filoLengthToVeil');
% end
