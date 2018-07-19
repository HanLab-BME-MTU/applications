function [out] = GCAAnalysisExtract_filoLength(analInfo,filoFilterSet,varargin)
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
%                      you want to extract (DEFAULT: 'Ext_') - Added 20160215 can also
%                      calculate percentage of bundle embedded by flagging
%                      'EmbedPercent'
%
% OUTPUT:
% filoLengthToVeil:  rx1 cell array where r (row) is the number of frames
%                    and each cell holds the distribution of filo lengths
%                    for the given frame
%
%

%% Check Parameters 

ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('filoPart', 'Ext_'); 
ip.addParameter('outPercent',false); % in the case of filoPart 
% 'Tot' will output the percentage embedded for each filopodia 
% instead of the total length of the actin bundle 
ip.addParameter('filterZeroPercent',true); 

ip.addParameter('umPerPixel',.216);

ip.parse(varargin{:});
%%

% 

if length(analInfo) == 1
    nFrames = 1; 
else 
    nFrames = length(analInfo)-1; 
end 
out = cell(nFrames,1);

for iFrame = 1:nFrames
    
    
    filoInfo = analInfo(iFrame).filoInfo;
    
    if ~isempty(filoInfo);
        % currently 2 columns of the filter - one for the external and one for the 
        % internal based on fitting criteria. 
        filterFrameAll= filoFilterSet{iFrame};
        
        if strcmpi(ip.Results.filoPart,'Int_'); 
            filterFrameC = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) == 1); 
        else 
            filterFrameC = filterFrameAll(:,1); 
       
        end %
        

     
        filoInfoFilt  = filoInfo(filterFrameC);
        
        % collect lengths: if just int or ext just use the respective filopodia
        % filter... 
        if ~strcmpi(ip.Results.filoPart,'Tot');
            lengths =  vertcat(filoInfoFilt(:).([(ip.Results.filoPart) 'length'])).*ip.Results.umPerPixel; % add as a parameter
        else
            % for now filter out all the internal filo that do not pass the
            % fitting criteria ( ie want to keep all external filo that
            % pass the criteria but not the embedded) 
            % 
            
             filterInt = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) ==0 ); % get the ID of all non-fits internally this filter is the length of 
%             % the original filoInfo detection 
             filterInt = filterInt(filterFrameC); % filter the above logical filter to make it the same   
            %filterInt = filterFrameC(:,2) ==0;  % get all the embedded that do not pass the filter criteria
            
            lengthsInt =  vertcat(filoInfoFilt(:).Int_length).*ip.Results.umPerPixel;
            lengthsExt = vertcat(filoInfoFilt(:).Ext_length).*ip.Results.umPerPixel;
            % convert NaN lengths of internal to zero
            lengthsInt(isnan(lengthsInt))=0;
            lengthsInt(filterInt) = 0;
            
            lengths = lengthsInt + lengthsExt;
            if ip.Results.outPercent; % 
                percent = lengthsInt./lengths; 
                if ip.Results.filterZeroPercent 
                    percent = percent(percent~=0); 
                end 
            end 
        end
        if ip.Results.outPercent; 
            out{iFrame} = percent;
        else
            out{iFrame} = lengths;
        end
    else
        out{iFrame} = [];
    end
    clear lengths
end
% if ~isempty(outDir)
% save([outDir filesep 'param_filoLengthToVeil.mat'],'filoLengthToVeil');
% end
