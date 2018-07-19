function [ tracksNA ] = changeTrackStateFormat( tracksNA )
%function [ tracksNA ] = changeTrackStateFormat( tracksNA ) changes states
%from character cell format to scalar
%   input:      tracksNA with charcter format status
%   output:      tracksNA with charcter scalar status
% state:
%       1: BA, or Before Adhesion
%       2: NA, Nascent Adhesion
%       3: FC, Focal Complex
%       4: FA, Focal Adhesion
%       5: ANA, After Adhesion
%       6: Out_of_Band
% Sangyoon Han June 108

if iscell(tracksNA(1).state)
    orgTracksNA=tracksNA;
    tracksNA = rmfield(tracksNA,'state');
    for ii=1:numel(tracksNA)
        for jj=tracksNA(ii).iFrame
            if strcmp(orgTracksNA(ii).state(jj),'BA')
                tracksNA(ii).state(jj) = 1;
            elseif strcmp(orgTracksNA(ii).state(jj),'NA')
                tracksNA(ii).state(jj) = 2;
            elseif strcmp(orgTracksNA(ii).state(jj),'FC')
                tracksNA(ii).state(jj) = 3;
            elseif strcmp(orgTracksNA(ii).state(jj),'FA')
                tracksNA(ii).state(jj) = 4;
            elseif strcmp(orgTracksNA(ii).state(jj),'ANA')
                tracksNA(ii).state(jj) = 5;
            elseif strcmp(orgTracksNA(ii).state(jj),'Out_of_Band')
                tracksNA(ii).state(jj) = 6;
            else
                tracksNA(ii).state(jj) = 6;
            end
        end
%         progressText(ii/numel(tracksNA),'changing status');
    end
elseif isscalar(tracksNA(1).state)
    disp('tracksNA.state is already in scalar format. Returning the original tracksNA...')
end
end

