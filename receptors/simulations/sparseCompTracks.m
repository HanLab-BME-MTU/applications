function tracksOut = sparseCompTracks(tracksIn)
%SPARSECOMPTRACKS converts compound tracks between sparse and non-sparse format
%
%SYNOPSIS tracksOut = sparseCompTracks(tracksIn)
%
%INPUT  tracksIn     : Tracks, in form of output of trackCloseGapsKalman. 
%                      They can be sparse or non-sparse.
%
%OUTPUT tracksOut    : Same as input tracks, just sparse if input is
%                      non-sparse, and non-sparse if input is sparse.
%
%Khuloud Jaqaman, December 2014

%% Input

if nargin < 1
    error('sparseCompTracks: Incorrect number of input arguments!')
end

%check if tracks are in sparse form
sparseForm = issparse(tracksIn(1).tracksFeatIndxCG);

%get number of tracks
numTracks = length(tracksIn);

%% Conversion

tracksOut = tracksIn;

if sparseForm %if input is sparse
    
    for iTrack = 1 : numTracks
        
        %convert from sparse to full
        tracksFeatIndxCG = full(tracksIn(iTrack).tracksFeatIndxCG);
        tracksCoordAmpCG = full(tracksIn(iTrack).tracksCoordAmpCG);
        
        %convert zeros to NaNs where necessary
        for iRow = 1 : size(tracksCoordAmpCG,1)
            colZero = find(tracksCoordAmpCG(iRow,:)==0);
            colZero = colZero(:)';
            xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
            colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
            tracksCoordAmpCG(iRow,colZero) = NaN;
        end
        
        %write in output structure
        tracksOut(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
        tracksOut(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
        
    end
    
else %if input is full
    
     for iTrack = 1 : numTracks
        
        %convert NaNs to zeros
        tracksCoordAmpCG = tracksIn(iTrack).tracksCoordAmpCG;
        tracksCoordAmpCG(isnan(tracksCoordAmpCG)) = 0;
        
        %convert from full to sparse
        tracksFeatIndxCG = sparse(tracksIn(iTrack).tracksFeatIndxCG);
        tracksCoordAmpCG = sparse(tracksCoordAmpCG);
        
        %write in output structure
        tracksOut(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
        tracksOut(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
        
    end
    
end


%% ~~~ the end ~~~

