function receptorTraj = compTrack2ReceptorTraj(receptorInfoLabeled,sample,track)

%COMPTRACK2RECEPTORTRAJ extracts the receptor trajectors present in a compound track
%
%SYNPOSIS receptors = compTrack2Receptors(receptorInfoLabeled,sample,track)
%
%INPUT  receptorInfoLabeled: Output of genReceptorInfoLabeled
%       sample             : Number designating which subsampling within
%                            the receptorInfoLabeled structure is to be
%                            used.
%       track              : Number designating which compound track is to
%                            be used
%
%OUTPUT receptors          : A matrix containing the trajectories only of
%                            those receptors present within the compound
%                            track
%
%Code started by Paul Blazek June 2015.

%% Main Body

%extract information from input
receptTraj = receptorInfoLabeled(sample).receptorTraj;
tracksFeatIndxCG = receptorInfoLabeled(sample).compTracks(track).tracksFeatIndxCG;
clust2recept = receptorInfoLabeled(sample).clust2receptAssign;
[~,~,numFrames] = size(receptTraj);

%convert sparse to full if necessary
if issparse(tracksFeatIndxCG)
    tracksFeatIndxCG = full(tracksFeatIndxCG);
end

%extract receptors for the given compound track, looking at all time frames
clustInCompTrack = tracksFeatIndxCG(tracksFeatIndxCG(:,1)~=0,1);
receptors = unique(clust2recept(clustInCompTrack,:,1));
for iFrame = 2 : numFrames
    clustInCompTrack = tracksFeatIndxCG(tracksFeatIndxCG(:,iFrame)~=0,iFrame);
    receptToAdd = clust2recept(clustInCompTrack,:,iFrame);
    receptors = unique([receptors; receptToAdd(:)]);
end
receptors = unique(receptors(:));
receptors = receptors(receptors~=0);

%get trajectories of these receptors
receptorTraj = receptTraj(receptors,:,:);