function tracksNewRef=addReferential(tracks,refArray,name)
% This function create duplicate of the master tracks in different
% referential and index them on the  master file.

%trackNewRef(length(refArray))=tracks(1);
tracksNewRef=cell(1,length(refArray));
for rIdx=1:length(refArray)
    tracksNewRef{rIdx}=refArray(rIdx).applyBaseToTrack(tracks);
end

for tIdx=1:length(tracks)
    try
    tracks(tIdx).addprop(name);
    catch 
    end
    setfield(tracks(tIdx),name,cellfun(@(t) t(tIdx),tracksNewRef,'unif',0));
end