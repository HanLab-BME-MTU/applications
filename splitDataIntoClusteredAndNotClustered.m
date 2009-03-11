function [lftInfo] = splitDataIntoClusteredAndNotClustered(clusterResults,lftInfo,condition);

%GET DATA FROM lftInfo
%status (complete, good gap, bad gap, cutoff)matrix
statMat = lftInfo.Mat_status;
% lifetime matrix
lftMat = lftInfo.Mat_lifetime;
% x-coordinate matrix
matX = lftInfo.Mat_xcoord;
% y-coordinate matrix
matY = lftInfo.Mat_ycoord;
% disapp status matrix
daMat = lftInfo.Mat_disapp;


%FIND CLUSTERED/UNCLUSTERED PITS IN lftInfo DATA
pitID = zeros(size(clusterResults,1),1);
for ipit = 1:size(clusterResults,1)

    if strcmp(condition,'clustered')
        if clusterResults(ipit,3) ~= 0
            [pitID(ipit),dummy] = find(statMat == 1 & daMat == 1 & lftMat == clusterResults(ipit,4) &...
                matX == clusterResults(ipit,1) & matY == clusterResults(ipit,2));
        elseif clusterResults(ipit,3) == 0
            pitID(ipit) == 0;
        end %
    elseif strcmp(condition,'unclustered')
        if clusterResults(ipit,3) == 0
            [pitID(ipit),dummy] = find(statMat == 1 & daMat == 1 & lftMat == clusterResults(ipit,4) &...
                matX == clusterResults(ipit,1) & matY == clusterResults(ipit,2));
        elseif clusterResults(ipit,3) ~= 0
            pitID(ipit) == 0;
        end %inside/outside
    end %clustered/unclustered
end%for loop
pitID = nonzeros(pitID);

%CONVERT lftInfo SO THAT ONLY CLUSTERED/UNCLUSTERED PITS ARE INCLUDED
lftInfo.Mat_status = lftInfo.Mat_status(pitID,:);
lftInfo.Mat_lifetime = lftInfo.Mat_lifetime(pitID,:);
lftInfo.Mat_xcoord = lftInfo.Mat_xcoord(pitID,:);
lftInfo.Mat_ycoord = lftInfo.Mat_ycoord(pitID,:);
lftInfo.Mat_disapp = lftInfo.Mat_disapp(pitID,:);
