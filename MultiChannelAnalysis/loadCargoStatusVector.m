function data = loadCargoStatusVector(data)

nMovies = length(data);
for i=1:nMovies
    load([data(i).source 'TrackInfoMatrices' filesep 'cargoStatus.mat']);   
    data(i).cargoStatus = cargoStatus;
    data(i).csValid = valid;
end;