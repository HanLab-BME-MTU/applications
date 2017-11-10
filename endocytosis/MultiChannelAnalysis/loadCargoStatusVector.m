function data = loadCargoStatusVector(data)

nMovies = length(data);
for i=1:nMovies
    if exist([data(i).source 'TrackInfoMatrices' filesep 'cargoStatus.mat']) == 2
    load([data(i).source 'TrackInfoMatrices' filesep 'cargoStatus.mat']);   
    data(i).cargoStatus = cargoStatus;
    data(i).csValid = valid;
    else
        data(i).cargoStatus = [];
        data(i).csValid = [];
        disp(['movie number ' num2str(i) ' is empty'])
    end
end;