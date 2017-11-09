MD=MovieData.load('/project/bioinformatics/Danuser_lab/shared/proudot/3d-vis/utrackPackage/scriptBased/UTrack-QD-v1/data/Cell12Toy2/analysis/movieData.mat');

loadStruct=load([MD.outputDirectory_ filesep 'detection/pointSourceAutoSigmaFit/tracks/tracksHandle.mat']);
tracks=loadStruct.tracks;

%% Track lifetime (edge property)
tracksLft=[tracks.lifetime]';

%% Track Median Speed (edge property)
s=[MD.pixelSize_,MD.pixelSize_,MD.pixelSizeZ_,MD.timeInterval_];
tracksMedSpeed= arrayfun(@(t) nanmedian(sum((   [s(1)*t.x(1:end-1);s(2)*t.y(1:end-1);s(3)*t.z(1:end-1)]- ...
    [s(1)*t.x(2:end);  s(2)*t.y(2:end); s(3)*t.z(2:end)  ]).^2).^0.5/s(4)) ,tracks);

%% Track Max Speed (edge property)
tracksMaxSpeed= arrayfun(@(t)    nanmax(sum((   [s(1)*t.x(1:end-1);s(2)*t.y(1:end-1);s(3)*t.z(1:end-1)]- ...
                                                    [s(1)*t.x(2:end);  s(2)*t.y(2:end); s(3)*t.z(2:end)  ]).^2).^0.5/s(4)) ,tracks,'unif',0);

%% Track diffCoeff (edge property)
tracksDiffCoeff=arrayfun(@(t) nanmean(sum([s(1)*t.x(1)-s(1)*t.x(2:end); s(2)*t.y(1)-s(2)*t.y(2:end); s(3)*t.z(1)-s(3)*t.z(2:end)].^2))/(6*t.lifetime*s(4)) ,tracks);
