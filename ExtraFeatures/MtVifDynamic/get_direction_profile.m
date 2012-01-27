function track_direction_profiles = get_direction_profile(projData,MD)
nFrame = MD.nFrames_;

SteerableOutputDir = [MD.outputDirectory_,'/VIF_SteerableDetector'];

xCoord = projData.xCoord;
yCoord = projData.yCoord;

nTrack = size(xCoord,1);
track_direction_profiles = cell(1,nTrack);

for iTrack  = 1 :nTrack
    
    x_tr = xCoord(iTrack, :);
    y_tr = yCoord(iTrack, :);
    
    [direction_array,angle_array,path, time_setps, profile] = track_angles_path(x_tr, y_tr, SteerableOutputDir); 
    
    track_direction_profiles{1,iTrack} = profile;
    
end
 