function [] = extractIntensityMat(data, dvector, tbuffer, channel, shiftvar, usemask, filenameRoot)
% read intensities at specified positions (e.g. at the locations of
% detected and tracked objects) from image series and save to a parameter file
% 
% SYNOPSIS extractIntensityMat(data, dvector, tbuffer, channel, shiftvar, usemask, filenameRoot);  
%
% INPUT:    data:       data structure, which contains the location of the
%                       movie data in the .source field
%           dvector:    distance vector for image readout (i.e. the
%                       intensity in the image is read out up to the
%                       distance dvector from the specified object
%                       locations)
%           tbuffer:    time buffer for reading intensity before appearance
%                       and after disappearance (i.e. intensity is read
%                       e.g. 10 frames before and 10 frames after the
%                       actual detected trajectory)
%           channel (optional): if the data structure contains a field that
%                       contains the location of the directory with the 
%                       image data from which the intensities are read, 
%                       specify that field here; the function will then
%                       automatically point you to this directory and save
%                       you lots of unnecessary mouse clicking
%                       DEFAULT is 'source'
%           shiftvar (optional): indicates whether or not the positions
%                       should be shifted by the vector value indicated in
%                       the field .colorShiftVector; 1=yes, 0=no; 
%                       DEFAULT is 0 (no)
%                       NOTE: This functionality is relevant if you are
%                       reading intensities from a dual-channel movie,
%                       where the objects are detected in one channel and
%                       you are reading intensities or parameters from a
%                       second channel that is shifted by a few pixels; see
%                       note below
%           usemask (optional): indicate whether or not a cell outline mask
%                       should be used for the reference points, provided
%                       it's available; 1=yes, 0=no
%                       DEFAULT is 0 (no)
%           filenameRoot (optional): root of the filename for results files
%                       that are written into your source directory, e.g.
%                       'ClathrinInt' will cause a file 'ClathrinInt.mat'
%                       and a file 'ClathrinInt_Ref.mat' to be written as
%                       results files
%                       DEFAULT is ''parameterMat', which will produce the 
%                       files parameterMat.mat' and 'parameterMat_Ref.mat'
    
%
% OUTPUT:   file with results is written into source directory; format of
%           the results file is a matrix with n rows and k columns, where n
%           is the number of objects/trajectories and k is the number of
%           frames and the matrix contains the intensities read from the
%           specified images at the location of the objects in the matrix
%           NOTE: this matrix has the same size as the lifetime
%           matrix, and has corresponding object rows
%
% last modified: Dinah Loerke, March 17, 2009
%
% NOTE regarding color Shift:
% If the data originate from different channels that are shifted, e.g. if
% you have a red and a green color channel that are shifted significantly 
% (several pixels) due to chromatic aberration in the optical system, then 
% your readout positions in the second channel should ideally be corrected
% by that shift. For this purpose, the shift information should be stored 
% in the data structure in a field called .colorShiftVector, in form of a 
% vector, e.g. [2,2] or [-10,5].
% What do the dimensions of the shiftvector mean? In the current
% implementation, a shift vector =[-10,-5] means that image2/channel2
% is shifted in such a way with respect to image1/channel1 that the
% point (1,1) in image1 (e.g. clathrin) is located at point (11,6) in 
% image2 (e.g. actin). Using this convention, we obtain the correct
% coordinates in the second channel by subtracting the shift vector (which
% amounts to adding if the shift is negative) from the positions in the
% first channel.
%
% Last modified: Francois Aguet, Feb 2010

if nargin<5 || isempty(shiftvar)
    shiftvar=0;
end
if nargin<6 || isempty(usemask)
    usemask = 0;
end
if nargin>6 && ~isempty(filenameRoot)
    defName = [filenameRoot,'.mat'];
    defName_ref = [filenameRoot,'_Ref.mat'];
else
    defName = 'parameterMat.mat';
    defName_ref = 'parameterMat_Ref.mat';
end

nMovies = length(data);
ImageStackList(1:nMovies) = struct('list', []);
for i = 1:nMovies
    
    % select first image for the intensity channel
    if nargin>3 && ~isempty(channel)
        ipath = data(i).channels{channel};
    else
        ipath = data(i).source;
    end
    imageName = dir([ipath filesep '*.tif']);
    ImageStackList(i).list = getFileStackNames([ipath imageName(1).name]);
end

% read intensities at designated positions (the positions of detected and
% tracked objects stored in the lifetimematrix)
for i=1:length(data)
    
    % print processing update 
    fprintf('movie #%02d',i);
    fprintf('\n');
    
    % determine tracking data path, and make mpm of all positions in the
    % lifetime matrix (using a number of buffer frames before and after the
    % end of trajectories)
    MPMpos_obj = extractAllPosSideBuffer([data(i).source 'LifetimeInfo'], tbuffer);
    
    % determine reference positions - many intensities or parameters
    % require a reference value within the cell; for this purpose, load the 
    % reference image, if one exists, and choose random positions
    cpath = [data(i).source 'SubregionsMask'];
    if (exist(cpath, 'dir')==7) && (usemask==1)
        mask = imread([cpath filesep 'mask0001.tif']);
        MPMpos_ref = extractAllPosSideBuffer_ref(MPMpos_obj, mask');
    else
        MPMpos_ref = extractAllPosSideBuffer_ref(MPMpos_obj, data(i).imagesize);
    end
    
    % put object and reference positions together, and make time vectors
    sy = size(MPMpos_obj,2);
    MPMall_pos = [MPMpos_obj(:,:,1);
                  MPMpos_ref(:,:,1)];
    tMat_t1 = [MPMpos_obj(:,1:2:sy,2);
               MPMpos_ref(:,1:2:sy,2)];
    tMat_t2 = [MPMpos_obj(:,1:2:sy,3);
               MPMpos_ref(:,1:2:sy,3)];
    
    % if the data originate from different channels that are shifted, then
    % that shift information for this particular movie is stored in a
    % field called .colorShiftVector
    % If such a field exists, then the object positions in MPMglobal have 
    % to be shifted by the shift vector to read the appropriate position in
    % the intensity/parameter images. In the current
    % implementation/convention for dimensionality, we obtain the correct
    % coordinates in the second channel by subtracting the shift vector (or
    % adding if the shift is negative) from the positions in the first 
    % channel. After the subtraction, positions lying outside the shifted 
    % image dimensions have to be set to nan again!
    
    if isfield(data,'colorShiftVector') && (shiftvar==1) && ~isempty(data(i).colorShiftVector)
        msy = size(MPMall_pos,2);
        shiftx = data(i).colorShiftVector(1);
        shifty = data(i).colorShiftVector(2);
        MPMshiftx = MPMall_pos(:,1:2:msy)-shifty;
        MPMshifty = MPMall_pos(:,2:2:msy)-shiftx;
        
        badpos = find( (MPMshiftx>data(i).imagesize(2)) | (MPMshiftx<1) | (MPMshifty>data(i).imagesize(1)) | (MPMshifty<1));
        MPMshiftx(badpos) = nan;
        MPMshifty(badpos) = nan;
        
        MPMall_pos(:,1:2:msy,1) = MPMshiftx;
        MPMall_pos(:,2:2:msy,1) = MPMshifty;
    end
    
    total_frame_num = length(ImageStackList(i).list);
    
%     for k=1:total_frame_num
%         
%         cimage = imread(imageStackList{k});
%         if k==1
%             intensityImageStack = zeros(size(cimage,1),size(cimage,2),total_frame_num);
%         end
%         intensityImageStack(:,:,k) = cimage;
%         
%     end
    
    nf = min(total_frame_num,size(MPMall_pos,2)/2);
    
    % read intensity images matrix
    iMat = extractIntensity_MPMfromImageStack(MPMall_pos, ImageStackList(i).list, dvector);
    
    % separate the imat results back into the objects and the reference
    % points
    bpoint = size(MPMpos_obj,1);
    epoint = size(MPMall_pos,1);
    
    iMat_obj = iMat(1:bpoint,1:nf,:);
    iMat_obj(:,:,2) = tMat_t1(1:bpoint,1:nf);
    iMat_obj(:,:,3) = tMat_t2(1:bpoint,1:nf);
    
    iMat_ref = iMat(bpoint+1:epoint,1:nf,:);
    iMat_ref(:,:,2) = tMat_t1(bpoint+1:epoint,1:nf);
    iMat_ref(:,:,3) = tMat_t2(bpoint+1:epoint,1:nf);    
    
    % save results into the .source directory under the specified names
    save([data(i).source defName], 'iMat_obj');
    save([data(i).source defName_ref],'iMat_ref'); 
end