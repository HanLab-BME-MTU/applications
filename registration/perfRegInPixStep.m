function ResT = perfRegInPixStep(mode, inputFileList, target_dir, T_filename_path)
% mode:     Is either 'pixelwise' or 'interpolative' transformation.
%           The default is 'interpolative'
%
% method:   Is either 'cumulative' or 'pairwise' transformation.
%           The default is 'pairwise' transformation


if nargin <1 || isempty(mode) || ~strcmpi(mode,'pixelwise')
    mode='interpolative';
    disp('A subpixel transformation (and interpolation) is performed');
elseif strcmpi(mode,'pixelwise')
    mode='pixelwise';
    disp('A pixelwise transformation (no interpolation) is performed');
end

if nargin <2 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.TIF;*.tif;*.jpg;*.png;*.*'}, ...
       'Select First Image to be registered');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

%get the target directory:
if nargin < 3 || isempty(target_dir)
    target_dir = uigetdir('','Select target directory for registered images');
end

%get the transformation file:
if nargin < 4 || isempty(T_filename_path)
    [T_filename, T_pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select the stage drift correction file');
     T_filename_path=[ T_pathname filesep T_filename];
     fileStruct=load(T_filename_path);
     T=fileStruct.T;
elseif ischar(T_filename_path)
     [T_pathname, ~, ~, ~]=getFilenameBody(T_filename_path);
     fileStruct=load(T_filename_path);
     T=fileStruct.T;
else
     T=T_filename_path;
end

nFiles = numel(inputFileList);
nTrans = length(T);


if nTrans==nFiles-1
    method='cumulative';
    fprintf('A cumulative transformation is performed');
else 
    method='pairwise';
    fprintf('A pairwise transformation is performed');
end


switch lower(mode)
    case 'interpolative'
        switch(method)
            case 'pairwise'
                %This is the default
                intT=T;    
            case 'cumulative'
                intT = cumsum(T);
        end                
    case 'pixelwise'
        switch(method)
            case 'pairwise'
                intT = round(T); 
            case 'cumulative'
                intT = round(cumsum(T));               
        end
end

        
if strcmpi(method,'pairwise') && nTrans<nFiles/2
    fprintf('\nPadded T with ZEROS (no transformation) length did not match number of images!');
    intT(nTrans+1:nFiles/2,:)=0;
    T(nTrans+1:nFiles/2,:)=0;
end

if strcmpi(method,'pairwise')
    I1 = double(imread(inputFileList{1}));
    if length(inputFileList)>=3
        I3 = double(imread(inputFileList{3}));
    else
        I3=I1;
    end
    if sum(sum(I3-I1))==0
        fprintf('\nIt is assumed, that images no: 1,3,5... are the reference frames');
        deltaIndex=2;
    else
        fprintf('\nThis is the wrong order, ref images should be no: 1,3,5');
        fprintf('\nYou should not use the pairwise mode!');
        fprintf('\nNothing has been done!');
        return;
    end
end

if strcmpi(method,'cumulative')
    deltaIndex=1;
end

maxX = ceil(max(abs(intT(:, 2))));
maxY = ceil(max(abs(intT(:, 1))));

I = double(imread(inputFileList{1}));
I = padarray(I, [maxY, maxX]);

[~, name, no] = getFilenameBody(inputFileList{1});
%the first image is the reference frame and remains the same thus save as
%it is:
imwrite(uint16(I),[target_dir, filesep, 'registered_', name, no,'.tif'],'tif','Compression','none');
    
projI = I;
projR = zeros(size(I));

transfIndex=1;
for frameIndex = 2:deltaIndex:nFiles
    % Read image.
    I = double(imread(inputFileList{frameIndex}));
    I = padarray(I, [maxY, maxX]);
    if strcmpi(method,'cumulative') || mod(frameIndex,2)==0
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(intT(transfIndex, :)) 1]);
        R = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);        
    elseif strcmpi(method,'pairwise') && mod(frameIndex,2)==1
        R=I;
    end
    [~, ~, no] = getFilenameBody(inputFileList{frameIndex});
    imwrite(uint16(R),[target_dir, filesep, 'registered_', name, no,'.tif'],'tiff','Compression','none');
    projI = projI + I;
    projR = projR + R;
    transfIndex=transfIndex+1;
end

%The reference frames have to be filled in:
if strcmpi(method,'pairwise')
    for frameIndex = 3:deltaIndex:nFiles
    	I = double(imread(inputFileList{frameIndex}));
        R = padarray(I, [maxY, maxX]);
        [~, ~, no] = getFilenameBody(inputFileList{frameIndex});
        imwrite(uint16(R),[target_dir, filesep, 'registered_', name, no,'.tif'],'tiff','Compression','none');
    end
end

figure(2)
colormap('jet');
subplot(1, 2, 1); imagesc(projI); title('Original stack projection');
subplot(1, 2, 2); imagesc(projR); title('Registered stack projection');

% ResT is the residual transformation which remains to be performed on
% all frames to register them with the FIRST image.

switch lower(mode)
    case 'interpolative'
        fprintf('\nNothing to save since subpixel interpolation has been applied!');          
    case 'pixelwise'
        switch(method)
            case 'pairwise'
                T = T-intT;
                save([T_pathname,filesep,'ResdualT.mat'], 'T');
            case 'cumulative'
                T = cumsum(T)-intT;
                save([T_pathname,filesep,'ResdualCumT.mat'], 'T');              
        end
end


% test if it really worked:
% [~, name, no] = getFilenameBody(inputFileList{1});
% I = double(imread([target_dir,filesep,'registered_', name, no,'.tif']));
% figure(11) 
% imagesc(I)
% 
% projI = I;
% for i = 2:n
%     [~, name, no] = getFilenameBody(inputFileList{i});
%     I = double(imread([target_dir, filesep,'registered_', name, no,'.tif']));    
%     projI = projI + I;
% end
% figure(3)
% imagesc(projI)

%if exist(pathForRegIm,'dir')==0
%    mkdir(pathForRegIm);
%end
