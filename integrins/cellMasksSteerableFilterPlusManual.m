function cellMasksSteerableFilterPlusManual(firstImageFile,firstMaskFile)

%get image information
[fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
dirImage=[fpath,filesep];
fNameImage=[fname,fno,fext];

%get mask information
[fpath,fnameMask]=getFilenameBody(firstMaskFile);
dirMask=[fpath,filesep];

%get all files in image stack
inFileList = getFileStackNames([dirImage,fNameImage]);
numFiles = length(inFileList);

%analyze each image
for iFile = 1 : numFiles
    
    try
    
    %read image
    image = imread(inFileList{iFile});
    
    %generate mask
    genMask = 1;
    while genMask
        mask = steerableFilterPlusManual(image);
        disp('Please close image when ready')
        userEntry = input('Move on to next image? "y" to go to next image, "n" to repeat analysis of this image ','s');
        genMask = strcmp(userEntry,'n');
    end
    
    %save mask
    [fpath,fname,fno,fext]=getFilenameBody(inFileList{iFile});
    filename = [dirMask,fnameMask,fno,fext];
    imwrite(mask,filename);
    
    catch
        
        disp('image ' num2str(iFile) ' failed');
    
    end
    
    
    
end