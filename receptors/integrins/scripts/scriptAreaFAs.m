
areaVec = [];

for i = 1 : length(movieStructAlphaVPax);
    
    disp(num2str(i))
    
    if movieStructAlphaVPax(i).activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName{1};
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisFAs\masks'])
        
        dirName = [topDir '\analysisFAs\masks\'];
        fName = ls('*0001.tif');
        
        %get all file names in stack
        outFileList = getFileStackNames([dirName,fName]);
        numFrames = length(outFileList);
        
        %read images
        currentImage = imread(outFileList{1});
        [isx,isy] = size(currentImage);
        imageStack = zeros(isx,isy,numFrames);
        imageStack(:,:,1) = currentImage;
        for iFrame = 2 : numFrames
            imageStack(:,:,iFrame) = imread(outFileList{iFrame});
        end
        imageStack = logical(imageStack);
        
        for iFrame = 1 : numFrames
            stats = regionprops(imageStack(:,:,iFrame),'Area');
            tmp = vertcat(stats.Area);
            areaVec = [areaVec; tmp];
        end
        
    end
    
end
