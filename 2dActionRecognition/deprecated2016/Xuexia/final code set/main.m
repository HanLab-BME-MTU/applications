function [vecBoundingBoxSize,vecFragmentVideoNumLabel,vecVideoClassLabel] = main(SELECTED_TIFF_XLSX_FILENAME,SELECTED_TIFF_LABELS_XLSX_FILENAME,TIME_WINDOW,MOTION_THRESHOLD)
%main script

[~,cellFilename] = xlsread(SELECTED_TIFF_XLSX_FILENAME);
vecVideoClassLabel = xlsread(SELECTED_TIFF_LABELS_XLSX_FILENAME);

%initialize storage matricies
vecBoundingBoxSize = [];
vecFragmentVideoNumLabel = [];
matPosHist = [];
matNegHist = [];

%counters for video and fragment
labelCounter = 0;
fragmentCounter = 0;
%loop through all files
for i = 1:length(cellFilename)
    %updates video counter
    labelCounter = labelCounter + 1;
    
    temp_filename = cellFilename{i};
    
    %splits the video into multiple fragments
    cellFragments = splitVideo(temp_filename,TIME_WINDOW);
    %loops through all fragments in current video
    for j = 1:length(cellFragments)
        %updates fragment counter
        fragmentCounter = fragmentCounter + 1;
        
        temp_fragment = cellFragments{j};
        
        %calculates the MIP transform (dim MxNx8)
        [temp_matPosMIP,temp_matNegMIP] = getMIP(temp_fragment,MOTION_THRESHOLD);
        %calculates the histogram feature vectors
        [temp_matPosHist,temp_matNegHist] = makeDescriptorMIP(temp_matPosMIP,temp_matNegMIP);
        
        %store values
        matPosHist = cat(length(size(temp_matPosMIP))+1,matPosHist,temp_matPosHist);
        matNegHist = cat(length(size(temp_matNegMIP))+1,matNegHist,temp_matNegHist);
        vecBoundingBoxSize = [vecBoundingBoxSize,size(temp_fragment,1)*size(temp_fragment,2)];
        vecFragmentVideoNumLabel = [vecFragmentVideoNumLabel,labelCounter];
    end 
end
end