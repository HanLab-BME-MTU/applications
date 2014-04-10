function [ imageDataStandardized ] = StandardizeIntravitalImageData( imageData )

    imageDataStandardized = cell(size(imageData));
    
    % standardize the intensity ranges
    for i = 1:numel(imageData)
       imageDataStandardized{i} = mat2gray(imageData{i}) * 4096;
    end
    
end