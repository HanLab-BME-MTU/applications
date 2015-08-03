function [ sortedList ,sortednum] = sortUnpaddedList(list)



%   Detailed explanation goes here
   %Initialize Cells for Sorting
    %path = path of file, body = body of filename, 
    %(tif etc) (all of these require a cell because they are strings)
    % num = number of filename (do not want in cell so can sort)
    
    nImTot = length(list); 
    
    pathCell = cell(nImTot,1);
    bodyCell = cell(nImTot,1);
    extCell = cell(nImTot,1);
    num = zeros(nImTot,1);
    
        %Sort List
        % For each frame get the image name from listOfImages
        for iFrame =  1:nImTot;
            name= [char(list(iFrame,2)) filesep char(list(iFrame,1))];

            %Call "getFilenameBody" from common dir to split filename into path,
            %body, no, and ext. Put path, body and ext into appropriate cell/number vector for that
            %frame
            [path body no ext ] = getFilenameBody(name);


            pathCell(iFrame) = cellstr(path);
            bodyCell(iFrame) = cellstr(body);
            extCell(iFrame) = cellstr(ext);

            % output "no" is a string so convert to number to sort
            num(iFrame)  = str2double(no);

        end
        %Sort number vector numerically
        sortednum = sort(num);

        %Convert number vector to cell
        sortednum_cell = num2cell(sortednum);

        %Create Sorted Image List
        sortedList = [pathCell, bodyCell, sortednum_cell, extCell];

    end

