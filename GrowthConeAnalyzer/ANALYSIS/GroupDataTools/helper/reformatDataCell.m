function [ dataMatPadded ] = reformatDataCell(inputCell)
% reformatDataCell
% Small Helper Function to Convert a Cell with r elements to a matrix with
% rxc matrix 
%%
% INPUT: 
%       inputCell:         any rx1 array cell (such as the output of any of 
%                          the GCAAnalysisParamInTime functions) where r (row) is
%                          typically the number of distributions to reformat (either per
%                          frame or per movie) 
% 
%OUTPUT:     
%       dataMatPadded:     rxc  matrix such that r is the max number number of
%                            observations in the inputCell and c (columns) is  the
%                            number of elements in the original cell 
%                            
%                            
%% 
% collect the max number of observations for all the distributions in the
% cell
toPad = max(cellfun(@(x) length(x),inputCell));

% pad each distribution with the appropriate N number to make a square
% matrix - output here is the same cell but padded with NaNs so all the same 
% length
forDataMat = cellfun(@(x) [x;nan(toPad-length(x),1)],inputCell,'uniformoutput',0); % this will be a dataMat rows = observations and columns = projects in group

% now can catenate the distribution data (each column a distribution) 
% into a data matrix as all rows are the same size. 
dataMatPadded = horzcat(forDataMat{:});


end

