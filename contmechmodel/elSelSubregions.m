function roi = elSelSubregions(img,varargin)
% elSelSubregions: Interactively select a set of subregions in an image.
%
% SYNOPSIS:
%    roi = elSelSubregions(img);
%    elSelSubregions(img,'outfile',filename)
%
% INPUT:
%    img : The image in which the subregions are selected.
%
%    OPTIONAL PAR/VALUE:
%      PAR               VALUE
%    'outfile' : A string that specifies the output file name where 
%                the selected subregions are saved.
%
% OUTPUT:
%    roi : A cell array of polygons that define the set of subregions. Each 
%          element of the cell array is a two-column vector [x y] that specifies

outfile = [];

if nargin > 1
   if mod(nargin-1,2) ~= 0
      error('Incorrect number of optional par/value pairs.');
   end

   for k = 1:2:nargin-1
      switch varargin{k}
         case 'outfile'
            outfile = varargin{k+1};
      end
   end
end

numSubRegions = input('How many subregions?\nYour answer: ');
fprintf(1,'\n');

if ~isempty(img)
   figure; imshow(img,[]);
end

figure(gcf); hold on;

roi = cell(1,numSubRegions);
for k = 1:numSubRegions
   title(sprintf('Please crop region No. %d',k));

   [bw,x,y] = roipoly;
   roi{k} = [x y];
end

title(sprintf('A total of %d regions are selected. Done.',numSubRegions));

if ~isempty(outfile)
   save(outfile,'roi');
end
