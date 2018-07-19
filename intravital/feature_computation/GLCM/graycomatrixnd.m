function [GLCMS,SI] = graycomatrixnd(varargin)
%graycomatrixnd Create gray-level co-occurrence matrix on an N-D image
%   GLCMS = graycomatrixnd(I) analyzes pairs of adjacent pixels in
%   a scaled version of I.  If I is a binary image, it is scaled to 2
%   levels. If I is an intensity image, it is scaled to 8 levels. In this
%   case, there are 8 x 8 = 64 possible ordered combinations of values for
%   each pixel pair. graycomatrixnd accumulates the total occurrence of each
%   such combination, producing a 8-by-8 output array, GLCMS. The row and
%   column subscripts in GLCMS correspond respectively to the first and second
%   (scaled) pixel-pair values.
%
%   GLCMS = graycomatrixnd(I,PARAM1,VALUE1,PARAM2,VALUE2,...) returns one or more
%   gray-level co-occurrence matrices, depending on the values of the optional
%   parameter/value pairs. Parameter names can be abbreviated, and case
%   does not matter.
% 
%   This implementation extends the in-built matlab function graycomatrix
%   to work on ND images and in addition it allows the user can provide an 
%   ROI Mask within which he wants the GLCM matrix to be computed.
% 
%   Parameters include:
%  
%       'Offset'         A numOffsets-by-numImageDimensions array of offset vectors 
%                        specifying the distance between the pixel-of-interest 
%                        and its neighbor. Each row array specifies the 
%                        relationship, or 'Offset', between a pair of pixels.
%                        Note that the first dimension corresponds to the rows.
%                        For example, if you want the number of occurrences 
%                        where the pixel of interest is one pixel to the
%                        left of its neighbor in a 3D image, then the 
%                        corresponding offset vector is [0, 1, 0].
% 
%                        Because this offset is often expressed as an angle, the
%                        following table lists the offset values that specify common
%                        angles for a 2D image, given the pixel distance D.
% 
%                        AngleXY     OFFSET
%                        -------     ------  
%                        0           [0 D]   
%                        45          [-D D]
%                        90          [-D 0]
%                        135         [-D -D]  
% 
%                        Default: [1] for 1D, [0 1] for 2D, [0 1 0] for 3D and so on
% 
%       'NumLevels'      An integer specifying the number of gray levels to use when
%                        scaling the grayscale values in I. For example, if
%                        'NumLevels' is 8, graycomatrixnd scales the values in I so
%                        they are integers between 1 and 8.  The number of gray levels
%                        determines the size of the gray-level co-occurrence matrix
%                        (GLCM).
% 
%                        'NumLevels' must be an integer. 'NumLevels' must be 2 if I
%                        is logical.
% 
%                        Default: 8 for numeric
%                                 2 for logical
% 
%       'GrayLimits'     A two-element vector, [LOW HIGH], that specifies how the
%                        grayscale values in I are linearly scaled into gray
%                        levels. Grayscale values less than or equal to LOW are
%                        scaled to 1. Grayscale values greater than or equal to
%                        HIGH are scaled to HIGH.  If 'GrayLimits' is set to [],
%                        graycomatrixnd uses the minimum and maximum grayscale values
%                        in I as limits, [min(I(:)) max(I(:))].
% 
%                        Default: the LOW and HIGH values specified by the
%                        class, e.g., [LOW HIGH] is [0 1] if I is double and
%                        [-32768 32767] if I is int16.
% 
%       'Symmetric'      A Boolean that creates a GLCM where the ordering of
%                        values in the pixel pairs is not considered. For
%                        example, when calculating the number of times the
%                        value 1 is adjacent to the value 2, graycomatrixnd
%                        counts both 1,2 and 2,1 pairings, if 'Symmetric' is
%                        set to true. When 'Symmetric' is set to false,
%                        graycomatrixnd only counts 1,2 or 2,1, depending on the
%                        value of 'offset'. The GLCM created in this way is
%                        symmetric across its diagonal, and is equivalent to
%                        the GLCM described by Haralick (1973).
% 
%                        The GLCM produced by the following syntax, 
% 
%                        graycomatrixnd(I, 'offset', [0 1], 'Symmetric', true)
% 
%                        is equivalent to the sum of the two GLCMs produced by
%                        these statements.
% 
%                        graycomatrixnd(I, 'offset', [0 1], 'Symmetric', false) 
%                        graycomatrixnd(I, 'offset', [0 -1], 'Symmetric', false) 
% 
%                        Default: false
% 
%       'ROIMask'        This should be a binary image of the same size as I and 
%                        specifies the Region of Interest (ROI) in which to compute 
%                        the GLCM. 
% 
%                        Default:
% 
%                        If this parameter is not provided the GLCM is computed over the
%                        the entire image.
% 
%                        Ex:
% 
%                           graycomatrixnd(I, 'ROIMask', imROIMask )                   
% 
% 
%   'maskOutNeighbors'   True/False
% 
%                        Has any effect only when the ROIMask is specified.
%                          
%                        Specifies whether or not to exclude a pixel-pair 
%                        if the neighboring pixel in the pair falls outside 
%                        the specified ROIMask.
% 
%                        Default: false
%  
%   [GLCMS,SI] = graycomatrixnd(...) returns the scaled image used to
%   calculate GLCM. The values in SI are between 1 and 'NumLevels'.
%
%   Class Support
%   -------------             
%   I can be any numeric or logical ND matrix.  I must be real, and nonsparse. 
%   SI is a double matrix having the same size as I.  GLCMS is an
%   'NumLevels'-by-'NumLevels'-by-P double array where P is the number of
%   offsets in OFFSET.
%  
%   Notes
%   -----
%   Another name for a gray-level co-occurrence matrix is a gray-level
%   spatial dependence matrix.
%
%   graycomatrixnd ignores pixels pairs if either of their values is NaN. It also
%   replaces Inf with the value 'NumLevels' and -Inf with the value 1.
%
%   graycomatrixnd ignores border pixels, if the corresponding neighbors
%   defined by 'Offset' fall outside the image boundaries.
%
%   References
%   ----------
%   Haralick, R.M., K. Shanmugan, and I. Dinstein, "Textural Features for
%   Image Classification", IEEE Transactions on Systems, Man, and
%   Cybernetics, Vol. SMC-3, 1973, pp. 610-621.
%
%   Haralick, R.M., and L.G. Shapiro. Computer and Robot Vision: Vol. 1,
%   Addison-Wesley, 1992, p. 459.
%
%   Example 1
%   ---------      
%   Calculate the gray-level co-occurrence matrix (GLCM) and return
%   the scaled version of the image, SI, used by graycomatrixnd to generate the
%   GLCM.
%
%        I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%       [GLCMS,SI] = graycomatrixnd(I,'NumLevels',9,'G',[])
%     
%   Example 2
%   ---------  
%   Calculate the gray-level co-occurrence matrix for a grayscale image.
%  
%       I = imread('circuit.tif');
%       GLCMS = graycomatrixnd(I,'Offset',[2 0])
%
%   Example 3
%   ---------
%   Calculate gray-level co-occurrences matrices for a grayscale image
%   using four different offsets.
%  
%       I = imread('cell.tif');
%       offsets = [0 1;-1 1;-1 0;-1 -1];
%       [GLCMS,SI] = graycomatrixnd(I,'Of',offsets); 
%
%   Example 4
%   ---------  
%   Calculate the symmetric gray-level co-occurrence matrix (the Haralick
%   definition) for a grayscale image.
%  
%       I = imread('circuit.tif');
%       GLCMS = graycomatrixnd(I,'Offset',[2 0],'Symmetric', true)
%  
%   See also: Compute_GLCM_Properties, graycoprops
% 
%   Author: Deepak Roy Chittajallu (Mar 11, 2013)
%

[I, Offset, NL, GL, makeSymmetric, imROIMask, flagMaskOutNeighbors ] = ParseInputs(varargin{:});

% get image dimensions
numDims = ndimsext(I);

% Scale I so that it contains integers between 1 and NL.
if GL(2) == GL(1)
  SI = ones(size(I));
else
  slope = (NL-1) / (GL(2) - GL(1));
  intercept = 1 - (slope*(GL(1)));
  SI = round(imlincomb(slope,I,intercept,'double'));
end

% Clip values if user had a value that is outside of the range, e.g., double
% image = [0 .5 2;0 1 1]; 2 is outside of [0,1]. The order of the following
% lines matters in the event that NL = 0.
SI(SI > NL) = NL;
SI(SI < 1) = 1;

numOffsets = size(Offset,1);

if NL ~= 0

  % Create vectors of row and column subscripts for every pixel and its neighbor.
  s = size(I);  
  ndgridargs = cell(1,numDims);
  for i = 1:numDims
    ndgridargs{i} = 1:s(i);
  end
  
  pixSubInd = cell(1,numDims);
  [pixSubInd{:}] = ndgrid(ndgridargs{:});
  
  for i = 1:numDims
      pixSubInd{i} = pixSubInd{i}(imROIMask > 0);
  end
  
  % Compute GLCMS
  GLCMS = zeros(NL,NL,numOffsets);
  
  for k = 1 : numOffsets
      
    GLCMS(:,:,k) = computeOneGLCM(pixSubInd, Offset(k,:), SI, NL, ...
                                  imROIMask, flagMaskOutNeighbors);
    
    if makeSymmetric 
        % Reflect glcm across the diagonal
        glcmTranspose = GLCMS(:,:,k).';
        GLCMS(:,:,k) = GLCMS(:,:,k) + glcmTranspose;
    end
    
  end

else
  GLCMS = zeros(0,0,numOffsets);
end

%-----------------------------------------------------------------------------
function oneGLCM = computeOneGLCM(pixSubInd, offset, si, nl, imROIMask, flagMaskOutNeighbors)
% computes GLCM given one Offset

  numDims = ndimsext(si);
  
  % Compute subindices of neighbors by adding the offset to subindices of
  % the current pixel
  neighPixSubInd = cell(1, numDims);
  for i = 1:numDims
      neighPixSubInd{i} = pixSubInd{i} + offset(i);
  end
  
  % Remove pixel and its neighbor if they have subscripts outside the image
  % boundary 
  s = size(si);
  flagBadPixel = false(size(pixSubInd,1), 1);
  
  for i = 1:numDims
    badind = find( pixSubInd{i} < 1 | neighPixSubInd{i} < 1 | ...
                   pixSubInd{i} > s(i) | neighPixSubInd{i} > s(i) );
    flagBadPixel(badind) = true;
  end
  
  for i = 1:numDims
       pixSubInd{i}(flagBadPixel) = [];
       neighPixSubInd{i}(flagBadPixel) = [];
  end
  
  % Create vectors containing the intensity values of each component of 
  % the pixel pair.
  pixLinInd = sub2ind(s, pixSubInd{:});
  neighPixLinInd = sub2ind(s, neighPixSubInd{:});
  
  v1 = si( pixLinInd );
  v2 = si( neighPixLinInd );
  
  % Make sure that v1 and v2 are column vectors.
  v1 = v1(:);
  v2 = v2(:);
  
  % If requested - remove pixel and its neighbor if any of them are outside
  % the roi mask
  if flagMaskOutNeighbors
    v1(~imROIMask(:)) = [];
    v2(~imROIMask(:)) = [];
  end
  
  % Remove pixel and its neighbor if their value is NaN.
  bad = isnan(v1) | isnan(v2);
  if any(bad)
    wid = sprintf('Images:%s:scaledImageContainsNan',mfilename);
    msg = 'GLCM does not count pixel pairs if either of their values is NaN.';
    warning(wid,'%s', msg);
  end
  Ind = [v1 v2];
  Ind(bad,:) = [];
  
  if isempty(Ind)
    oneGLCM = zeros(nl);
  else
    % Tabulate the occurrences of pixel pairs having v1 and v2.
    oneGLCM = accumarray(Ind, 1, [nl nl]);
  end
  
%-----------------------------------------------------------------------------
function [numDims] = ndimsext(M)

    numDims = ndims(M);
    if numel(M)==max(size(M))
        numDims = 1;
    end        

%-----------------------------------------------------------------------------
function [I, offset, nl, gl, sym, imROIMask, flagMaskOutNeighbors] = ParseInputs(varargin)

narginchk(1,13);

% Check I
I = varargin{1};
validateattributes(I,{'logical','numeric'},{'real','nonsparse'}, ...
                   mfilename,'I',1);

numDims = ndimsext(I);

% Assign Defaults

    % Offset
    if numDims > 1
        offset = zeros(1,numDims);
        offset(2) = 1;
    else    
        offset = 1;
    end

    % number of quantized gray levels
    if islogical(I)
      nl = 2;
    else
      nl = 8;
    end

    % gray limit
    gl = getrangefromclass(I);
    
    % should the glcm be symmetric
    sym = false;
    
    % roi mask - only pixels within this mask will contribute to the GLCM
    imROIMask = ones(size(I));

    % whether or not to ignore pixel-pairs with neighbor outside ROI mask
    flagMaskOutNeighbors = false;
    
% Parse Input Arguments
if nargin ~= 1
 
  paramStrings = {'Offset', 'NumLevels', 'GrayLimits', ...
                  'Symmetric', 'ROIMask', 'maskOutNeighbors'};
  
  for k = 2:2:nargin

    param = lower(varargin{k});
    inputStr = validatestring(param, paramStrings, mfilename, 'PARAM', k);
    idx = k + 1;  %Advance index to the VALUE portion of the input.
    if idx > nargin
      eid = sprintf('Images:%s:missingParameterValue', mfilename);
      msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
      error(eid,'%s', msg);        
    end
    
    switch (inputStr)
     
        case 'Offset'

            offset = varargin{idx};
            validateattributes(offset,{'logical','numeric'},...
                               {'2d','nonempty','integer','real'},...
                               mfilename, 'OFFSET', idx);

            if (numDims > 1 && size(offset,2) ~= numDims) || ...
               (numDims == 1 && numel(offset) ~= numDims)    

                eid = sprintf('Images:%s:invalidOffsetSize',mfilename);
                msg = 'OFFSET must be an numOffsets x numImageDimensions array.';
                error(eid,'%s',msg);

            end

            offset = double(offset);

        case 'NumLevels'

            nl = varargin{idx};
            validateattributes(nl,{'logical','numeric'},...
                               {'real','integer','nonnegative','nonempty','nonsparse'},...
                               mfilename, 'NL', idx);

            if numel(nl) > 1

                eid = sprintf('Images:%s:invalidNumLevels',mfilename);
                msg = 'NL cannot contain more than one element.';
                error(eid,'%s',msg);

            elseif islogical(I) && nl ~= 2

                eid = sprintf('Images:%s:invalidNumLevelsForBinary',mfilename);
                msg = 'NL must be two for a binary image.';
                error(eid,'%s',msg);

            end

            nl = double(nl);      

        case 'GrayLimits'

            gl = varargin{idx};

            % step 1: checking for classes
            validateattributes(gl, {'logical','numeric'}, ...
                                   {},mfilename, 'GL', idx);
            if isempty(gl)
                gl = [min(I(:)) max(I(:))];
            end

            % step 2: checking for attributes
            validateattributes(gl,{'logical','numeric'}, ...
                                  {'vector','real'}, mfilename, 'GL', idx);

            if numel(gl) ~= 2
                error(message('images:graycomatrix:invalidGrayLimitsSize'));
            end

            gl = double(gl);

        case 'Symmetric'

            sym = varargin{idx};
            validateattributes(sym, {'logical'}, ...
                                    {'scalar'}, mfilename, 'Symmetric', idx);

        case 'ROIMask'

            validateattributes(varargin{idx}, ...
                               {'numeric', 'logical'}, ...
                               {'nonsparse', 'nonempty'}, mfilename, ...
                               'ROIMask', idx );

            % check size
            if ~isempty( find( size(varargin{idx}) ~= size(I) ) )            
                error('%s : Input Image and ROIMask should be of the same dimension', ...
                      mfilename);            
            end

            imROIMask = varargin{idx};        
        
        case 'maskOutNeighbors'
            
            flagMaskOutNeighbors = varargin{idx};
            validateattributes(sym, {'logical'}, ...
                                    {'scalar'}, mfilename, ...
                                    'maskOutNeighbors', idx);
            
        
    end
    
  end
  
end
