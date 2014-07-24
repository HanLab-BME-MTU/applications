function [imFilterResult] = callMatitkFilter( strFilterCode, cellArrParams, imInput1, imInput2, matSeedPoints, pixelSpacing )

    if nargin < 1
        % spit out all the filters currently present in matitk_custom and some help on how to use a filter
        matitk_custom( '?' ); 
        return;
    end

    if nargin < 2        
        % user wants help on a specific filter in matitk -- so spit out 
        % what ever documentation that is written for it inside matitk
        matitk_custom( strFilterCode ); 
        return;
    end
    
    if ~iscell( cellArrParams )
        error( 'ERROR: cellArrParams must be a cell array of filter parameters' );
    end
    
    if ~exist('imInput2', 'var') || isempty( imInput2 )
       imInput2 = []; 
    else
        if ndims( imInput1 ) ~= ndims( imInput2 )
            error( 'ERROR: both input1 and input2 must be of the same dimension' );
        end
    end

    if ~exist( 'matSeedPoints' ) || isempty( matSeedPoints )
        matSeedPoints = [];
    else
        if ~isnumeric( matSeedPoints ) || size( matSeedPoints, 2 ) ~= ndims( imInput1 )
            error( 'ERROR: dimension of the seed points must match the dimension of the input image' );
        end
    end
    
    if ~exist( 'pixelSpacing', 'var') || isempty( pixelSpacing )
        pixelSpacing = ones(1,ndims(imInput1));
    end
    
    
    switch ndims( imInput1 ) 
       
        case 2
            
            % Matitk currently only support 3D images as input (example of 
            % how stupid a toolbox can be). So convertinputs from 2D to 3D 
            % by replicating image into two slices. 
            imInput1 = cat( 3, imInput1, imInput1 );          
            pixelSpacing(3) = 1;
            
            if ~isempty( imInput2 )
                imInput2 = cat( 3, imInput2, imInput2 );          
            end
            
            if ~isempty(matSeedPoints)
                matSeedPoints(:,3) = 1;
            end            
            
            imFilterResult = matitk_custom( strFilterCode, cellArrParams, imInput1, imInput2, matSeedPoints, pixelSpacing );
            imFilterResult = imFilterResult(:,:,1);
            
        case 3
        
            imFilterResult = matitk_custom( strFilterCode, cellArrParams, imInput1, imInput2, matSeedPoints, pixelSpacing );            
            
        otherwise
            
            error( 'ERROR: matitk currently does not support images of this dimension' );
        
    end
    
end