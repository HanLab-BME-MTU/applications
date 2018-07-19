function imarisDataType = getImarisDataType(matlabDataType)
%GETIMARISDATATYPE returns the imaris data type corresponding to the input matlab data type (class)
%
% imarisDataType = getImarisDataType(matlabDataType);
%
%   Input:
%
%       matlabDataType - Character array specifying Matlab data type (class), e.g. double, single,
%       uint16, uint8 etc
%
%   Output:
%       
%       imarisDataType - Character array specifying Imaris data type, e.g.
%       eTypeFloat, eTypeSingle etc.
%

%Hunter Elliott, 4/2013


switch matlabDataType
    
    case 'single'        
        imarisDataType = 'eTypeFloat';        
    case 'uint8'        
        imarisDataType = 'eTypeUInt8';
    case 'uint16'        
        imarisDataType = 'eTypeUInt16';

    otherwise
        error('GETIMARISDATATYPE:unrecognizedType',['Error! "' matlabDataType '" is not supported!'])
end


