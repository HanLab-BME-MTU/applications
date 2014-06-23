function HelloWorldXT(aImarisApplicationID)
% Hello World! XT example
%
% <CustomTools>
%  <Menu>
%   <Submenu name="ImarisXT course">
%    <Item name="Hello World!" icon="Matlab">
%   <Command>MatlabXT::HelloWorldXT(%i)</Command>
%  </Item>
% </Submenu>
% </Menu>
% </CustomTools>

% Instantiate an IceImarisConnector object
conn = IceImarisConnector(aImarisApplicationID);

% Display version info in a dialog
uiwait(msgbox(['... from IceImarisConnector ', conn.version(), ...
    ' and ', char(conn.mImarisApplication.GetVersion())], ...
    'Hello World!','modal'));

