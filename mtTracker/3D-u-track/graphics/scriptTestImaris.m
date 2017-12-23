%MD=MovieData.load('C:\Users\S158741\data\raw\EB3-KT\cell1_12_half_volume_double_time\deskew\MovieData.mat');
%MDCrop=crop3D(MD,ones(size(MD.getChannel(1).loadStack(1))),'keepFrame',1:10)


pack=MDCrop.searchPackageName('testImaris','selectIdx','last');
testImaris(MDCrop,'package',pack,'iceConnector',origIceConn,'forceRunIdx',2);
