function Display3DRGBMaskOverlayInImaris( imInput, imRGBMask, voxelSpacing )

    imarisApp = imarisStartNew(nargout==0);
    imarisScene = imarisApp.mFactory.CreateDataContainer;
    imarisApp.mSurpassScene = imarisScene;
    imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
    imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene

    volSize = size( imRGBMask );
    volData = imarisApp.mFactory.CreateDataSet;
    volData.Create( 'eTypeUint16', volSize(1), volSize(2),  volSize(3), 4, 1);
    volData.mExtendMinX = 0;
    volData.mExtendMinY = 0;
    volData.mExtendMinZ = 0;
    volData.mExtendMaxX = volData.mSizeX * voxelSpacing(1);
    volData.mExtendMaxY = volData.mSizeY * voxelSpacing(2);
    volData.mExtendMaxZ = volData.mSizeZ * voxelSpacing(3);
    volData.mUnit = 'um';
    
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    volData.SetDataVolume( uint16( imInput ), 0, 0);
    volData.SetChannelColor( 0, 1, 1, 1, 0 );
    volData.SetChannelRange( 0, ImageIntensityRange(1), ImageIntensityRange(2) );
    
    channelColor = [1, 0, 0; 0 1 0; 0 0 1];
    imRGBMask = uint16( mat2gray( imRGBMask ) * 255.0 );
    for i = 1:3
        volData.SetDataVolume( imRGBMask(:,:,:,i), i, 0);
        volData.SetChannelColor( i, channelColor(i,1), channelColor(i,2), channelColor(i,3), 0 );
        volData.SetChannelRange( i, 0, 255 );
    end
    
    imarisApp.mDataSet = volData;
    imarisApp.mSurpassCamera.Fit;

    %Strangely enough, this is both necessary and sufficient to display the mDataSet in the newly created scene.
    imarisScene.AddChild(imarisApp.mFactory.CreateVolume);
   
end