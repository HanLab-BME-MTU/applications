classdef Imaris < handle
    
    % ---------------------
    % This class allows you to interact conveniently with Imaris over the 
    % COM interface. The new Java interface is not yet supported!
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'private',SetAccess = 'private')
        imarisApp;
        cameraOrientation;
        cameraPosition;
        cameraFocus;
        cameraHeight;
        snapshotsCounter = 0;
    end
    
    properties (GetAccess = 'public',SetAccess = 'public')
        displayEnabled = true;
        snapshotsEnabled = false;
        snapshotsPath;
    end
    
    methods
        function obj = Imaris(varargin)
            if numel(varargin) > 0
                obj.displayEnabled = varargin{1};
                if numel(varargin) > 1
                    obj.snapshotsEnabled = varargin{2};
                    obj.snapshotsPath = varargin{3};
                end
            end
            
            if obj.displayEnabled
                % Launch application
                obj.imarisApp = actxserver('Imaris.Application');
                obj.imarisApp.mVisible = 1;
                
                % Create the surpass scene in Imaris
                surpassScene = obj.imarisApp.mFactory.CreateDataContainer;
                surpassScene.mName = 'Imaris: Imaris Scene';
                obj.imarisApp.mSurpassScene = surpassScene;
                
                % Create some dummy data to be able to use the "Save" option of Imaris
                type = 'eTypeUInt8';
                sizeX = 1000;
                sizeY = 1000;
                sizeZ = 1000;
                sizeC = 1;
                sizeT = 1;
                obj.imarisApp.mDataSet.Create(type,sizeX,sizeY,sizeZ,sizeC,sizeT);
            end
        end
        
        clearScene(obj);
        displayCurve(obj,points,name);
        displayDashedCurve(obj);
        displayPoints(obj,points,spotSize,color,name);
        displaySegements(obj);
        displayVectors(obj);
        fitAndSaveCamera(obj);
        fitCamera(obj);
        loadCamera(obj);
        resetScene(obj);
        saveCamera(obj);
        setupScene(obj);
        takeSnapshot(obj);
        takeSnapshotAndResetScene(obj);
        keepWindowOpen(obj);
        
    end % methods
    
end % classdef


