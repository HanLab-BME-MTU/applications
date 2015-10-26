classdef OrientationSpaceFilter < handle
    %SteerableVanGinkelFilter is a class object that represents a polar
    %seperable frequency domain filter
    %
    % Based on the thesis by Michael van Ginkel. Chapter 3
    % "Image Analysis using Orientation Sapce based on Steerable Filters".
    % Delft University of Technology. October 7th, 2002.
    %
    % f_c: maximum frequency for the radial filter
    % b_f: frequency bandwidth for the radial filter
    % K: number of rotation angles through 360 degrees

    % Mark Kittisopikul, August 22nd, 2015
    % Jaqaman Lab
    % UT Southwestern
    
    properties (SetAccess = immutable)
        % radial central frequency
        f_c
        % radial frequency bandwidth
        b_f
        % angular order
        K
        % basis angles
        angles
    end
    
    properties (Transient)
        % filter size, should correspond with image size
        size
        % Filter itself
        F
    end
    
    properties (Dependent = true)
        % f_c
        centralFrequency
        % b_f
        frequencyBandwidth
        % K
        order
    end
    
    methods
        function obj = SteerableVanGinkelFilter(f_c,b_f,K)
            obj.f_c = f_c;
            obj.b_f = b_f;
            obj.K = K;
            
            n = 2*K + 1;
            obj.angles = 0:pi/n:pi-pi/n;
            
        end
        function f_c = get.centralFrequency(obj)
            f_c = obj.f_c;
        end
        function b_f = get.frequencyBandwidth(obj)
            b_f = obj.b_f;
        end
        function K = get.order(obj)
            K = obj.K;
        end
        function R = mtimes(obj,I)
            R = getResponse(obj,I);
        end
        function R = getResponse(obj,I)
            If = fft2(I);
            ridgeResponse = obj.applyRidgeFilter(If);
            edgeResponse = obj.applyEdgeFilter(If);
            angularResponse = ridgeResponse + edgeResponse;
            R = OrientationSpaceResponse(obj,angularResponse);
        end
        function R = getRidgeResponse(obj,I)
            If = fft2(I);
            ridgeResponse = obj.applyRidgeFilter(If);
            R = OrientationSpaceResponse(obj,ridgeResponse);
        end
        function R = getEdgeResponse(obj,I)
            If = fft2(I);
            edgeResponse = obj.applyEdgeFilter(If);
            R = OrientationSpaceResponse(obj,edgeResponse);
        end
        function imshow(obj,n,varargin)
            if(nargin < 2 || isempty(n))
                n = 1;
            end
            if(nargin < 3)
                varargin{1} = [];
            end
            imshow(fftshift(obj.F(:,:,n)),varargin{:});
        end
    end
    methods
        function setupFilter(obj,siz)
            if( isempty(obj.size) || siz ~= obj.size || isempty(obj.F))
                obj.size = siz;
                obj.F = orientationSpace.kernel(obj.f_c, obj.b_f, obj.K, obj.angles, obj.size);
            end
        end
        function ridgeResponse = applyRidgeFilter(obj,If)
            obj.setupFilter(size(If));
            ridgeResponse = real(ifft2(bsxfun(@times,If,real(obj.F))));
        end
        function edgeResponse = applyEdgeFilter(obj,If)
            obj.setupFilter(size(If));
            edgeResponse = 1j*real(ifft2(bsxfun(@times,If.*-1j,imag(obj.F))));
        end
    end
end

