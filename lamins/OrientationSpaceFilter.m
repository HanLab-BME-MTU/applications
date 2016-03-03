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
        % Angular gaussians useful for manipulating response
        angularGaussians
    end
    
    properties (Dependent = true)
        % f_c
        centralFrequency
        % b_f
        frequencyBandwidth
        % K
        order
        % number of angular filter templates
        n
    end
    
    methods
        function obj = OrientationSpaceFilter(f_c,b_f,K)
            obj.f_c = f_c;
            obj.b_f = b_f;
            obj.K = K;
            
            n = 2*K + 1;
            obj.angles = 0:pi/n:pi-pi/n;
            
        end
        function ridgeFilter = real(obj)
            ridgeFilter = OrientationSpaceRidgeFilter(obj.f_c,obj.b_f,obj.K);
            % The filter itself does not change (for the moment)
            ridgeFilter.F = obj.F;
            ridgeFilter.size = obj.size;
        end
        function edgeFilter = imag(obj)
            edgeFilter = OrientationSpaceEdgeFilter(obj.f_c,obj.b_f,obj.K);
            % The filter itself does not change (for the moment)
            edgeFilter.F = obj.F;
            edgeFilter.size = obj.size;
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
        function n = get.n(obj)
            n = obj.K*2+1;
        end
        function R = mtimes(obj,I)
            % Convolution
            if(isa(obj,'OrientationSpaceFilter'))
                R = getResponse(obj,I);
            elseif(isa(I,'OrientationSpaceFilter'))
                % The convolution is commutative, swap the parameters
                R = getResponse(I,obj);
            end
        end
        function R = getResponse(obj,I)
            If = fft2(I);
            ridgeResponse = obj.applyRidgeFilter(If);
            edgeResponse = obj.applyEdgeFilter(If);
            angularResponse = ridgeResponse + edgeResponse;
            R(numel(obj)) = OrientationSpaceResponse;
            for o=1:numel(obj)
                R(o) = OrientationSpaceResponse(obj(o),angularResponse(:,:,:,o));
            end
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
        function A = getAngularGaussians(obj)
            if(isempty(obj.angularGaussians))
                N = obj.n;
                x = 0:N-1;
                xx = bsxfun(@minus,x,x');
                xx = wraparoundN(xx,-N/2,N/2);
                obj.angularGaussians = exp(-xx.^2/2);
            end
            A = obj.angularGaussians;
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
            for o=1:numel(obj)
                if( isempty(obj(o).size) || any(siz ~= obj(o).size) || isempty(obj(o).F))
                    obj(o).size = siz;
                    obj(o).F = orientationSpace.kernel(obj(o).f_c, obj(o).b_f, obj(o).K, obj(o).angles, obj(o).size);
                end
            end
        end
        function ridgeResponse = applyRidgeFilter(obj,If)
            obj.setupFilter(size(If)); %#ok<CPROP>
            ridgeResponse = real(ifft2(bsxfun(@times,If,real(cat(4,obj.F)))));
        end
        function edgeResponse = applyEdgeFilter(obj,If)
            obj.setupFilter(size(If)); %#ok<CPROP>
            edgeResponse = 1j*real(ifft2(bsxfun(@times,If.*-1j,imag(cat(4,obj.F)))));
        end
    end
end

