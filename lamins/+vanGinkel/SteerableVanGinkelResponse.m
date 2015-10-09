classdef SteerableVanGinkelResponse < handle
    %SteerableVanGinkelResponse Response object for SteerableVanGinkelFilter
    %
    %
    
    properties
        filter
        angularResponse
        n
    end
    
    properties (Transient)
        angularGaussians
        matrix
    end
    
    properties (Transient, Access = protected)
        cache
    end
    
    properties (Dependent)
        res
        nms
        theta
        a
    end
    
    
    
    methods
        function obj = SteerableVanGinkelResponse(filter,angularResponse)
            obj.filter = filter;
            obj.angularResponse = angularResponse;
            obj.n = size(angularResponse,3);
        end
        
        % Response, defaults to maximum response with initial basis
        function set.res(obj,res)
            obj.cache.res = res;
        end
        function res = get.res(obj)
            if(~isfield(obj.cache,'res') || isempty(obj.cache.res))
                obj.cache.res = getMaxResponse(obj);
            end
            res = obj.cache.res;
        end
        
        % Non-Maximal Suppression, defaults to nms based on initial basis
        function set.nms(obj,nms)
            obj.cache.nms = nms;
        end
        function nms = get.nms(obj)
            if(~isfield(obj.cache,'nms') || isempty(obj.cache.nms))
                [~,obj.cache.nms] = nonMaximumSuppression(obj.res,obj.theta);
            end
            nms = obj.cache.nms;
        end
        
        % Orientation, defaults to best orientation from basis
        function set.theta(obj,theta)
            obj.cache.theta = theta;
        end
        function theta = get.theta(obj)
            if(~isfield(obj.cache,'theta') || isempty(obj.cache.theta))
                [~,obj.cache.theta] = getMaxResponse(obj);
            end
            theta = obj.cache.theta;
        end
        
        % Shortcut for angularResponse
        function set.a(obj,a)
            obj.angularResponse = a;
        end
        function a = get.a(obj)
            a = obj.angularResponse;
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
        function M = getMatrix(obj)
            if(isempty(obj.matrix))
                A = obj.getAngularGaussians;
                aR = obj.angularResponse;
                obj.matrix = reshape(aR,size(aR,1)*size(aR,2),size(aR,3))/A;
            end
            M = obj.matrix;
        end
        function [response,samples] = getResponseAtPoint(obj,r,c,samples)
            if(nargin < 4)
                samples = obj.n;
            end
            siz = size(obj.angularResponse);
            linIdx = sub2ind(siz([1 2]),r,c);
            if(isscalar(samples) && samples == obj.n)
                a = reshape(obj.angularResponse,[],obj.n);
                response = squeeze(a(linIdx,:));
            else
                if(isscalar(samples))
                    if(mod(samples,1) == 0)
                        samples = (0:samples-1)'*obj.n/samples;
                    else
                        samples = pi/samples;
                        samples = (0:samples-1)'*obj.n/samples;
                    end
                end
                M = obj.getMatrix;
                sampling = bsxfun(@minus,0:obj.n-1,samples);
                sampling = wraparoundN(sampling,-obj.n/2,obj.n/2)';
                sampling = exp(-sampling.^2/2);
                response = M(linIdx,:)*sampling;
            end
        end
        function response = getResponseAtOrientation(obj,angle)
        %getResponseAtOrientation(orientationAngle)
            % Returns the response image plane at the orientation specified
            % in radians
            M = obj.getMatrix;
            if(isinteger(angle) && angle ~= 0)
                %do nothing
                angle = double(angle);
            else
                angle = double(angle)/pi*obj.n;
            end
            xt = bsxfun(@minus,0:obj.n-1,angle)';
            xt = wraparoundN(xt,-obj.n/2,obj.n/2);
            response = M*exp(-xt.^2/2);
            response = reshape(response,size(obj.angularResponse,1),size(obj.angularResponse,2));
            if(abs(wraparoundN(angle,-pi,pi))/pi > 0.5)
                response = real(response) - 1j*imag(response);
            end
        end
        function [response,theta] = getMaxResponse(obj,nn)
            import vanGinkel.*;
            if(nargin < 2)
                nn = obj.n;
            end
            if(isfinite(nn))
                [response,theta] = obj.getMaxFiniteResponse(nn);
            else
                [theta,response] = vanGinkelMaxima(obj.angularResponse);
            end
        end
        function [response,theta] = getMaxFiniteResponse(obj,nn)
            import vanGinkel.*;
            a = obj.angularResponse;
            if(obj.n ~= nn)
                vanGinkelUpsample(a,pi/nn);
            end
            [response,theta] = max(real(a),[],3);
            [response_i,theta_i] = max(cat(3,imag(a),-imag(a)),[],3);
            response = response +1j*response_i;
            theta = theta + 1j*theta_i;
        end
        function varargout = subsref(obj,S)
            switch(S(1).type)
%                 case '.'
                case '()'
                    varargout{1} = obj.getResponseAtPoint(S(1).subs{:});
                case '{}'
                    varargout{1} = obj.getResponseAtOrientation(S(1).subs{:});
                otherwise
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
            end
        end
        function h = imshow(obj,varargin)
            h = imshow(obj.getMaxResponse,varargin{:});
        end
        function h = plot(obj,angles,r,c,varargin)
            [Y,samples] = obj.getResponseAtPoint(r,c,angles);
            h = plot(samples/obj.n,Y,varargin{:});
        end
    end
    
end

