function [ h2 ] = imThetaOffsetLine( h , varargin)
%imThetaOffsetLine Creates an imline offset only by theta

dim.x = xlim;
dim.y = ylim;

center = [mean(dim.x) mean(dim.y)];

hf = gcf;
h2 = lamins.simsim.imRadialLine(varargin{:});
deltaTheta = pi/6;
offsetLine(h,h2,pi/6);

addNewPositionCallback(h,@doOffset);
addNewPositionCallback(h2,@offsetLineCallback);

arrayfun(@(h2) iptaddcallback(h2,'ButtonDownFcn',@resumeOnDoubleClick),get(h2,'Children'));

% while(isvalid(h2))
%     wait(h2);
%     resume(h);
% end

     function resumeOnDoubleClick(varargin)
        if strcmp(get(hf,'SelectionType'),'open');
            resume(h);
        end
     end

    function offsetLineCallback(pos)
        [theta(1),rho(1)] = cart2pol(pos(1,1)-center(1),pos(1,2)-center(2));
        [theta(2),rho(2)] = cart2pol(pos(2,1)-center(1),pos(2,2)-center(2));
        
        theta = mean(theta);
        
        pos = getPosition(h);
        [theta(2)] = cart2pol(pos(1,1)-center(1),pos(1,2)-center(2));
        [theta(2)] = cart2pol(pos(2,1)-center(1),pos(2,2)-center(2));
        
        [pos(1,1), pos(1,2)] = pol2cart(theta(2),rho(1));
        [pos(2,1), pos(2,2)] = pol2cart(theta(2),rho(2));
        setConstrainedPosition(h,pos+repmat(center,2,1));
        
        deltaTheta = theta(2) - theta(1);
    end
    function doOffset(varargin)
        offsetLine(h,h2,deltaTheta);
    end
    function offsetLine(hRef,hOffset,deltaTheta)
        if(~isvalid(hOffset))
            return;
        end
        pos = getPosition(hRef);
        [theta(1),rho(1)] = cart2pol(pos(1,1)-center(1),pos(1,2)-center(2));
        [theta(2),rho(2)] = cart2pol(pos(2,1)-center(1),pos(2,2)-center(2));
        
        theta = mean(theta) - deltaTheta;
        
        [pos(1,1), pos(1,2)] = pol2cart(theta,rho(1));
        [pos(2,1), pos(2,2)] = pol2cart(theta,rho(2));
        setPosition(hOffset,pos+repmat(center,2,1));
    end


end

