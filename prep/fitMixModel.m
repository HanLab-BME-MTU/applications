function [estimates] = fitMixModel(image, xyvec, sig, amps, b)
% fit mixture model to image with Jacobian, using info from cands as start
% INPUT: image
%        xyvec: (numspec x 2) vector containing startpoints for x,y, 
%               coordinates of speckles 
%       sig: sigma of Guassian, corresponding to (1/3)*PSF width
%       amps: amplitudes (deltaI) vector
%       b; background value


[xs,ys]=size(image);
[numspec,nn]=size(xyvec);

%speckle positions
mspos=xyvec;
%speckle background level
msIBkg=zeros(numspec,1);
msIBkg(:)=b;
%speckle intensity elevation over background
msdeltaI=zeros(numspec,1);
msdeltaI(:)=amps(:);
mssigma=sig;

% IN ORDER TO FIT USING SPARSE MATRIX, FIRST IDENTIFY RELEVANT AREAS 
% i.e., the ones SURROUNDING SPECKLES 
% make mask of image with value 1 in a (3*sig radius) circle around 
% speckles and value 0 elsewhere
% step 1: mat1=xs*ys x 2 matrix with x,y coordinates of all consecutive
% pixels in image
mat1=zeros(xs*ys,2);
[mat1(:,1),mat1(:,2)]=ind2sub([xs,ys],1:(xs*ys));
mat2=mat1;
mat3=zeros(xs*ys,1);

%mask is copy of image with 1 for speckle vicinity; max projection of all 
%single mat3's (calculated in loop below)
mask=zeros(xs*ys,1);
%masksin = sparse matrix with 1/0 values for coordinates; each speckle is
%one column
masksin=sparse(xs*ys,numspec);
for n=1:numspec
    mat2(:,1:2)=[mat1(:,1)-mspos(n,1) mat1(:,2)-mspos(n,2)];
    mat3=sqrt( ( mat2(:,1).^2 + (mat2(:,2)).^2 ));
    % set=1 for distances below 3*sig, 0 for greater distances
    % masksinmasksin=sparse matrix with val>0 in mat3 defined
    mat3(mat3>(3*sig))=0;
    mat3(mat3>0)=1;
%     imshow(reshape(mat3,xs,ys),[]);
%     pause(0.001);
    masksin(find(mat3),n)=1;
    mask=mask+mat3;
end

%all positions in masksin where value is nonzero
nonzMaskVec=find(masksin);
valVec=nonzMaskVec;
valVec(:)=1;
mask(mask>0)=1;

% startpoint for fit
start_point = [mspos msdeltaI];
%lower_upper_boundary=[lowbound; upbound]
lowbound=[mspos-2 0.8*msdeltaI]; 
upbound=[mspos+2 1.2*msdeltaI];  

options = optimset('TolFun',1e-3,'TolX',1e-3,'Display','iter','MaxIter',15,'Jacobian','on');
checker=0;
checker2=0;

[estimates,renorm,residual,exitflag,output]= lsqnonlin(@Gauss2Dfun, start_point, lowbound, upbound, options);

% in the present implementation, the function Gauss2Dfun fits the x-
% and y-coordinates of speckles and the speckle amplitude. 

% the returned function is the function difference F-image, NOT the squared
% error; the square minimization is done implicitly in the lsqnonlin
% function

    function [F,J] = Gauss2Dfun(params)
        x0 = params(:,1);
        y0 = params(:,2);
        bav = msIBkg;
        A = params(:,3);
        w0 = mssigma;
        [ix,iy]=size(params);
        % mat 1 is defined above (i.e. outside the function) and contains the x,y-
        % coordinates of all (:)-consecutive pixels in the image
        imxyvar=mat1;
        %loop over number of speckles
        for z=1:ix
            zz=( (imxyvar(:,1)-x0(z)).^2 + (imxyvar(:,2)-y0(z)).^2 );
            FittedCurveS = A(z) * exp( -(1/(2*(w0^2))) *  zz );
            if(z==1)
                FittedCurve=FittedCurveS;
            else
                FittedCurve=FittedCurve+FittedCurveS;
            end
        end

        % multiply by mask (also pre-defined outside the actual fitting
        % function to save time) in order to take into account for the fit
        % only the relevant areas around speckles
        % and not try to fit unspecific background
        F = (FittedCurve(:)+bav(1)).*mask(:) - image(:).*mask(:);
        
        checker2=checker2+1;
        
        % now calculate the Jacobian matrix!!
        
        if (nargout>1)
            % fun returns vector of m=(xs*ys) components (corresponding to 
            % image size), and parameter vector x has n=(numspec x 3) components, 
            % namely the x and y coordinates for numspec speckles and their amplitudes; thus,
            % Jacobian is m x n matrix containing partial derivatives
            [sx,sy]=size(start_point);
            J=sparse( (xs*ys),(sx*sy) );
            % reserve space in the sparse Jacobian:
            % for each differentiating parameter, fill up the J-coordinates
            % corresponding to speckle coordinates in the image with the 
            % value one. The speckle coordinate positions (3-sig circle) are 
            % stored in the matrix masksin calculated outside the function.
            % The first parameter x=(x1,x2,x3,...xnumspec) corresponds to 
            % the columns 1:numspec in the Jacobian
            J(:,1:numspec)=masksin;
            % same for the second parameter: y=(y1,y2,y3,....ynumspec)
            % corresponds to columns (numspec+1):end
            J(:,(numspec+1):2*numspec)=masksin;
            % ...and same for third parameter, amplitude I0
            J(:,(2*numspec+1):3*numspec)=masksin;
            
            %if desired, display the Jacobian matrix outline, reshaped as image matrix,
            %side-by-side for the numspec x 2 different elements
%             if(checker==0)
%                 spy(reshape(J,xs,numspec*2*ys));
%                 pause(0.2);
%                 checker=1;
%             end   %of if
            
            %now calculate actual values for Jacobian for each speckle
            for s=1:numspec
                %general procedure:
                %make left- and right-shifted image, subtract and divide by
                %two, according to
                %dx/dt=[(dx+Dt)-(dx-Dt)] / (2*Dt)
                                
                %v (1/0)-value position vector for speckle s, s-th column of masksin
                v=masksin(:,s);
                %vind index vector of v, positions for which v>0
                %xx,yy x and y-coordinates of nonzero points for this
                %speckle
                vind=find(v);
                %ind=1:(xs*ys);
                [xx,yy]=ind2sub([xs,ys],vind);
                
                %shift distance
                sh=0.3;
                
                %left shifted image of speckle s
                %not entire image is calculated, but just projection onto
                %nonzeros in v (due to structure of xx and yy above)
                zzleft=( (xx-x0(s)-sh).^2 + (yy-y0(s)).^2 );
                varleft = A(z) * exp( -(1/(2*(w0^2))) * ( zzleft) );
                timageleft=varleft-nonzeros(image(:).*v) + bav(1);
                %right shifted image of speckle s
                zzright=( (xx-x0(s)+sh).^2 + (yy-y0(s)).^2 );
                varright = A(z) * exp( -(1/(2*(w0^2))) * ( zzright) );
                timageright=varright-nonzeros(image(:).*v) + bav(1);
                %central difference image x (left-right)
                diff1=(-timageright+timageleft)/2;
                J(vind,s)=diff1;
                
                %up-shifted image of speckle s
                zztop=( (xx-x0(s)).^2 + (yy-y0(s)+sh).^2 );
                vartop = A(z) * exp( -(1/(2*(w0^2))) * ( zztop) );
                timagetop=vartop-nonzeros(image(:).*v) + bav(1);
                %down-shifted image of speckle s
                zzbottom=( (xx-x0(s)).^2 + (yy-y0(s)-sh).^2 );
                varbottom = A(z) * exp( -(1/(2*(w0^2))) * ( zzbottom) );
                timagebottom=varbottom-nonzeros(image(:).*v) + bav(1);
                %central difference image y (bottom-top)
                diff2=(-timagetop+timagebottom)/2;
                J(vind,s+numspec)=diff2;
                
                %differentiation df/da (amplitude)
                zzamp=( (xx-x0(s)).^2 + (yy-y0(s)).^2 );
                diff3=exp( -(1/(2*(w0^2))) * ( zzamp) );
                J(vind,s+2*numspec)=diff3;
                                
            end  %of for s=1:numspec
                       
        end % of if (nargout>1)
        
    end %of function F

output    %output of fit, to check for errors
disp(exitflag);

end  %of function  
