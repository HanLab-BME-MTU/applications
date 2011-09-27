function distMat = bwMaxDirectDist3D(mask)

%Ill explain later... For now, just trust me.

[M,N,P] = size(mask);
%mask = repmat(mask,[1 1 1 6]);

%WTF - For some reason it's MUCH faster with double or single than integer class!!
distMat = zeros([size(mask) 6],'single');

distMat(1,:,:,1) = M;
distMat(M,:,:,2) = M;
distMat(:,1,:,3) = N;
distMat(:,N,:,4) = N;
distMat(:,:,1,5) = P;
distMat(:,:,P,6) = P;

%This is even slower....
% for m = 2:M       
%     %distMat(m,~mask(m,:,:),1) = distMat(m-1,~mask(m,:),1) + 1;       
%     distMat(m,:,:,1) = distMat(m-1,:,:,1) + 1;
%     distMat(mask) = 0;
% end

%Gotta be a better way to do this!!! Tooo sloww....
for n = 1:N
    for m = 2:M       
        distMat(m,n,~mask(m,n,:),1) = distMat(m-1,n,~mask(m,n,:),1) + 1;       
    end
    for m = M-1:-1:1       
        distMat(m,n,~mask(m,n,:),2) = distMat(m+1,n,~mask(m,n,:),2) + 1;       
    end    
end
for m = 1:M
    for n = 2:N       
        distMat(m,n,~mask(m,n,:),3) = distMat(m,n-1,~mask(m,n,:),3) + 1;       
    end
    for n = N-1:-1:1       
        distMat(m,n,~mask(m,n,:),4) = distMat(m,n+1,~mask(m,n,:),4) + 1;       
    end
end
for m = 1:M
    for p = 2:P       
        distMat(m,~mask(m,:,p),p,5) = distMat(m,~mask(m,:,p),p-1,5) + 1;       
    end
    for p = P-1:-1:1       
        distMat(m,~mask(m,:,p),p,6) = distMat(m,~mask(m,:,p),p+1,6) + 1;       
    end
end

%FASTER WAY TO DO THIS????
% maxMat = cat(4,ones(size(mask))*M,ones(size(mask))*M,ones(size(mask))*N,ones(size(mask))*N,ones(size(mask))*P,ones(size(mask))*P);
% distMat(distMat > maxMat) = maxMat(distMat > maxMat);
