function newG=EnhanceGrid(oldG)

    [Ny,Nx] = size(oldG);
    Tx=2*Nx-1;
    Ty=2*Ny-1;
    newG=zeros(Ty,Tx);

    newG(1:2:Ty,1:2:Tx)=oldG;
    newG(2:2:(Ty-1),:)=(newG(1:2:(Ty-2),:)+newG(3:2:Ty,:))/2;
    newG(:,2:2:(Tx-1))=(newG(:,1:2:(Tx-2))+newG(:,3:2:Tx))/2;
    
