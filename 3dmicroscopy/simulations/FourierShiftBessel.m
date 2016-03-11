

load('BesselEFields.mat','BesselEField2');

numBeams=7;

matSize=size(BesselEField2);

CoherentBesselEField=zeros(matSize(1),3*matSize(1));

CoherentBesselEField(:,matSize(1)+1:2.*matSize(1))= CoherentBesselEField(:,matSize(1)+1:2.*matSize(1))+BesselEField2;

counter=0;

for i=1:.10:200
    
    counter=counter+1;
    
    BesselFinal=zeros(matSize(1),matSize(1).*3);
    
    for j=-3:1:3
        
        BesselShift=FourierShift2D(CoherentBesselEField,[0 j.*i]);
        
        BesselFinal=BesselFinal+BesselShift;
        
    end
    
    BesselFinal=abs(BesselFinal).^2;

    filename = [ 'myDataFile' num2str(counter) '.mat' ];
    
    save(filename,'BesselFinal', '-v7.3')
        
    clear BesselFinal
    
end

