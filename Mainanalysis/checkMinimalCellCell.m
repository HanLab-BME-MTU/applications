                      
function coord1= checkMinimalCellCell(coordnuc,altercoor,MinDistCellCell)                      
                      
                      


%ensure a minimal distance between them
count=1;
while count < length(coordnuc)
       paff=[];
       paff= min(sqrt((coordnuc(count+1:end,1)-coordnuc(count,1)).^2 + (coordnuc(count+1:end,2)-coordnuc(count,2)).^2));
       if paff < MinDistCellCell
               coordnuc(count,:)=[];
               count=count-1;
       end
       count=count+1;
end
clear count;    
clear paff;


%ensure minimal distance between them
count=1;
while count < length(altercoor)
       paffalt=[];
       paffalt= min(sqrt((altercoor(count+1:end,1)-altercoor(count,1)).^2 + (altercoor(count+1:end,2)-altercoor(count,2)).^2));

       if paffalt < MinDistCellCell
               altercoor(count,:)=[];
               count=count-1;
       end
       count=count+1;
end
clear paffalt;
clear count;


%ensure minimal distance between dark, nasty-looking cells and
%third eyes (look above(findnucloitrack,halosfind))
namesnumbers=[0,0];
if isempty(altercoor) == 0
    for h=1:size(altercoor,1)
            uff=[];
             
            uff= min(sqrt((coordnuc(:,1)-altercoor(h,1)).^2+(coordnuc(:,2)-altercoor(h,2)).^2));
            %note that here the minimal distance is larger than
            %between two cells found by the same routine,
            %because... ahhhmmm... just to make sure
            if uff  >  (1.5*MinDistCellCell)
                  namesnumbers(end+1,1)=altercoor(h,1);
                  namesnumbers(end,2)=altercoor(h,2);
            end
     end
     namesnumbers(1,:)=[];
     
     if isempty(namesnumbers) == 0
            coord1=cat(1,coordnuc,namesnumbers);
     else
            coord1=coordnuc;
     end
         
else
      coord1=coordnuc;
end
clear altercoor;    
clear namesnumbers;   
clear coordnuc;   
clear uff;


