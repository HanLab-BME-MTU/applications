% Conv helper

conv1=convTool(256,7,2,1)
max1=MaxPoolTool(conv1,2,2,1)

function outW2=convTool(width,kW,dW,padW)
% if nargin == 0
%     width = 256; %(image input)
% end
    
% padW = 2;    %(padding)
% kW = 4;      %(kernel size kWxkW)
% dW = 2;      %(step) 
outW2  = floor((width  + 2*padW - kW) / dW + 1);
outW=num2str(outW2);
disp(['Conv2DInput: [',num2str(width),'x',num2str(width),'] -->  [',outW,'x',outW,']']) 
end


function outW2=MaxPoolTool(width,kW,dW,padW)
    outW2  = floor((width  + 2*padW - kW) / dW + 1);
    outW=num2str(outW2);
    disp(['MaxPoolInput: [',num2str(width),'x',num2str(width),'] -->  [',outW,'x',outW,']']) 
    
end
