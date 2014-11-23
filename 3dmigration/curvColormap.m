function cMap = curvColormap(nCol,type,colors)
if nargin <1 || isempty(nCol);
    nCol = 256;
else
    nCol = 2*round(nCol/2);
end

if nargin < 2 || isempty(type)
    type = '3col';
end

if nargin < 3 || isempty(colors)
    %Colormap for emphasizing zero-ish values
    colors = [.4 .4 .4; 1 0 0 ; .4 .4 .2];
end
nPer = nCol/2;
%quick and dirty way to have a common colormap accross functions
switch type
    
    
    case 'iso'
        a = isomorphicColormap('r',nCol);
        a = bsxfun(@minus,a,a(1,:)); 
        b = isomorphicColormap('b',nCol);
        b = bsxfun(@minus,b,b(1,:));

        cMap = cat(1,a(round(end/2):-1:1,:),b(1:round(end/2),:));
        
    case 'kMid'
        
        
        a = [zeros(nPer,2) linspace(1,0,nPer)'];
        b = [linspace(0,1,nPer)' zeros(nPer,2)];
        
        cMap = cat(1,a,b);
        
    case '3col'
        x = linspace(-2,2,nCol);
        y = x ./ sqrt(1+ x .^2);%If it's linear the majority of the scale is taken up by red.
        %plot(x,y)
        cMap = interp1([min(y) 0 max(y)]',colors,y); 
        
        
end
