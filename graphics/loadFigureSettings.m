function fset = loadFigureSettings(mode)

% light blue
fset.cfB = hsv2rgb([5/9 0.3 1]);
fset.ceB = hsv2rgb([5/9 1 1]);

% dark blue
fset.cfB2 = hsv2rgb([0.6 0.3 1]);
fset.ceB2 = hsv2rgb([0.6 1 1]);


fset.cfR = hsv2rgb([0 0.3 1]);
fset.ceR = hsv2rgb([0 1 1]);

fset.cfG = hsv2rgb([1/3 0.3 1]);
fset.ceG = hsv2rgb([1/3 1 1]);

fset.cfY = hsv2rgb([1/7 0.3 1]);
fset.ceY = hsv2rgb([1/7 1 1]);

fset.cf0 = hsv2rgb([0 0 0.7]);
fset.ce0 = hsv2rgb([0 0 0.3]);


fset.fontName = {'FontName', 'Helvetica'};

switch mode
    case 'print'
        fset.ifont = [fset.fontName, 'FontSize', 7];
        fset.tfont = [fset.fontName, 'FontSize', 8];
        fset.sfont = [fset.fontName, 'FontSize', 10];
        fset.lfont = [fset.fontName, 'FontSize', 12];
        fset.axOpts = ['Layer', 'top', 'TickDir', 'out', 'LineWidth', 1.5, fset.sfont];
    otherwise
        fset.ifont = [fset.fontName, 'FontSize', 12];
        fset.tfont = [fset.fontName, 'FontSize', 16];
        fset.sfont = [fset.fontName, 'FontSize', 20];
        fset.lfont = [fset.fontName, 'FontSize', 24];
        fset.axOpts = ['Layer', 'top', 'TickDir', 'out', 'LineWidth', 2, fset.sfont];
end

ce = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
% ce = [0 1 0; 1 0 0; 1 0.5 0; 1 1 0; 0 1 1; 0.5 0 1; 0 0 1; 0 0.5 1];
cf = rgb2hsv(ce);
cf(:,2) = 0.3;
cf = hsv2rgb(cf);
fset.ceTrackClasses = ce;
fset.cfTrackClasses = cf;