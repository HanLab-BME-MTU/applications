% FEMLAB Model M-file
% Generated 18-Feb-2004 14:50:23 by FEMLAB 2.3.0.153.

flclear fem
% FEMLAB Version
clear vrsn;
vrsn.name='FEMLAB 2.3';
vrsn.major=0;
vrsn.build=153;
fem.version=vrsn;

% Recorded command sequence

% New geometry 1
fem.sdim={'x','y'};

% Geometry
clear s c p
R1=rect2(-200,200,-150,150,0);
E1=ellip2(0,0,50,50,0);
objs={R1,E1};
names={'R1','E1'};
s.objs=objs;
s.name=names;

objs={};
names={};
c.objs=objs;
c.name=names;

objs={};
names={};
p.objs=objs;
p.name=names;

drawstruct=struct('s',s,'c',c,'p',p);
fem.draw=drawstruct;
fem.geom=geomcsg(fem);

clear appl

% Application mode 1
appl{1}.mode=flsmps('dim',{'u','v','u_t','v_t'},'sdim',{'x','y'},'submode', ...
'std','tdiff','on');
appl{1}.dim={'u','v','u_t','v_t'};
appl{1}.form='coefficient';
appl{1}.border='off';
appl{1}.name='ps';
appl{1}.var={'omega','2*pi*100'};
appl{1}.assign={'E';'E';'Fx';'Fx';'FxAmp';'FxAmp';'FxPh';'FxPh';'Fxphdep'; ...
'Fxphdep';'Fy';'Fy';'FyAmp';'FyAmp';'FyPh';'FyPh';'Fyphdep';'Fyphdep';'Hx'; ...
'Hx';'Hy';'Hy';'Rx';'Rx';'Ry';'Ry';'Tax';'Tax';'Tay';'Tay';'Temp';'Temp'; ...
'Tempref';'Tempref';'Tflag';'Tflag';'Tlx';'Tlx';'Tly';'Tly';'alpha';'alpha'; ...
'alphadM';'alphadM';'betadK';'betadK';'disp';'disp';'e1';'e1';'e2';'e2'; ...
'e3';'e3';'eigfreq';'eigfreq';'eiglambda';'eiglambda';'ex';'ex';'ex_amp'; ...
'ex_amp';'ex_ph';'ex_ph';'exy';'exy';'exy_amp';'exy_amp';'exy_ph';'exy_ph'; ...
'ey';'ey';'ey_amp';'ey_amp';'ey_ph';'ey_ph';'ez';'ez';'ez_amp';'ez_amp'; ...
'ez_ph';'ez_ph';'mises';'mises';'nu';'nu';'omega';'omega';'rho';'rho';'s1'; ...
's1';'s2';'s2';'s3';'s3';'sx';'sx';'sx_amp';'sx_amp';'sx_ph';'sx_ph';'sxy'; ...
'sxy';'sxy_amp';'sxy_amp';'sxy_ph';'sxy_ph';'sy';'sy';'sy_amp';'sy_amp'; ...
'sy_ph';'sy_ph';'thickness';'thickness';'tresca';'tresca';'u_amp';'u_amp'; ...
'u_ph';'u_ph';'u_t';'u_t';'u_t_amp';'u_t_amp';'u_t_ph';'u_t_ph';'u_tt'; ...
'u_tt';'u_tt_amp';'u_tt_amp';'u_tt_ph';'u_tt_ph';'v_amp';'v_amp';'v_ph'; ...
'v_ph';'v_t';'v_t';'v_t_amp';'v_t_amp';'v_t_ph';'v_t_ph';'v_tt';'v_tt'; ...
'v_tt_amp';'v_tt_amp';'v_tt_ph';'v_tt_ph'};
appl{1}.elemdefault='Lag2';
appl{1}.shape={'shlag(2,''u'')','shlag(2,''v'')'};
appl{1}.sshape=2;
appl{1}.eigtype='freq';
appl{1}.equ.E={'2+sin(2*pi*x/50)','2+sin(2*pi*x/50)'};
appl{1}.equ.nu={'0.33','0.33'};
appl{1}.equ.alpha={'0','0'};
appl{1}.equ.rho={'0','0'};
appl{1}.equ.alphadM={'0','0'};
appl{1}.equ.betadK={'0','0'};
appl{1}.equ.Fx={'0','1e-5*sqrt(x^2+y^2)'};
appl{1}.equ.Fy={'0','-1e-5*sqrt(x^2+y^2)'};
appl{1}.equ.FxAmp={'1','1'};
appl{1}.equ.FyAmp={'1','1'};
appl{1}.equ.FxPh={'0','0'};
appl{1}.equ.FyPh={'0','0'};
appl{1}.equ.thickness={'1','1'};
appl{1}.equ.loadtype={'area','area'};
appl{1}.equ.constrtype={'standard','standard'};
appl{1}.equ.H={{{'0'},{'0'};{'0'},{'0'}},{{'0'},{'0'};{'0'},{'0'}}};
appl{1}.equ.R={{{'0'};{'0'}},{{'0'};{'0'}}};
appl{1}.equ.Hx={0,0};
appl{1}.equ.Hy={0,0};
appl{1}.equ.Rx={'0','0'};
appl{1}.equ.Ry={'0','0'};
appl{1}.equ.Tflag={0,0};
appl{1}.equ.Temp={'0','0'};
appl{1}.equ.Tempref={'0','0'};
appl{1}.equ.gporder={{4;4},{4;4}};
appl{1}.equ.cporder={{2;2},{2;2}};
appl{1}.equ.shape={[1 2],[1 2]};
appl{1}.equ.init={{{'0'};{'0'}},{{'0'};{'0'}}};
appl{1}.equ.usage={1,1};
appl{1}.equ.ind=[1 2];
appl{1}.bnd.coord={'global','global','global'};
appl{1}.bnd.constrtype={'standard','standard','standard'};
appl{1}.bnd.Fx={'1e-2*nx','0','0'};
appl{1}.bnd.Fy={'1e-2*ny','0','0'};
appl{1}.bnd.FxAmp={'1','1','1'};
appl{1}.bnd.FyAmp={'1','1','1'};
appl{1}.bnd.FxPh={'0','0','0'};
appl{1}.bnd.FyPh={'0','0','0'};
appl{1}.bnd.Hx={0,1,0};
appl{1}.bnd.Hy={0,1,0};
appl{1}.bnd.Rx={'0','-2*nx','0'};
appl{1}.bnd.Ry={'0','-2*ny','0'};
appl{1}.bnd.H={{{'0'},{'0'};{'0'},{'0'}},{{'0'},{'0'};{'0'},{'0'}},{{'0'}, ...
{'0'};{'0'},{'0'}}};
appl{1}.bnd.R={{{'0'};{'0'}},{{'0'};{'0'}},{{'0'};{'0'}}};
appl{1}.bnd.thickness={'1','1','1'};
appl{1}.bnd.loadtype={'length','length','length'};
appl{1}.bnd.gporder={{0;0},{0;0},{0;0}};
appl{1}.bnd.cporder={{0;0},{0;0},{0;0}};
appl{1}.bnd.shape={0,0,0};
appl{1}.bnd.ind=[1 2 2 1 3 3 3 3];
appl{1}.pnt.constrtype={'standard'};
appl{1}.pnt.Fx={'0'};
appl{1}.pnt.Fy={'0'};
appl{1}.pnt.FxAmp={'1'};
appl{1}.pnt.FyAmp={'1'};
appl{1}.pnt.FxPh={'0'};
appl{1}.pnt.FyPh={'0'};
appl{1}.pnt.Hx={0};
appl{1}.pnt.Hy={0};
appl{1}.pnt.Rx={'0'};
appl{1}.pnt.Ry={'0'};
appl{1}.pnt.H={{{'0'},{'0'};{'0'},{'0'}}};
appl{1}.pnt.R={{{'0'};{'0'}}};
appl{1}.pnt.shape={0};
appl{1}.pnt.ind=ones(1,8);

fem.appl=appl;

% Material library
clear lib
lib.Mat1.E='2';
lib.Mat1.nu='0.33';
lib.Mat1.alpha='0';
lib.Mat1.rho='0';
lib.Mat1.type='material';
fem.lib=lib;
