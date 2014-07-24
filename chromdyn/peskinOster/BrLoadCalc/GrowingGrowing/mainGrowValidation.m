% Main: growing optimization validation
clear
delta= 8e-3/13;





opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',300000000,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);


for i=1:122
    
    [velData1,velData2,velDataNoise51,velDataNoise52,velDataNoise451,velDataNoise452,cAngle1,cAngle2,omega1,omega2,v01,v02,xObj]=DataGenerationGrow(3000);
    
	x0=[0 0 0 0];
    opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',300,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);
	[fct3000(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:3000),velDataNoise452(1:3000),cAngle1(1:3000),cAngle2(1:3000),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velData1,velData2,cAngle1,cAngle2,delta,v01,v02);
	maxFct3000(i)=max(fct3000(:,i));sumFct3000(i)=sum(fct3000(:,i))/3000;
	
	[xRes3000(i,:),objVal3000(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:3000),velDataNoise452(1:3000),cAngle1(1:3000),cAngle2(1:3000),delta,v01,v02);
	maxObj3000(i)=max(objVal3000(:,i));sumObj3000(i)=sum(objVal3000(:,i))/3000;
    
	[fct2000(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:2000),velDataNoise452(1:2000),cAngle1(1:2000),cAngle2(1:2000),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velData1,velData2,cAngle1,cAngle2,delta,v01,v02);
	maxFct2000(i)=max(fct2000(:,i));sumFct2000(i)=sum(fct2000(:,i))/2000;
	
	[xRes2000(i,:),objVal2000(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:2000),velDataNoise452(1:2000),cAngle1(1:2000),cAngle2(1:2000),delta,v01,v02);
	maxObj2000(i)=max(objVal2000(:,i));sumObj2000(i)=sum(objVal2000(:,i))/2000;

    
	[fct1000(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:1000),velDataNoise452(1:1000),cAngle1(1:1000),cAngle2(1:1000),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct1000(i)=max(fct1000(:,i));sumFct1000(i)=sum(fct1000(:,i))/1000;
	
	[xRes1000(i,:),objVal1000(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:1000),velDataNoise452(1:1000),cAngle1(1:1000),cAngle2(1:1000),delta,v01,v02);
	maxObj1000(i)=max(objVal1000(:,i));sumObj1000(i)=sum(objVal1000(:,i))/1000;
    
    
	[fct500(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:500),velDataNoise452(1:500),cAngle1(1:500),cAngle2(1:500),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct500(i)=max(fct500(:,i));sumFct500(i)=sum(fct500(:,i))/500;
	
	[xRes500(i,:),objVal500(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:500),velDataNoise452(1:500),cAngle1(1:500),cAngle2(1:500),delta,v01,v02);
	maxObj500(i)=max(objVal500(:,i));sumObj500(i)=sum(objVal500(:,i))/500;
    
    [fct250(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:250),velDataNoise452(1:250),cAngle1(1:250),cAngle2(1:250),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[2500 2500 2500 2500 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct250(i)=max(fct250(:,i));sumFct250(i)=sum(fct250(:,i))/250;
	
	[xRes250(i,:),objVal250(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[2500 2500 2500 2500 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:250),velDataNoise452(1:250),cAngle1(1:250),cAngle2(1:250),delta,v01,v02);
	maxObj250(i)=max(objVal250(:,i));sumObj250(i)=sum(objVal250(:,i))/250;
	    
    
    [fct100(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:100),velDataNoise452(1:100),cAngle1(1:100),cAngle2(1:100),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct100(i)=max(fct100(:,i));sumFct100(i)=sum(fct100(:,i))/100;
	
	[xRes100(i,:),objVal100(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:100),velDataNoise452(1:100),cAngle1(1:100),cAngle2(1:100),delta,v01,v02);
	maxObj100(i)=max(objVal100(:,i));sumObj100(i)=sum(objVal100(:,i))/100;

    
	[fct50(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:50),velDataNoise452(1:50),cAngle1(1:50),cAngle2(1:50),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct50(i)=max(fct50(:,i));sumFct50(i)=sum(fct50(:,i))/50;
	
	[xRes50(i,:),objVal50(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:50),velDataNoise452(1:50),cAngle1(1:50),cAngle2(1:50),delta,v01,v02);
	maxObj50(i)=max(objVal50(:,i));sumObj50(i)=sum(objVal50(:,i))/50;
    
	[fct25(:,i)]=brBothGrowingOptimFct(xObj,velDataNoise451(1:25),velDataNoise452(1:25),cAngle1(1:25),cAngle2(1:25),delta);
	%[xRes(i,:),objVal(i,:)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
	maxFct25(i)=max(fct25(:,i));sumFct25(i)=sum(fct25(:,i))/25;
	
	[xRes25(i,:),objVal25(:,i)]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451(1:25),velDataNoise452(1:25),cAngle1(1:25),cAngle2(1:25),delta,v01,v02);
	maxObj25(i)=max(objVal25(:,i));sumObj25(i)=sum(objVal25(:,i))/25;
    
%     a(1,:)=[mean(maxFct3000(find(maxFct3000<2000))) std(maxFct3000(find(maxFct3000<2000)))];
    a(1,:)=[mean(maxFct3000(find(maxFct3000<2000))) std(maxFct1000(find(maxFct3000<2000)))];
  
    a(2,:)=[mean(maxFct2000(find(maxFct2000<2000))) std(maxFct500(find(maxFct2000<2000)))];
    
    a(3,:)=[mean(maxFct1000(find(maxFct1000<2000))) std(maxFct1000(find(maxFct1000<2000)))];
  
    a(4,:)=[mean(maxFct500(find(maxFct500<2000))) std(maxFct500(find(maxFct500<2000)))];
  
    a(5,:)=[mean(maxFct250(find(maxFct250<2000))) std(maxFct250(find(maxFct250<2000)))];
  
    a(6,:)=[mean(maxFct100(find(maxFct100<2000))) std(maxFct100(find(maxFct100<2000)))];
  
    a(7,:)=[mean(maxFct50(find(maxFct50<2000))) std(maxFct50(find(maxFct50<2000)))];
  
    a(8,:)=[mean(maxFct25(find(maxFct25<2000))) std(maxFct25(find(maxFct25<2000)))];
  
    maxFct=a;
    
    a(1,:)=[mean(maxObj3000(find(maxObj3000<2000))) std(maxObj3000(find(maxObj3000<2000)))];
    
    a(2,:)=[mean(maxObj2000(find(maxObj2000<2000))) std(maxObj2000(find(maxObj2000<2000)))];
  
    a(3,:)=[mean(maxObj1000(find(maxObj1000<2000))) std(maxObj1000(find(maxObj1000<2000)))];
  
    a(4,:)=[mean(maxObj500(find(maxObj500<2000))) std(maxObj500(find(maxObj500<2000)))];
  
    a(5,:)=[mean(maxObj250(find(maxObj250<2000))) std(maxObj250(find(maxObj250<2000)))];
  
    a(6,:)=[mean(maxObj100(find(maxObj100<2000))) std(maxObj100(find(maxObj100<2000)))];
  
    a(7,:)=[mean(maxObj50(find(maxObj50<2000))) std(maxObj50(find(maxObj50<2000)))];
  
    a(8,:)=[mean(maxObj25(find(maxObj25<2000))) std(maxObj25(find(maxObj25<2000)))];
  
    maxObj=a;
  
    a(1,:)=[mean(sumObj3000(find(sumObj3000<2000))) std(sumObj3000(find(sumObj3000<2000)))];
    
    a(2,:)=[mean(sumObj2000(find(sumObj2000<2000))) std(sumObj2000(find(sumObj2000<2000)))];
    
    a(3,:)=[mean(sumObj1000(find(sumObj1000<2000))) std(sumObj1000(find(sumObj1000<2000)))];
 
    a(4,:)=[mean(sumObj500(find(sumObj500<2000))) std(sumObj500(find(sumObj500<2000)))];
  
    a(5,:)=[mean(sumObj250(find(sumObj250<2000))) std(sumObj250(find(sumObj250<2000)))];
  
    a(6,:)=[mean(sumObj100(find(sumObj100<2000))) std(sumObj100(find(sumObj100<2000)))];
  
    a(7,:)=[mean(sumObj50(find(sumObj50<2000))) std(sumObj50(find(sumObj50<2000)))];
  
    a(8,:)=[mean(sumObj25(find(sumObj25<2000))) std(sumObj25(find(sumObj25<2000)))];
  
    sumObj=a;
    
  
    a(1,:)=[mean(sumFct3000(find(sumFct3000<2000))) std(sumFct3000(find(sumFct3000<2000)))];
    
    a(2,:)=[mean(sumFct2000(find(sumFct2000<2000))) std(sumFct2000(find(sumFct2000<2000)))];
    
    a(3,:)=[mean(sumFct1000(find(sumFct1000<2000))) std(sumFct1000(find(sumFct1000<2000)))];
  
    a(4,:)=[mean(sumFct500(find(sumFct500<2000))) std(sumFct500(find(sumFct500<2000)))];
  
    a(5,:)=[mean(sumFct250(find(sumFct250<2000))) std(sumFct250(find(sumFct250<2000)))];
  
    a(6,:)=[mean(sumFct100(find(sumFct100<2000))) std(sumFct100(find(sumFct100<2000)))];
  
    a(7,:)=[mean(sumFct50(find(sumFct50<2000))) std(sumFct50(find(sumFct50<2000)))];
  
    a(8,:)=[mean(sumFct25(find(sumFct25<2000))) std(sumFct25(find(sumFct25<2000)))];

    sumFct=a;
	
% 	% Plot section
% 	
% 	xObj=[25 0.15 45 0.3];
% 	

% 	hold on;
% 	
% 	omega=[0:0.05:5];
% 	plot(omega1,velDataNoise4551,'dr');
% 	plot(omega2,velDataNoise452,'dk');
% 	
% 	plot(omega,peskin(xObj(1:2),omega,delta),'-.r');
% 	plot(omega,peskin(xObj(3:4),omega,delta),'-.k');
% 	plot(omega,peskin(xResN2(i,1:2),omega,delta),'-r');
% 	plot(omega,peskin(xResN2(i,3:4),omega,delta),'-k');
% 	
% 	text(3,1e-3,['sum of objVal = ' num2str(sum(objValN2(:,i))/length(velData1))]);
% 	
% 	text(3,1.5e-3,['max of objVal = ' num2str(max(objValN2(:,i)))]);
% 	a=num2str((xObj(1) -xResN2(i,1))./xObj(1));
% 	b= num2str((xObj(2) -xResN2(i,2))./xObj(2));
% 	c=num2str((xObj(3) -xResN2(i,3))./xObj(3));
% 	d=num2str((xObj(4) -xResN2(i,4))./xObj(4));
% 	text(2,3e-3,['rel error: {' a ';' b ';' c ';' d '}' ]);
% 	
% 	velData1=velData1(1:25);
% 	velData2=velData2(1:25);
% 	velData1=velDataNoise1(1:25);
% 	velDataNoise452=velDataNoise452(1:25);
% 
% 	
% 	figure
% 	hold on;

% 	
% 
% 	[xResN3,objValN3]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1(1:25),cAngle2(1:25),delta,v01,v02);
% 	omega=[0:0.05:5];
% 	plot(omega1(1:25),velDataNoise451,'dr');
% 	plot(omega2(1:25),velDataNoise452,'dk');
% 	
% 	plot(omega,peskin(xRes(1:2),omega,delta),'-.r');
% 	plot(omega,peskin(xRes(3:4),omega,delta),'-.k');
% 	plot(omega,peskin(xResN3(1:2),omega,delta),'-r');
% 	plot(omega,peskin(xResN3(3:4),omega,delta),'-k');
% 	
% 	text(3,1e-3,['sum of objVal = ' num2str(sum(objValN3)/length(velData1))]);
% 	
% 	text(3,1.5e-3,['max of objVal = ' num2str(max(objValN3))]);
% 	a=num2str((xObj(1) -xResN3(1))./xObj(1));
% 	b= num2str((xObj(2) -xResN3(2))./xObj(2));
% 	c=num2str((xObj(3) -xResN3(3))./xObj(3));
% 	d=num2str((xObj(4) -xResN3(4))./xObj(4));
% 	text(2,3e-3,['rel error: {' a ';' b ';' c ';' d '}' ]);
% 	
%     %----------------- 5%
% 	[xResN52,objValN52]=fminimax(@brBothGrowingOptimFct45,xObj,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v01,v02);
% 	
% 	
% 	% Plot section
% 	
% 	xObj=[25 0.15 45 0.3];
% 	figure
% 	hold on;
% 	
% 	omega=[0:0.05:5];
% 	plot(omega1,velDataNoise451,'dr');
% 	plot(omega2,velDataNoise52,'dk');
% 	
% 	plot(omega,peskin(xRes(1:2),omega,delta),'-.r');
% 	plot(omega,peskin(xRes(3:4),omega,delta),'-.k');
% 	plot(omega,peskin(xResN52(1:2),omega,delta),'-r');
% 	plot(omega,peskin(xResN52(3:4),omega,delta),'-k');
% 	
% 	text(3,1e-3,['sum of objVal = ' num2str(sum(objValN52)/length(velData1))]);
% 	
% 	text(3,1.5e-3,['max of objVal = ' num2str(max(objValN52))]);
% 	a=num2str((xObj(1) -xResN52(1))./xObj(1));
% 	b= num2str((xObj(2) -xResN52(2))./xObj(2));
% 	c=num2str((xObj(3) -xResN52(3))./xObj(3));
% 	d=num2str((xObj(4) -xResN52(4))./xObj(4));
% 	text(2,3e-3,['rel error: {' a ';' b ';' c ';' d '}' ]);
% 	
% 	velData1=velData1(1:25);
% 	velData2=velData2(1:25);
% 	velDataNoise51=velDataNoise51(1:25);
% 	velDataNoise52=velDataNoise52(1:25);
% 	cAngle1=cAngle1(1:25);
% 	cAngle2=cAngle2(1:25);
% 	
% 	figure
% 	hold on;
% 	
% 	omega1=omega1(1:25);
% 	omega2=omega2(1:25);
% 	[xResN53,objValN53]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst45,opts,velDataNoise51,velDataNoise52,cAngle1,cAngle2,delta,v01,v02);
% 	omega=[0:0.05:5];
% 	plot(omega1,velDataNoise51,'dr');
% 	plot(omega2,velDataNoise52,'dk');
% 	
% 	plot(omega,peskin(xRes(1:2),omega,delta),'-.r');
% 	plot(omega,peskin(xRes(3:4),omega,delta),'-.k');
% 	plot(omega,peskin(xResN53(1:2),omega,delta),'-r');
% 	plot(omega,peskin(xResN53(3:4),omega,delta),'-k');
% 	
% 	text(3,1e-3,['sum of objVal = ' num2str(sum(objValN53)/length(velData1))]);
% 	
% 	text(3,1.5e-3,['max of objVal = ' num2str(max(objValN53))]);
% 	a=num2str((xObj(1) -xResN53(1))./xObj(1));
% 	b= num2str((xObj(2) -xResN53(2))./xObj(2));
% 	c=num2str((xObj(3) -xResN53(3))./xObj(3));
% 	d=num2str((xObj(4) -xResN53(4))./xObj(4));
% 	text(2,3e-3,['rel error: {' a ';' b ';' c ';' d '}' ]);
	
    
end

