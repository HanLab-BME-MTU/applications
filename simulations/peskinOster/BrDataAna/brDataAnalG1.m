% -----------------brDataAnal

clear;close all;

load('velocityG1.mat');

velS=abs(vel(find(vel<0)));
velSCons=abs(velCons(find(velCons<0)));


vUnload=mean(velS);
stdVel=std(velS);


discret=[min(velS):1e-4:max(velS)];
for i=1:length(discret)
    a(i)=length(find(velS<=discret(i)));
end

plot(discret,a);

hold on;

discretCons=[min(velSCons):1e-4:max(velSCons)];
for i=1:length(discretCons)
    aCons(i)=length(find(velSCons<=discretCons(i)));
end

plot(discretCons,aCons,'r');





