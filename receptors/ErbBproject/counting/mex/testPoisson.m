
x=randi(100,256,256)+100;

tic
for k=1:1000
    a=poissrnd(x);
end
toc

toc-tic

tic
for k=1:1000
    a=poissonRV(x);
end
toc

toc-tic