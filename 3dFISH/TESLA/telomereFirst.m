%% 2016/10/26, telomere first
% From Jungsik

load U2OSlanesstatistics3.mat

datamat = U2OSlanesstatistics3;

size(datamat)
sum(sum(~isnan(datamat(:))))

vec = datamat(~isnan(datamat));
length(unique(vec))

% descriptive
figure
boxplot(datamat)

summary(datamat(:))

%%

datamat(datamat > 10) = NaN;

summary(datamat(:))



%% histo at the resolution of 1kb
%discretesupport = 0:0.1:17.6;

datDiscrete1 = round(datamat, 0);

mov(1:33)= struct('cdata',[],'colormap',[]);
figure

for k = 1:32
    
    h = histogram(datamat(:,k));
    h.Normalization = 'count';
    h.BinEdges = -0.5:1:17.5;
    %h.BinEdges = [0:9, 10, 18];
    mov(k)=getframe(gcf);
end

implay(mov)


%% bootstrapping & support coverage rate

totnumLane = 32
numBoots = 500
bootstrappedCoverage = nan(numBoots, totnumLane);

% total sample support
denSupport = unique(datDiscrete1(~isnan(datDiscrete1)));

for numReplicate = 1:totnumLane

    for b = 1:numBoots
        bind = randsample(totnumLane, numReplicate, true);
        bind1 = sort(bind);
        x = datDiscrete1(:, bind1);
        denS1 = unique(x(~isnan(x)));
        p1 = numel(denS1)/numel(denSupport);
        bootstrappedCoverage(b, numReplicate) = p1;
    end

end

bootCoverageRate = mean(bootstrappedCoverage);

% Confidence Bound at 95% confidence level 
uCL = quantile(bootstrappedCoverage, 0.975);
lCL = quantile(bootstrappedCoverage, 0.025);


pmat = [bootCoverageRate', uCL', lCL'];

figure, 
plot(pmat)
legend('avg. Coverage', 'upper CL', 'lower CL', 'Location', 'southeast') 

h = refline([0, 0.9]);
h.Color = [.5 .5 .5];
h.LineStyle = '--';
h = refline([0, 0.8]);
h.Color = [.5 .5 .5];
h.LineStyle = '--';




%%  end
