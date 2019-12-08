clear all

load('simdat');

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 400])

%% means vs variances plot
subplot(1,2,2);
hold on
for i = 1:5
    counts = squeeze(simcounts(i,:,:));
    N = size(counts,1);
    variances = var(counts);
    means = mean(counts);
    plot(means,variances,'*');
end

legendsize = 16;
titlesize = 25;

xlim([0 200]);
ylim([0 200]);
xlabel('Mean captured molecules, $\bar{n}$','Interpreter','LaTeX','FontSize',titlesize);
ylabel('Variance of captured molecules, $\sigma_n^2$','Interpreter','LaTeX','FontSize',titlesize);
plot(0:200,0:200,'--','LineWidth',2) % y = x reference
legend({'$\epsilon$ = 0.1', '$\epsilon$ = 1', '$\epsilon$ = 10', '$\epsilon$ = 30', '$\epsilon$ = 100'},...
    'Interpreter','LaTeX','Location', 'NorthWest','FontSize',legendsize);

%% histogram
subplot(1,2,1);
h = squeeze(simcounts(2,:,75));
hold on
histogram(h,'Normalization','pdf')
orange = [0.8500 0.3250 0.0980]; % rgb triple for matlab default orange
plot(100:200,poisspdf(100:200,mean(h)),'o','Color',orange,'MarkerSize',4,'MarkerFaceColor',orange)
pdist = plot(100:200,poisspdf(100:200,mean(h)),'LineWidth',2,'Color',orange)

xlabel('Number of molecules captured, $n$','Interpreter','LaTeX','FontSize',titlesize);
ylabel('P($n$)','Interpreter','LaTeX','FontSize',titlesize);
legend(pdist,{'Poiss($\bar{n}$)'},'Interpreter','LaTeX','FontSize',legendsize)

print(fig,'simplots','-dpng')