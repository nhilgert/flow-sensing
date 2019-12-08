clear all
clc

fig = figure('Renderer', 'painters', 'Position', [10 10 1125 450]);

b1000err = dlmread('./errors/b1000err.dat');
epsilon = 10.^(-2:0.02:0);
epsilon = epsilon(1:end-1);

subplot(1,2,2)

loglog(epsilon,b1000err,'LineWidth',1.5);
xlim([0.01 1])

line([.1 .1], [.000001 3.344e-3],'LineWidth',1,'Color','k','LineStyle','--');
line([.001 .1], [3.344e-3 3.344e-3],'LineWidth',1,'Color','k','LineStyle','--');

% fontsizes
titlesize = 25;
legendsize = 15;
ax = gca;
ax.FontSize = 14;
xlabel("P\'eclet number, $\epsilon$",'Interpreter','LaTeX','FontSize',titlesize);
ylabel('Approximation error','Interpreter','LaTeX','FontSize',titlesize);


subplot(1,2,1)
hold on
chi = dlmread('./perturbed51.dat');

numerical = dlmread('./numericalsolns/at51/b1000pec51.dat');

plot([0:0.001:pi],numerical,'LineWidth',1)
plot([0:0.001:pi],chi,'LineWidth',2.5,'LineStyle','--')

legend({'Numerical', 'Perturbative'},...
    'Interpreter','LaTeX','FontSize',legendsize,'Location','NorthEast')

titlesize = 25;
legendsize = 15;
ax = gca;
ax.FontSize = 14;
xlabel('$\theta$','Interpreter','LaTeX','FontSize',titlesize);
ylabel("Relative surface concentration",'Interpreter','LaTeX','FontSize',titlesize);
xlim([0,pi]);
xticks([0:pi/4:pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})

print(fig,'pdeplot','-dpng')
