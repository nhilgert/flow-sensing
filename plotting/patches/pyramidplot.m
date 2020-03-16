clear all

load('pyramid')
fsize = 18;
imagesc(cdata);
set(gca, 'XTick', [0:306:1533]-125, 'XTickLabel', [-1:10] , 'FontSize', fsize) % 10 ticks
set(gca, 'YTick', [1:1563/9:1563]+70, 'YTickLabel', [4:-1:-4])
set(gca, 'TickLength', [0 0]);

c = colorbar('Ticks',0:0.25:1,'TickLabels',{'0','0.5','1.0','1.5','2.0'});

fsize = 30;
xlabel('$\ell$','interpreter','latex','FontSize',fsize);
ylabel('$m$','interpreter','latex','FontSize',fsize);

set(gcf, 'Position', [10 10 800 600]);

saveas(gcf, 'pyramidplot.png')