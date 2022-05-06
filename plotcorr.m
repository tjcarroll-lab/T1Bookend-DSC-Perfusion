function fh = plotcorr(x,y,fh)
%plotcorr Plot a correlation plot
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2020-11-02

if nargin < 3
    fh = figure;
    %plotcolor = [0 0.447 0.741];
    plotcolor = [0 0 0];
else
    plotcolor = [0.85 0.325 0.098];
end
if size(x,1) < size(x,2)
    x = x';
end
if size(y,1) < size(y,2)
    y = y';
end
figure(fh);
hold on;
plot(x,y,'o','markersize',6,'markerfacecolor',plotcolor,'markeredgecolor',plotcolor);
hold on;
[b,bint,r,rint,stats] = regress(y,[x,ones(size(x,1),1)]);
%[b,bint,r,rint,stats] = regress(y,x);bint
%b(2) = 0;
plot(get(gca,'xlim'),b(1).*get(gca,'xlim')+b(2),'--','color',plotcolor,'linewidth',2);
% plot(get(gca,'xlim'),bint(1,1).*get(gca,'xlim')+bint(2,1),':','color',color1,'linewidth',2);
% plot(get(gca,'xlim'),bint(1,2).*get(gca,'xlim')+bint(2,2),':','color',color1,'linewidth',2);
% if plotsimple
%     leg1 = sprintf('not DD corrected');
% else
%     leg1 = sprintf('%s: %.3g*x + %.3g; r^2 = %.3g; p = %.3g','no DD',b(1),b(2),stats(1),stats(3));
% end
set(gca,'FontSize',12);
if stats(3) < 0.01
    pvalue = '< 0.01';
else
    pvalue = ['= ' sprintf('%.2f',stats(3))];
end
CCC = f_CCC([x y],0.05);
leg1 = {sprintf('y = %.2f*x + %.2f',b(1),b(2)),...
    sprintf('r = %.2f; p-value %s',sqrt(stats(1)),pvalue),... % Do pearson r (not r^2)
    %sprintf('\\rho_c = %.2f',CCC.est),... % Lin's CCC
    };
txtbox = annotation('textbox',[.35 .55 .3 .3],'string',leg1,'FitBoxToText','on','linestyle','none','fontsize',14,'horizontalalignment','center','fontangle','italic','color',plotcolor,'fontweight','bold');
end

