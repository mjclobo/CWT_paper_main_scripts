% This script will run cwt analysis on a number of artificicial
% data sets of increasing complexity.
% These results will be used to validate this analysis in a paper,
% before using the package on CR WL data.

addpath('../cwtMultiPackage/CWT_Multi')

%% add custom filter lengths for the constituent analysis
% lon, lat for Astoria
lon = 123+46.1/60;            % W
lat = 46+12.4/60;             % N

%% run cwt routine on artficial data
load('./data/Artificial/artifNS.mat')
% load('./data/Artificial/M2NSsquare.mat')

coFiltLength = [1081,1080,1080,1081,1080,1081,1081,1081,1081,1080,1080,1081,1080,1081,1080,1080,1080,1080,1080,1080,1080,1080,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours
coNames = ["TwoQ1","Q1","O1","NO1","K1","J1","OO1","Eps2","Mu2","N2","M2","L2","S2","Eta2","MO3","M3","MK3","SK3","MN4","M4","MS4","S4"];

% spFiltLength = [336,181,145,73,65,65,65,65];
% spNames = ["low4d","low2d","D1","D2","D3","D4","d6","d8"];

% coFiltLength = 1080*ones(22,1);

[constits,species,ref,admit] = cwtMulti(dat.dtime,dat.wl,'rel_phase','1970-01-01'); %,'coFiltLength',coFiltLength);

%% Plots
on = ones(length(constits.decTimesAll),1);
ph = ones(length(dat.dtime),1);

figure(); p=tiledlayout(3,1);

ax1 = nexttile();
plot(dat.dtime,dat.M2.amp,'--','color',[0 0.4470 0.7410],'linewidth',2.0,'DisplayName','M_2 true')
hold(ax1,'on')
plot(constits.decTimesAll,constits.M2.amps,'-','color',[0 0.4470 0.7410],'linewidth',2.0,'DisplayName','M_2 CWT\_Multi')
plot(constits.decTimesAll,on*(17386/90809),'--','color',[0.8500 0.3250 0.0980],'linewidth',2.0,'DisplayName','N_2 True')
plot(constits.decTimesAll,constits.N2.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2.0,'DisplayName','N_2 CWT\_Multi')
plot(constits.decTimesAll,on*(42248/90809),'--','color',[0.9290 0.6940 0.1250],'linewidth',2.0,'DisplayName','S_2 True')
plot(constits.decTimesAll,constits.S2.amps,'-','color',[0.9290 0.6940 0.1250],'linewidth',2.0,'DisplayName','S_2 CWT\_Multi')
hold(ax1,'off')
set(gca,'xtick',[])
grid on
xlim([dat.dtime(1) dat.dtime(end)])
ylim([0 1.2])
ylabel('Amplitude')
legend('location','southeast')
title('Semidiurnal constituent amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2 = nexttile();
plot(dat.dtime,dat.D2.amps,'--','color',[0 0.4470 0.7410],'linewidth',2.0,'DisplayName','D_2 true')
hold(ax2,'on')
plot(species.decTimesAll,species.D2.amps,'-','color',[0 0.4470 0.7410],'linewidth',2.0,'DisplayName','D_2 CWT\_Multi (species)')
hold(ax2,'off')
set(gca,'xtick',[])
grid on
xlim([dat.dtime(1) dat.dtime(end)])
ylabel('Amplitude')
legend('location','northeast')
title('Semidiurnal species amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')


ax3=nexttile;
plot(dat.dtime,ph*dat.phases(1)*(180/pi),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 true')
hold(ax3,'on')
plot(constits.decTimesAll,constits.M2.phases,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plot(constits.decTimesAll,on*dat.phases(3)*(180/pi),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','N_2 true')
plot(constits.decTimesAll,constits.N2.phases,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plot(constits.decTimesAll,on*dat.phases(2)*(180/pi),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','S_2 true')
plot(constits.decTimesAll,constits.S2.phases,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','S_2 CWT\_Multi')
hold(ax3,'off')
grid on
axis tight
ylim([0 360])
ylabel('Phase (rad)')
% xlabel('Arbitrary time')
legend('location','southeast')
title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 0.8, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

set(gcf,'Position',[100 100 1500 800])
saveas(p,'./figs/artificial_nonstationary_D2_amps_phases.png')

%% %%%%%%%%%%%%%% PLOTTING FUNCTION
function plotAnn(txt)
%     xlim(xlim+[-10,10]);ylim(ylim+[-10,10]);
    text(min(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end



