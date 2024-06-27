% This script will run cwt analysis on a number of artificicial
% data sets of increasing complexity.
% These results will be used to validate this analysis in a paper,
% before using the package on CR WL data.

addpath('../cwtMultiPackage/CWT_Multi')

%% add custom filter lengths for the constituent analysis
% lon, lat for Astoria
lon = 123+46.1/60;            % W
lat = 46+12.4/60;             % N

%% including K2
% coFiltLength    = [1081,1081,337,1081,337,1081,1081,1081,1081,337,337,1081,1081,337,1081,1081,1081,1081,1081,1081,1081,1081,1081,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours
% 
% constitsD1Freqs = [28.00621204, 26.868350, 25.81933871, 24.84120241, 23.93447213, 23.09848146, 22.30608083].^-1;
% constitsD1Names = ["TwoQ1","Q1","O1","NO1","K1","J1","OO1"];
% 
% constitsD2Freqs = [13.1272668, 12.8717576, 12.65834751, 12.4206012, 12.19162085, 12, 1/0.0835614924, 11.75452172].^-1;
% constitsD2Names = ["Eps2","Mu2","N2","M2","L2","S2","K2","Eta2"];   % Nu2 is number 10
% 
% constitsD3Freqs = [8.3863068, 8.280400802, 8.177140247, 7.9936].^-1;
% constitsD3Names = ["MO3","M3","MK3","SK3"];
% 
% constitsD4PlusFreqs = [6.2691737, 6.21030, 6.103339, 6.0, 4.93088021306, 4.140200399, 3.52964079728,...
%     3.10515029954, 2.748563985947, 2.48412023963, 2.25054027184, 2.07010019969].^-1;
% constitsD4PlusNames = ["MN4","M4","MS4","S4","MK5","M6","MK7","M8","MK9","M10","MK11","M12"];
% 
% coFreqs = {constitsD1Freqs,constitsD2Freqs,constitsD3Freqs,constitsD4PlusFreqs};
% coNames = horzcat(constitsD1Names,constitsD2Names,constitsD3Names,constitsD4PlusNames);

%% run cwt routine on artficial data
load('./data/Artificial/artifS_6min.mat')
% load('./data/Artificial/M2NSsquare.mat')

dat.dtime.TimeZone = 'UTC';

%%
% coFiltLength = [363,362,362,362,362,362,363,363,363,363,363,362,362,363,362,363,363,362,362,362,362,362,362,362,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours
% coFiltLength = [1081,1080,1080,1081,1080,1081,1081,1081,1081,1080,1080,1081,1080,1081,1080,1080,1080,1080,1080,1080,1080,1080,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours
% coNames = ["TwoQ1","Q1","O1","NO1","K1","J1","OO1","Eps2","Mu2","N2","M2","L2","S2","Eta2","MO3","M3","MK3","SK3","MN4","M4","MS4","S4"];

% spFiltLength = [336,181,145,73,65,65,65,65];
% spNames = ["low4d","low2d","d1","d2","d3","d4","d6","d8"];

% coFiltLength = 1080*ones(22,1);

% [k2constits,k2species,k2ref,k2admit] = cwtMulti(dat.dtime,dat.wl,'coFreqs',coFreqs,coFiltLength,coNames);
% %     'pfAmps','pfResid','pfResidSpectra','pfResp');

% [constits,species,ref,admit] = cwt_multi(dat.dtime,dat.wl,lon,lat);

[constits,species,ref,admit] = cwtMulti(dat.dtime,dat.wlnoise); % ,'coFiltLength',coFiltLength);

% load('./data/Artificial/artifSNoT2R2.mat')
% 
% [modk2constits,~,~,~] = cwtMulti(dat.dtime,dat.wl,'coFreqs',coFreqs,coFiltLength,coNames);

%% time stuff
dat.dtime.TimeZone='UTC';
constits.decTimesAll.TimeZone = 'UTC';
%% run UTide on artificial data
% addpath('../utide_tool/')
% 
% dat.datenums = datenum(dat.dtime);
% dat.dates = dat.dtime;
% dat.lat = lat; dat.lon = lon;
% 
% % Define moving parameters etc.
% window = 30*1.0;         % window in days
% incr = 7;%1/24;     % increment in days 
% rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
% good_pct = 0.8;     % specify what percentage of data in window must be not NaN
% 
% maxFLength = 5*24;   % same cutoff period as low-pass filtering done in CWT\_Multi
% 
% dat.dates = dat.dates.' ; dat.wl = dat.wl.' ; dat.datenums = dat.datenums.'; dat.dtime = dat.dtime.';
% [cd,OUT] = cwt_utide(dat,window,incr,good_pct,rayleigh,maxFLength,[dat.dtime(1) dat.dtime(end)]);%16P1/17K1

%% Plots
% pH = ones(size(dat.t));
% 
% p=tiledlayout(4,1);
% 
% ax1=nexttile;
% plot(dat.dtime,pH*dat.amps(1),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 true')
% hold(ax1,'on')
% plot(constits.decTimesAll,constits.M2.amps,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(2),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 true')
% plot(constits.decTimesAll,constits.S2.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(3),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 true')
% plot(constits.decTimesAll,constits.N2.amps,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% hold(ax1,'off')
% grid on
% axis tight
% ylim([0 1.4])
% ylabel('Amplitude')
% % xlabel('Arbitrary time')
% set(gca,'xtick',[])
% legend('location','eastoutside')
% title('Semidiurnal amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% ax = gca;
% ax.FontSize = 14;
% plotAnn('(a)')
% 
% ax2=nexttile;
% plot(dat.dtime,pH*dat.amps(7),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 true')
% hold(ax2,'on')
% plot(constits.decTimesAll,constits.Q1.amps,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(6),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 true')
% plot(constits.decTimesAll,constits.O1.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(4),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 true')
% plot(constits.decTimesAll,constits.K1.amps,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% hold(ax2,'off')
% grid on
% axis tight
% ylim([0 0.85])
% ylabel('Amplitude')
% % xlabel('Arbitrary time')
% set(gca,'xtick',[])
% legend('location','eastoutside')
% title('Diurnal amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% ax = gca;
% ax.FontSize = 14;
% plotAnn('(b)')
% 
% ax3=nexttile;
% plot(dat.dtime,pH*dat.phases(1)*(180/pi),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 true')
% hold(ax3,'on')
% plot(constits.decTimesAll,constits.M2.phases,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(2)*(180/pi),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 true')
% plot(constits.decTimesAll,constits.S2.phases,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(3)*(180/pi),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 true')
% plot(constits.decTimesAll,constits.N2.phases,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% hold(ax3,'off')
% grid on
% axis tight
% ylim([0 370])
% ylabel('Phase (rad)')
% % xlabel('Arbitrary time')
% set(gca,'xtick',[])
% legend('location','eastoutside')
% title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% % title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 0.8, 0])
% ax = gca;
% ax.FontSize = 14;
% plotAnn('(c)')
% 
% ax4=nexttile();
% plot(dat.dtime,pH*dat.phases(7)*(180/pi),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 true')
% hold on
% plot(constits.decTimesAll,unwrap(constits.Q1.phases),'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(6)*(180/pi),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 true')
% plot(constits.decTimesAll,constits.O1.phases,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(4)*(180/pi),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 true')
% plot(constits.decTimesAll,constits.K1.phases,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% hold off
% grid on
% axis tight
% ylim([0 390])
% ylabel('Phase (rad)')
% % xlabel('Arbitrary time')
% legend('location','eastoutside')
% title('Diurnal phases', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% ax = gca;
% ax.FontSize = 14;
% plotAnn('(d)')
% 
% set(gcf,'Position',[100 100 2500 1200])
% % saveas(p,'./figs/artificial_stationary_all_amps.png')



%% k2 s2 comparison

p = figure();
plot(dat.dtime,pH*dat.amps(2),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 true')
hold on
plot(constits.decTimesAll,constits.S2.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 default filters')
plot(constits.decTimesAll,k2constits.S2.amps,'-','linewidth',2,'DisplayName','S_2 w/ K_2 filter')
plot(constits.decTimesAll,modk2constits.S2.amps,'-','linewidth',2,'DisplayName','S_2 w/ K_2 filter (no T_2 or R_2)')
hold off
grid on
axis tight
ylim([0.3,0.65]) % ylim([0.3 0.625])
ylabel('Amplitude')
% xlabel('Arbitrary time')
legend('location','northeast')
title('Artificial data S2 amplitude case study', 'Units', 'normalized', 'Position', [0.5, 1.025, 0])
ax = gca;
ax.FontSize = 14;
set(gcf,'Position',[100 100 1800 500])
% saveas(p,'./figs/artificial_s2_case_study.png')



%% Alternate main panel w/ sums
load('./data/Artificial/artifS.mat')
dat.dtime.TimeZone = 'UTC';
constits.decTimesAll.TimeZone = 'UTC';

figure();
pH = ones(size(dat.t));

p=tiledlayout(4,1);

ax1=nexttile;
plot(dat.dtime,pH*dat.amps(1),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 true')
hold(ax1,'on')
plot(constits.decTimesAll,constits.M2.amps,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(2),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 true')
full_S2 = abs(dat.S2K2compWL);
plot(dat.dtime,full_S2,'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 + K_2 true')
plot(constits.decTimesAll,constits.S2.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(3),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 true')
full_N2 = abs(dat.N2NU2compWL);
plot(dat.dtime,full_N2,'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 + NU_2 true')
plot(constits.decTimesAll,constits.N2.amps,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 CWT\_Multi')
hold(ax1,'off')
grid on
axis tight
ylim([0 1.4])
ylabel('Amplitude')
% xlabel('Arbitrary time')
set(gca,'xtick',[])
legend('location','eastoutside')
title('Semidiurnal amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(dat.dtime,pH*dat.amps(7),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 true')
full_Q1 = abs(dat.Q1RHO1compWL);
plot(dat.dtime,full_Q1,'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 + RHO_1 true')
hold(ax2,'on')
plot(constits.decTimesAll,constits.Q1.amps,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plot(dat.dtime,pH*dat.amps(6),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 true')
plot(constits.decTimesAll,constits.O1.amps,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.amps(4),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 true')
full_K1 = abs(dat.K1P1compWL);
plot(dat.dtime,full_K1,'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 + P_1 true')
plot(constits.decTimesAll,constits.K1.amps,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% full_K1 = pH*(dat.amps(4) + dat.amps(26) + dat.amps(22) + dat.amps(12) + dat.amps(23) );
hold(ax2,'off')
grid on
axis tight
ylim([0 0.85])
ylabel('Amplitude')
% xlabel('Arbitrary time')
set(gca,'xtick',[])
legend('location','eastoutside')
title('Diurnal amplitudes', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(dat.dtime,pH*dat.phases(1)*(180/pi),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 true')
hold(ax3,'on')
plot(constits.decTimesAll,constits.M2.phases,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(2)*(180/pi),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 true')
full_S2_phase = mod(2*pi*dat.freqs(2)*dat.t-angle(dat.S2K2compWL),2*pi)*(180/pi); 
plot(dat.dtime,full_S2_phase,'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 + K_2 true')
plot(constits.decTimesAll,constits.S2.phases,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(3)*(180/pi),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 true')
full_N2_phase = mod(2*pi*dat.freqs(3)*dat.t-angle(dat.N2NU2compWL),2*pi)*(180/pi); 
plot(dat.dtime,full_N2_phase,'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 + NU_2 true')
plot(constits.decTimesAll,constits.N2.phases,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','N_2 CWT\_Multi')
hold(ax3,'off')
grid on
axis tight
ylim([0 370])
ylabel('Phase (rad)')
% xlabel('Arbitrary time')
set(gca,'xtick',[])
legend('location','eastoutside')
title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
% title('Semidiurnal phases', 'Units', 'normalized', 'Position', [0.5, 0.8, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile();
% plot(dat.dtime,pH*dat.phases(7)*(180/pi),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 true')
full_Q1_phase = mod(2*pi*dat.freqs(7)*dat.t-angle(dat.Q1RHO1compWL),2*pi)*(180/pi); 
plot(dat.dtime,full_Q1_phase,'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 + RHO_1 true')
hold on
plot(constits.decTimesAll,constits.Q1.phases,'-','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plot(dat.dtime,pH*dat.phases(6)*(180/pi),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 true')
plot(constits.decTimesAll,constits.O1.phases,'-','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plot(dat.dtime,pH*dat.phases(4)*(180/pi),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 true')
full_K1_phase = mod(2*pi*dat.freqs(4)*dat.t-angle(dat.K1P1compWL),2*pi)*(180/pi); 
plot(dat.dtime,full_K1_phase,'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 + P_1 true')
plot(constits.decTimesAll,constits.K1.phases,'-','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','K_1 CWT\_Multi')

hold off
grid on
axis tight
ylim([0 390])
ylabel('Phase (rad)')
% xlabel('Arbitrary time')
legend('location','eastoutside')
title('Diurnal phases', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

set(gcf,'Position',[100 100 2500 1200])
% saveas(p,'./figs/artificial_stationary_all_amps.png')


%% %%%%%%%%%%%%%% PLOTTING FUNCTION
function plotAnn(txt)
%     xlim(xlim+[-10,10]);ylim(ylim+[-10,10]);
    text(min(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end


