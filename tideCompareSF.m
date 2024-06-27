% Comparing CWT_Multi to UTide
% clear
% [constits,species,ref,admit] = cwtMulti(dat4.dates,dat4.wl);
% [constits,species,ref,admit] = cwtMulti(FR_dat.dates,FR_dat.wl,'dynamicInference');
%% load data
addpath('../cwtMultiPackage/CWT_Multi')

load('./data/dayflow_1997to2020_NDOI_RIO_WEST_XGEO.mat')
dayFlow.dtime = datetime(dayFlow.dn,'ConvertFrom','datenum','Format','HH:mm:ss.SSS');

load('./FR_dat.mat')
load('./sf_data.mat')

%% align datasets to common start/end
startDate = datetime(2014,06,01); % change this back to 2013,06 to avoid edge effects
endDate = datetime(2018,06,01);  % 2018,06,01

% endDate = datetime(2018,05,31);  % 2018,06,01

startDate.TimeZone='UTC'; endDate.TimeZone='UTC';
FR_dat.dates.TimeZone='UTC';

FR_dat.start = find(FR_dat.dates==startDate);
FR_dat.end = find(FR_dat.dates==endDate);

FR_dat.dates = FR_dat.dates(FR_dat.start:FR_dat.end);
FR_dat.wl = FR_dat.wl(FR_dat.start:FR_dat.end);

sf_data.start = find(sf_data.dates==startDate);
sf_data.end = find(sf_data.dates==endDate);

sf_data.dates = sf_data.dates(sf_data.start:sf_data.end);
sf_data.wl = sf_data.wl(sf_data.start:sf_data.end);

nan_ind = unique(sort([find(isnan(FR_dat.wl));find(isnan(sf_data.wl))]));

FR_dat.wl(nan_ind) = NaN;
sf_data.wl(nan_ind) = NaN;

% set dates for plotting
end_date = datetime(2017,10,01,'TimeZone','UTC');

%% UTide analysis and reconstruction
addpath('./utide_tool/')

sf_data.datenums = datenum(sf_data.dates);
FR_dat.datenums = datenum(FR_dat.dates);
sf_data.dtime = sf_data.dates;
FR_dat.dtime = FR_dat.dates;

sf_data.dates.TimeZone='UTC';
FR_dat.dates.TimeZone='UTC';

rayleigh = 0.9;

%% a single call to UTide for whole time period
% sf_coef_all = ut_solv(sf_data.datenums, sf_data.wl,[],sf_data.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
% 
% fr_coef_all = ut_solv(FR_dat.datenums, FR_dat.wl,[],FR_dat.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
%         
% % export to csv for David
% tt = table(zeros(length(sf_coef_all.name),1));
% tt = addvars(tt,sf_coef_all.name,'NewVariableNames','name');
% tt = addvars(tt,sf_coef_all.A,'NewVariableNames','amp');
% tt = addvars(tt,sf_coef_all.A_ci,'NewVariableNames','amp_CI');
% tt = addvars(tt,sf_coef_all.g,'NewVariableNames','phase');
% tt = addvars(tt,sf_coef_all.g_ci,'NewVariableNames','phase_CI');
% 
% tt = removevars(tt,'Var1');
% 
% writetable(tt,"SF_Jun2014_to_Jun2018_Utide.csv")
% 
% tt = table(zeros(length(sf_coef_all.name),1));
% tt = addvars(tt,fr_coef_all.name,'NewVariableNames','name');
% tt = addvars(tt,fr_coef_all.A,'NewVariableNames','amp');
% tt = addvars(tt,fr_coef_all.A_ci,'NewVariableNames','amp_CI');
% tt = addvars(tt,fr_coef_all.g,'NewVariableNames','phase');
% tt = addvars(tt,fr_coef_all.g_ci,'NewVariableNames','phase_CI');
% 
% tt = removevars(tt,'Var1');
% 
% writetable(tt,"FR_Jun2014_to_Jun2018_Utide.csv")


%% UTide moving window analysis and reconstruction
% Define moving parameters etc.
window = 30*1.0;         % window in days
incr = 20/24;%1/24;     % increment in days 
% rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.8;     % specify what percentage of data in window must be not NaN

maxFLength = 5*24;   % same cutoff period as low-pass filtering done in CWT\_Multi

stat_window=[datetime(2015,01,01,'TimeZone','UTC') end_date];

[cd,OUT] = cwt_utide(FR_dat,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

[cdSF,OUTSF] = cwt_utide(sf_data,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

% %%
% OUT.A = zeros(length(OUT.names),1);
% OUT.A_ci = zeros(length(OUT.names),1);
% OUT.g = zeros(length(OUT.names),1);
% OUT.g_ci = zeros(length(OUT.names),1);
% OUT.SNR_avg = zeros(length(OUT.names),1);
% 
% n=1;
% for k=1:length(OUT.names)
%     nombre = invalidNameIn(OUT.names{k});
%     OUT.A(n) = nanmean(cd.(nombre).amp);
%     OUT.A_ci(n) = nanmean(cd.(nombre).ampci);
%     OUT.g(n) = nanmean(cd.(nombre).phase);
%     OUT.g_ci(n) = nanmean(cd.(nombre).phaseci);
%     OUT.SNR_avg(n) = cd.(nombre).avgSNR;
%     n=n+1;
% end
% 
% OUTSF.A = zeros(length(OUTSF.names),1);
% OUTSF.A_ci = zeros(length(OUTSF.names),1);
% OUTSF.g = zeros(length(OUTSF.names),1);
% OUTSF.g_ci = zeros(length(OUTSF.names),1);
% OUTSF.SNR_avg = zeros(length(OUTSF.names),1);
% 
% n=1;
% for k=1:length(OUTSF.names)
%     nombre = invalidNameIn(OUTSF.names{k});
%     OUTSF.A(n) = nanmean(cdSF.(nombre).amp);
%     OUTSF.A_ci(n) = nanmean(cdSF.(nombre).ampci);
%     OUTSF.g(n) = nanmean(cdSF.(nombre).phase);
%     OUTSF.g_ci(n) = nanmean(cdSF.(nombre).phaseci);
%     OUTSF.SNR_avg(n) = cdSF.(nombre).avgSNR;
%     n=n+1;
% end
% 
% %% export avg. values for UTide output 
% % export to csv for David
% tt = table(zeros(length(OUTSF.names),1));
% tt = addvars(tt,OUTSF.names,'NewVariableNames','name');
% tt = addvars(tt,OUTSF.A,'NewVariableNames','amp');
% tt = addvars(tt,OUTSF.A_ci,'NewVariableNames','amp_CI');
% tt = addvars(tt,OUTSF.g,'NewVariableNames','phase');
% tt = addvars(tt,OUTSF.g_ci,'NewVariableNames','phase_CI');
% tt = addvars(tt,OUTSF.SNR_avg,'NewVariableNames','avg_SNR');
% 
% tt = removevars(tt,'Var1');
% 
% writetable(tt,"SF_Jun2014_to_Jun2018_Utide.csv")
% 
% tt = table(zeros(length(OUT.names),1));
% tt = addvars(tt,OUT.names,'NewVariableNames','name');
% tt = addvars(tt,OUT.A,'NewVariableNames','amp');
% tt = addvars(tt,OUT.A_ci,'NewVariableNames','amp_CI');
% tt = addvars(tt,OUT.g,'NewVariableNames','phase');
% tt = addvars(tt,OUT.g_ci,'NewVariableNames','phase_CI');
% tt = addvars(tt,OUT.SNR_avg,'NewVariableNames','avg_SNR');
% 
% tt = removevars(tt,'Var1');
% 
% writetable(tt,"FR_Jun2014_to_Jun2018_Utide.csv")


%% change filter lengths (hoping to improve rmse and residual sprectra)
% altLength = 1080;
% altLengthToo = 363;
% altLengthForSF = 1080;
% altShortLength=65;
% coFiltLength = [altLength,altLength,601,altLengthToo,altLength,altLengthToo,altLength,altLength,altLength,...
%     altLength,altLength,altLengthToo,altLengthToo,altLength,altLengthToo,altLength,altLength,altLength,altLengthToo,...
%     altLength,altLength,altLengthToo,altLength,altLength,altShortLength,altShortLength,altShortLength,altShortLength,...
%     altShortLength,altShortLength,altShortLength,altShortLength];

% constitsD1Names = ["Alp1","TwoQ1","Q1","O1","NO1","K1","J1","OO1","Ups1"];
% constitsD2Names = ["Eps2","Mu2","N2","M2","L2","S2","Eta2"];
% constitsD3Names = ["MO3","M3","MK3","SK3"];
% constitsD4PlusNames = ["MN4","M4","MS4","S4","MK5","M6","MK7","M8","MK9","M10","MK11","M12"];
% constitsNames = horzcat(constitsD1Names,constitsD2Names,constitsD3Names,constitsD4PlusNames);

%% run cwt routine on Vancouver/Astoria data

disp("Beginning CWT\_Multi analysis")

% run CWT analysis/dynamic inference on False River data
[constits,species,ref,admit] = cwtMulti(FR_dat.dates,FR_dat.wl,'dynamicInference');%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

% Run admittance on False River data, using SF Bay data as reference
[adConstits,adSpecies,adRef,adAdmit] = cwtMulti(FR_dat.dates,FR_dat.wl,'performAdmittance',sf_data.lon,sf_data.lat,'refStation',sf_data.wl);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);
[~,~,tpRef,tpAdmit] = cwtMulti(FR_dat.dates,FR_dat.wl,'performAdmittance',sf_data.lon,sf_data.lat);

% run full CWT analysis on SF Bay data
[SFconstits,SFspecies,~,~] = cwtMulti(sf_data.dates,sf_data.wl);%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

constits.decTimesAll.TimeZone='UTC';
adSpecies.decTimesAll.TimeZone='UTC';
species.decTimesAll.TimeZone='UTC';
constits.alltimes.TimeZone='UTC';
cd.datetimes.TimeZone='UTC';
SFconstits.decTimesAll.TimeZone='UTC';
%% test on new SF Bay data
% load('./sf_data.mat')
% [SFconstits,SFspecies,~,~] = cwtMulti(sf_data.dates,sf_data.wl,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

[~,start_date_ind] = min(abs(constits.decTimesAll - datetime(2015,01,01,'TimeZone','UTC')));

%% reconstruction comparison
% p=figure();
% plot(constits.alltimes,constits.reconstruction.reconHi,'DisplayName','CWT\_Multi reconstruction')
% % yyaxis right
% hold on
% plot(constits.alltimes,OUT.wl,'DisplayName','UTide reconstruction')
% legend('location','northeast')
% grid on
% xlim([datetime(2015,01,01) end_date])
% ax = gca;
% ax.FontSize = 14;
% set(gcf,'Position',[100 100 1500 600])
% saveas(p,'./figs/reconCompSF.png')
% 
% p=figure();
% plot(constits.alltimes,OUT.wl-constits.reconstruction.reconHi)
% grid on
% xlim([datetime(2015,01,01) end_date])
% title('UTide reconstruction - CWT\_Multi reconstruction')
% ax = gca;
% ax.FontSize = 14;
% set(gcf,'Position',[100 100 1500 600])
% saveas(p,'./figs/reconResidCompSF.png')
% 
% p = figure();
% a = constits.dataIn_hi(1:end-1)-OUT.wl(2:end).';
% plot(constits.alltimes(1:end-1),a)
% hold on
% plot(constits.alltimes,constits.dataIn_hi-constits.reconstruction.reconHi.')
% hold off


%% defining plotting colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];

%% Water level plots
figure(); p=tiledlayout(2,1);

ax1=nexttile;
plot(sf_data.dates,sf_data.wl,'color',blue,'linewidth',2,'DisplayName','water level')
hold(ax1,'on')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
grid on
legend('location','northeast')
title('SF Bay water level', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(FR_dat.dates,FR_dat.wl,'color',blue,'linewidth',2,'DisplayName','water level')
hold(ax2,'on')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
grid on
legend('location','northeast')
title('False River water level', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
% set(gca,'xtick',[])
axis tight
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/SF_FR_water_levels.png')

%% plot admittance
if isempty(adAdmit)==0
    p=figure();
    x = constits.decTimesAll;
    x.TimeZone='UTC';
    plot(x,adAdmit.constits.M2.amps,'linewidth',2,'DisplayName','admittance: FR/SF water level')
    hold on
    plot(x,tpAdmit.constits.M2.amps,'linewidth',2,'DisplayName','admittance: FR/SF tidal potential')
    xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
    xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
    title('M_2 admittance', 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
    hold off
    grid on
    ylim([0.15 0.7])
    ylabel('Admittance amplitude')
    xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
    ax = gca;
    ax.FontSize = 14;
    legend('location','northeast')
    set(gcf,'Position',[100 100 1500 600])
    saveas(p,'./figs/SF_FR_M2_admittance.png')
end

%% wavelet "energy" plot
p=figure();
yyaxis right
plot(abs(SFconstits.M2.filter),'k--','linewidth',2,'DisplayName','M_2 wavelet (abs; R)')
ylabel('Magnitude')
ylim([5.3e-7 1.5e-4])
yyaxis left
plot(real(SFconstits.M2.filter),'xb-','linewidth',1,'DisplayName','M_2 wavelet (real; L)')
hold on
plot(imag(SFconstits.M2.filter),'sr-','linewidth',1,'DisplayName','M_2 wavelet (imag; L)')
hold off
n=floor(length(real(SFconstits.M2.filter))/4);
lower_bound=n;upper_bound=length(real(SFconstits.M2.filter))-n;
xline(lower_bound,'linewidth',2,'color',[0 0 0]+0.5,'HandleVisibility','off')
xline(upper_bound,'linewidth',2,'color',[0 0 0]+0.5,'DisplayName','80% energy window')
ylim([-1.75e-4 2e-4])
xlim([1 length(SFconstits.M2.filter)])
ylabel('Amplitude')
xlabel('Index no.')
legend('location','northeast')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.FontSize = 14;
set(gcf,'Position',[100 100 1500 600])
% saveas(p,'./figs/M2_wavelet.png')

%% two dynamically inferred ONLY constituents (P1 and K2) for False River
% DI AST
figure(); p=tiledlayout(2,1);

% P1 dyn inf

ax1=nexttile;
plot(constits.decTimesAll,constits.dynamicInference.P1.sixMoAmp(2:end-1),'linewidth',2,'DisplayName','P_1 CWT\_Multi (6 month filter)')
hold(ax1,'on')
plot(constits.decTimesAll,constits.dynamicInference.P1.shortAmps,'linewidth',2,'DisplayName','P_1 CWT\_Multi (dynamically inferred 2 weeks filter)')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
title('P_1 dynamic inference', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
legend('location','northeast')
grid on
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(constits.decTimesAll,constits.dynamicInference.K2.sixMoAmp(2:end-1),'linewidth',2,'DisplayName','K_2 CWT\_Multi (6 month filter)')
hold(ax2,'on')
plot(constits.decTimesAll,constits.dynamicInference.K2.shortAmps,'linewidth',2,'DisplayName','K_2 CWT\_Multi (dynamically inferred 2 weeks filter)')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
grid on
legend('location','northeast')
title('K_2 dynamic inference', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

set(gcf,'Position',[100 100 2500 1750])
saveas(p,'./figs/FR_dynamic_inference_D1_amps.png')

%% species and river flow results
figure(); p=tiledlayout(3,1);
dayFlow.dtime.TimeZone='UTC';
ax1=nexttile;
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2,'color','black')
hold(ax1,'on')
qw{1}=plot(nan,'color','black','linewidth',2);
qw{2}=plot(nan,'color',[0 0.4470 0.7410],'linewidth',2); qw{3}=plot(nan,'color',[0.85 0.325 0.098],'linewidth',2);
qw{4}=plot(nan,'color',[0.929 0.694 0.125],'linewidth',2);qw{5}=plot(nan,'color',[0.494,0.184,0.556],'linewidth',2);
qw{6}=plot(nan,'r--','linewidth',2); qw{7}=plot(nan,'g--','linewidth',2);
hold(ax1,'off')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
title('River flow', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylabel('NDOI discharge (m^3/s)')
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

legend([qw{:}], {'River flow','D_1','D_2','D_3','D_4','construction start','flood start'},'Location','northeast')


% xlim([constits.decTimesAll(1) constits.decTimesAll(end)])

ax2=nexttile;
plot(adSpecies.decTimesAll,adSpecies.D1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax2,'on')
plot(adSpecies.decTimesAll,adSpecies.D2.amps,'linewidth',2,'DisplayName','D_2')
plot(adSpecies.decTimesAll,adSpecies.D3.amps,'linewidth',2,'DisplayName','D_3')
plot(adSpecies.decTimesAll,adSpecies.D4.amps,'linewidth',2,'DisplayName','D_4')
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2,'color','black','DisplayName','River flow')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
ylabel('Amplitude (m)')
grid on
% legend('Location','northoutside') %'Position',[0.8,1.8,0.4,0.3])  %'location','northeast')
title('False River species results', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ylim([0 0.45])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(adSpecies.decTimesAll,adRef.species.D1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax3,'on')
plot(species.decTimesAll,adRef.species.D2.amps,'linewidth',2,'DisplayName','D_2')
plot(species.decTimesAll,adRef.species.D3.amps,'linewidth',2,'DisplayName','D_3')
plot(species.decTimesAll,adRef.species.D4.amps,'linewidth',2,'DisplayName','D_4')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
ylabel('Amplitude (m)')
grid on
% legend('location','northeast')
title('SF Bay species results', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
set(gcf,'Position',[100 100 2400 2000])
plotAnn('(c)')

saveas(p,'./figs/SF_FR_flow_and_species.png')

%% residual spectrum
figure(); p=tiledlayout(2,1);
sp_names = ["D_1","D_2","D_3","D_4","D_5","D_6","D_7","D_8","D_9","D_{10}","D_{11}","D_{12}"];

ax1=nexttile;
loglog(constits.input_freqs.^-1,constits.input_pwr,'r--','linewidth',1,'DisplayName','input data')
hold(ax1,'on')
plot(OUT.resid_freqs.^-1,OUT.resid_pwr,'b-','linewidth',2,'DisplayName','UTide residual')
plot(constits.fband.freqs.^-1,constits.fband.pwrs,'k-','linewidth',2,'DisplayName','CWT\_Multi residual')
for i=1:length(species.input.omegas(3:end))
    x = (species.input.omegas(i+2)/2/pi)^-1/3600/24;
    xline(x,'k--','linewidth',2,'alpha',0.5,'HandleVisibility','off')
    text(x, min(ylim),sp_names(i), 'Horiz','left', 'Vert','bottom','FontSize',16) 
end
hold(ax1,'off')
grid on
xlim([0 5])
ylim auto
yma = max(constits.fband.pwrs)*1.2;
ymi = min(constits.fband.pwrs);
set(gca,'xdir','reverse')

legend('Location','northeast')
% axis tight
title('False River spectra', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ax = gca;
ax.FontSize = 14;
plotAnnSpec('(a)')

ax2=nexttile;
loglog(constits.input_freqs.^-1,SFconstits.input_pwr,'r--','linewidth',1,'DisplayName','input data')
hold(ax2,'on')
plot(OUT.resid_freqs.^-1,OUTSF.resid_pwr,'b-','linewidth',2,'DisplayName','UTide residual')
plot(constits.fband.freqs.^-1,SFconstits.fband.pwrs,'k-','linewidth',2,'DisplayName','CWT\_Multi residual')
for i=1:length(species.input.omegas(3:end))
    x = (species.input.omegas(i+2)/2/pi)^-1/3600/24;
    xline(x,'k--','linewidth',2,'alpha',0.5,'HandleVisibility','off')
    text(x, min(ylim),sp_names(i), 'Horiz','left', 'Vert','bottom','FontSize',16) 
end
hold(ax2,'off')
grid on
xlim([0 5])
ylim auto
yma = max(constits.fband.pwrs)*1.2;
ymi = min(constits.fband.pwrs);
set(gca,'xdir','reverse')
plotAnnSpec('(b)')
xlabel('Period [days]')
legend('Location','northeast')
% axis tight
title('SF Bay spectra', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
ax = gca;
ax.FontSize = 14;

ylabel(p,'PSD [m^2 cpd^{-1}]','FontSize',18)
set(gcf,'Position',[100 100 1800 1000])
saveas(p,'./figs/SF_FR_spectra.png')

%% compare reconstruction & recon statistics
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(constits.alltimes,constits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax1,'on')
plot(constits.alltimes,OUT.wlAll,'DisplayName','UTide reconstruction')
hold(ax1,'off')
grid on
legend('location','northeast')
ylabel('Water level [m]')
title('False River reconstruction', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(constits.alltimes,FR_dat.wl.' - constits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.0 1.0])
yyaxis right
plot(OUT.dtimesTrim,OUT.resid,'DisplayName','UTide error')
ylim([-0.5 1.5])
grid on
legend('location','northeast')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
title('False River reconstruction error', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(constits.alltimes,SFconstits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax3,'on')
plot(OUTSF.dtimes,OUTSF.wlAll,'DisplayName','UTide reconstruction')
hold(ax3,'off')
grid on
legend('location','northeast')
title('SF Bay reconstruction', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ylabel('Water level [m]')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(constits.alltimes,sf_data.wl.' - SFconstits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.0 1.0])
yyaxis right
plot(OUTSF.dtimesTrim,OUTSF.resid,'DisplayName','UTide error')
ylim([-0.5 1.5])
grid on
legend('location','northeast')
title('SF Bay reconstruction error', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
ax = gca;
ax.FontSize = 14;
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plotAnn('(d)')

set(gcf,'Position',[100 100 2500 2000])
saveas(p,'./figs/SF_FR_UTide_CWT_reconstruction.png')

%% FALSE RIVER DIURNAL AMPLITUDE
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2)
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
title('River flow', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylabel('NDOI discharge (m^3/s)')
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax1=nexttile;
plot(cd.datetimes,cd.Q1.amp,'color',blue,'linewidth',2,'DisplayName','Q_1 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
plotCI(constits.decTimesAll(start_date_ind:end),constits.Q1.amps(start_date_ind:end),constits.fband.ci(3),2)
hold(ax1,'off')
legend('location','northeast')
title('False River Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')


ax2=nexttile;
plot(cd.datetimes,cd.O1.amp,'color',blue,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.O1.amps(start_date_ind:end),constits.fband.ci(4),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River O_1', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylabel('Amplitude (m)','FontSize',18)
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

% K1
ax3=nexttile;
plot(cd.datetimes,cd.K1.amp,'color',blue,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.K1.amps(start_date_ind:end),constits.fband.ci(6),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('False River K_1', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
legend('location','northeast')
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/FR_UTide_CWT_D1_amps.png')

%% FALSE RIVER SEMIDIURNAL AMPLITUDE

figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(cd.datetimes,cd.M2.amp,'color',blue,'linewidth',2,'DisplayName','M_2 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.M2.amps(start_date_ind:end),constits.fband.ci(13),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax1,'off')
ylim([0.23 0.33])
legend('location','northeast')
title('False River M_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cd.N2.amp,'color',blue,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.N2.amps(start_date_ind:end),constits.fband.ci(12),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River N_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,cd.S2.amp,'color',blue,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.S2.amps(start_date_ind:end),constits.fband.ci(15),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('False River S_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,cd.L2.amp,'color',blue,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.L2.amps,'color',red,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.L2.amps(start_date_ind:end),constits.fband.ci(14),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('False River L_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/FR_UTide_CWT_D2_amps.png')


%% FALSE RIVER OVERTIDE AMPLITUDES
figure();p=tiledlayout(2,1);

ax1=nexttile;
plot(cd.datetimes,cd.MK3.amp,'color',blue,'linewidth',2,'DisplayName','MK_3 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.MK3.amps(start_date_ind:end),constits.fband.ci(19),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('False River MK_3', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

% plot results (M4)
ax2=nexttile;
plot(cd.datetimes,cd.M4.amp,'color',blue,'linewidth',2,'DisplayName','M_4 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(constits.decTimesAll,constits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.M4.amps(start_date_ind:end),constits.fband.ci(22),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River M_4', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 0 1500 1000])
saveas(p,'./figs/FR_UTide_CWT_overtide_amps.png')

%% SF BAY DIURNAL AMPLITUDES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(cd.datetimes,cdSF.Q1.amp,'color',blue,'linewidth',2,'DisplayName','Q_1 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.Q1.amps(start_date_ind:end),SFconstits.fband.ci(3),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('SF Bay Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cdSF.O1.amp,'color',blue,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.O1.amps(start_date_ind:end),SFconstits.fband.ci(4),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1
ax3=nexttile;
plot(cd.datetimes,cdSF.K1.amp,'color',blue,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.K1.amps(start_date_ind:end),SFconstits.fband.ci(6),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('SF Bay K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
legend('location','northeast')
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/SF_UTide_CWT_D1_amps.png')


%% SF BAY SEMIDIURNAL AMPLITUDES
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(cd.datetimes,cdSF.M2.amp,'color',blue,'linewidth',2,'DisplayName','M_2 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M2.amps(start_date_ind:end),SFconstits.fband.ci(13),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax1,'off')
ylim([0.49 0.63])
legend('location','northeast')
title('SF Bay M_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cdSF.N2.amp,'color',blue,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.N2.amps(start_date_ind:end),SFconstits.fband.ci(12),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
ylim([0.085 0.17])
legend('location','northeast')
title('SF Bay N_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,cdSF.S2.amp,'color',blue,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
axis tight
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.S2.amps(start_date_ind:end),SFconstits.fband.ci(15),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('SF Bay S_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,cdSF.L2.amp,'color',blue,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.L2.amps,'color',red,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.L2.amps(start_date_ind:end),SFconstits.fband.ci(14),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('SF Bay L_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/SF_UTide_CWT_D2_amps.png')


%% SF BAY OVERTIDE AMPLITUDES
figure();p=tiledlayout(2,1);

ax1=nexttile;
plot(cd.datetimes,cdSF.MK3.amp,'color',blue,'linewidth',2,'DisplayName','MK_3 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.MK3.amps(start_date_ind:end),SFconstits.fband.ci(19),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('SF Bay MK_3', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

% plot results (M4)
ax2=nexttile;
plot(cd.datetimes,cdSF.M4.amp,'color',blue,'linewidth',2,'DisplayName','M_4 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
plot(SFconstits.decTimesAll,SFconstits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M4.amps(start_date_ind:end),SFconstits.fband.ci(22),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay M_4', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 0 1500 1000])
saveas(p,'./figs/SF_UTide_CWT_overtide_amps.png')

%% SF BAY DIURNAL PHASES 

figure(); p=tiledlayout(3,1);

ax1=nexttile;
plotPhase(SFconstits.decTimesAll,SFconstits.Q1.phases,cd.datetimes,OUTSF.Q1.phases_alt,ax1,SFconstits.fband.ci_ph(3),2,'Q_1',-145)
hold(ax1,'on')
% plot(cd.datetimes,OUTSF.Q1.phases_alt,'linewidth',2,'DisplayName','Q_1 UTide')
% plot(SFconstits.decTimesAll,SFconstits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.Q1.phases(start_date_ind:end),SFconstits.fband.ci_ph(3),2)
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('SF Bay Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(cd.datetimes,OUTSF.O1.phases_alt,'linewidth',2,'DisplayName','O_1 UTide')
% plot(SFconstits.decTimesAll,SFconstits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.O1.phases(start_date_ind:end),SFconstits.fband.ci_ph(4),2)
plotPhase(SFconstits.decTimesAll,SFconstits.O1.phases,cd.datetimes,OUTSF.O1.phases_alt,ax2,SFconstits.fband.ci_ph(4),2,'O_1',-145)
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.amps,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1
ax3=nexttile;
% plot(cd.datetimes,OUTSF.K1.phases_alt,'linewidth',2,'DisplayName','K_1 UTide')
plotPhase(SFconstits.decTimesAll,SFconstits.K1.phases,cd.datetimes,OUTSF.K1.phases_alt,ax3,SFconstits.fband.ci_ph(6),2,'K_1',-145)
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(SFconstits.decTimesAll,SFconstits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.K1.phases(start_date_ind:end),SFconstits.fband.ci_ph(6),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('SF Bay K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
legend('location','northeast')
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/SF_UTide_CWT_D1_phases.png')

%% SF BAY SEMIDIURNAL PHASES
figure(); p=tiledlayout(4,1);

ax1=nexttile;
% plot(cd.datetimes,OUTSF.M2.phases_alt,'linewidth',2,'DisplayName','M_2 UTide')
plotPhase(SFconstits.decTimesAll,SFconstits.M2.phases,cd.datetimes,OUTSF.M2.phases_alt,ax1,SFconstits.fband.ci_ph(13),2,'M_2',-145)

hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(SFconstits.decTimesAll,SFconstits.M2.phases,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M2.phases(start_date_ind:end),SFconstits.fband.ci_ph(13),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.phases,'r--','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('SF Bay M_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(cd.datetimes,OUTSF.N2.phases_alt,'linewidth',2,'DisplayName','N_2 UTide')
plotPhase(SFconstits.decTimesAll,SFconstits.N2.phases,cd.datetimes,OUTSF.N2.phases_alt,ax2,SFconstits.fband.ci_ph(12),2,'N_2',-145)
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(SFconstits.decTimesAll,SFconstits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.N2.phases(start_date_ind:end),SFconstits.fband.ci_ph(12),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.phases,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay N_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
% plot(cd.datetimes,OUTSF.S2.phases_alt,'linewidth',2,'DisplayName','S_2 UTide')
plotPhase(SFconstits.decTimesAll,SFconstits.S2.phases,cd.datetimes,OUTSF.S2.phases_alt,ax3,SFconstits.fband.ci_ph(15),2,'S_2',-145)
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(SFconstits.decTimesAll,SFconstits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.S2.phases(start_date_ind:end),SFconstits.fband.ci_ph(15),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('SF Bay S_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
% plot(cd.datetimes,OUTSF.L2.phases_alt,'linewidth',2,'DisplayName','L_2 UTide')
plotPhase(SFconstits.decTimesAll,SFconstits.L2.phases,cd.datetimes,OUTSF.L2.phases_alt,ax4,SFconstits.fband.ci_ph(14),2,'L_2',-145)
hold(ax4,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(SFconstits.decTimesAll,SFconstits.L2.phases,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
% plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.L2.phases(start_date_ind:end),SFconstits.fband.ci_ph(14),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('SF Bay L_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/SF_UTide_CWT_D2_phases.png')



%% FALSE RIVER DIURNAL PHASES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
% plot(cd.datetimes,OUT.Q1.phases_alt,'linewidth',2,'DisplayName','Q_1 UTide')
plotPhase(constits.decTimesAll,constits.Q1.phases-40,cd.datetimes,OUT.Q1.phases_alt-40,ax1,constits.fband.ci_ph(3),2,'Q_1',-145)
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.Q1.phases(start_date_ind:end),constits.fband.ci_ph(3),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('False River Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(cd.datetimes,OUT.O1.phases_alt,'linewidth',2,'DisplayName','O_1 UTide')
plotPhase(constits.decTimesAll,constits.O1.phases,cd.datetimes,OUT.O1.phases_alt,ax2,constits.fband.ci_ph(4),2,'O_1',-145)
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.O1.phases(start_date_ind:end),constits.fband.ci_ph(4),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.phases,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1 
ax3=nexttile;
% plot(cd.datetimes,OUT.K1.phases_alt,'linewidth',2,'DisplayName','K_1 UTide')
plotPhase(constits.decTimesAll,constits.K1.phases,cd.datetimes,OUT.K1.phases_alt,ax3,constits.fband.ci_ph(6),2,'K_1',-145)
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.K1.phases(start_date_ind:end),constits.fband.ci_ph(6),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('False River K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
legend('location','northeast')
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/FR_UTide_CWT_D1_phases.png')


%% FALSE RIVER SEMIDIURNAL PHASES

figure(); p=tiledlayout(4,1);

ax1=nexttile;
% plot(cd.datetimes,OUT.M2.phases_alt,'linewidth',2,'DisplayName','M_2 UTide')
plotPhase(constits.decTimesAll,constits.M2.phases,cd.datetimes,OUT.M2.phases_alt,ax1,constits.fband.ci_ph(13),2,'M_2',-145)
hold(ax1,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.M2.phases,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.M2.phases(start_date_ind:end),constits.fband.ci_ph(13),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.phases,'r--','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('False River M_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(cd.datetimes,OUT.N2.phases_alt,'linewidth',2,'DisplayName','N_2 UTide')
plotPhase(constits.decTimesAll,constits.N2.phases,cd.datetimes,OUT.N2.phases_alt,ax2,constits.fband.ci_ph(12),2,'N_2',-145)
hold(ax2,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.N2.phases(start_date_ind:end),constits.fband.ci_ph(12),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.D2.phases,'r--','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River N_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
% plot(cd.datetimes,OUT.S2.phases_alt,'linewidth',2,'DisplayName','S_2 UTide')
plotPhase(constits.decTimesAll,constits.S2.phases,cd.datetimes,OUT.S2.phases_alt,ax3,constits.fband.ci_ph(15),2,'S_2',-145)
hold(ax3,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.S2.phases(start_date_ind:end),constits.fband.ci_ph(15),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('False River S_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
set(gca,'xtick',[])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
% plot(cd.datetimes,OUT.L2.phases_alt,'linewidth',2,'DisplayName','L_2 UTide')
plotPhase(constits.decTimesAll,constits.L2.phases,cd.datetimes,OUT.L2.phases_alt,ax4,constits.fband.ci_ph(14),2,'L_2',-145)
hold(ax4,'on')
xlim([datetime(2015,01,01,'TimeZone','UTC') end_date])
% plot(constits.decTimesAll,constits.L2.phases,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
% plotCI(constits.decTimesAll(start_date_ind:end),constits.L2.phases(start_date_ind:end),constits.fband.ci_ph(14),2)
xline(datetime(2015,05,08,'TimeZone','UTC'),'r--','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07,'TimeZone','UTC'),'g--','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('False River L_2', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
grid on
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/FR_UTide_CWT_D2_phases.png')

%% BARE FILTER PLOT
p=figure();
yyaxis right
% plot(abs(SFconstits.M2.filter),'k--','linewidth',2,'DisplayName','M_2 wavelet (abs; R)')
ylabel('Magnitude')
ylim([5.3e-7 1.5e-4])
yyaxis left
plot(real(SFconstits.K1.filter),'xb-','linewidth',1,'DisplayName','M_2 wavelet (real; L)')
hold on
plot(imag(SFconstits.K1.filter),'sr-','linewidth',1,'DisplayName','M_2 wavelet (imag; L)')
hold off
n=floor(length(real(SFconstits.K1.filter))/4);
lower_bound=n;upper_bound=length(real(SFconstits.M2.filter))-n;
% xline(lower_bound,'linewidth',2,'color',[0 0 0]+0.5,'HandleVisibility','off')
% xline(upper_bound,'linewidth',2,'color',[0 0 0]+0.5,'DisplayName','80% energy window')
ylim([-1.75e-4 2e-4])
xlim([1 length(SFconstits.M2.filter)])
ylabel('Amplitude')
xlabel('Index no.')
% legend('location','northeast')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.FontSize = 14;
set(gcf,'Position',[100 100 1500 600])
% saveas(p,'./figs/M2_wavelet.png')

%%
figure(); p=tiledlayout(1,1);
ax1 = nexttile;
% plotPhase(constits.decTimesAll,constits.O1.phases,cd.datetimes,OUT.O1.phases_alt,ax1,constits.fband.ci_ph(4),2,'O_1',-145)
plotPhase(constits.decTimesAll,constits.Q1.phases-40,cd.datetimes,OUT.Q1.phases_alt-40,ax1,constits.fband.ci_ph(3),2,'Q_1',-145)

% plotPhase(SFconstits.decTimesAll,SFconstits.L2.phases,cd.datetimes,OUTSF.L2.phases_alt,ax1,SFconstits.fband.ci_ph(14),2,'L_2',-145)


% plot(constits.decTimesAll,matts_unwrap(constits.Q1.phases,330))
% hold on
% plot(cd.datetimes,matts_unwrap(OUT.Q1.phases_alt,330))
% hold off
% ylim([-1000 1000])


%% %%%%%%%%%%%%%% PLOTTING FUNCTIONS
function plotAnn(txt)
    text(min(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end

function plotAnnSpec(txt)
    text(max(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end

function plotCI(x,y,cf,color_number)
%     alpha = 1.0; % was 0.3
%     if color_number==1
%         cfn = [0 0.4470 0.7410 alpha];
%     elseif color_number==2
%         cfn = [0.8500 0.3250 0.0980 alpha];
%     elseif color_number==3
%         cfn = [0.9290 0.6940 0.1250 alpha];
%     else
%         error("must pass a valid color number (1,2,3).")
%     end
%     
%     nMid = 45;
%     nHorz = 5;
% 
%     x1      = x(nMid-nHorz:nMid+nHorz);
%     xMid    = linspace(x(nMid),x(nMid)+seconds(10),2*nHorz+1);
%     yAbove  = (y(nMid)+cf)*ones(2*nHorz+1,1);
%     yBelow  = (y(nMid)-cf)*ones(2*nHorz+1,1);
%     yMid    = linspace(y(nMid)-cf,y(nMid)+cf,2*nHorz+1);
%     
%     plot(x1,yAbove,'-','color',cfn,'linewidth',2,'HandleVisibility','off')
%     plot(x1,yBelow,'-','color',cfn,'linewidth',2,'HandleVisibility','off')
%     plot(xMid,yMid,'-','color',cfn,'linewidth',2,'HandleVisibility','off')
    
    % Alternately, just plot the text values
    txt = strcat('\pm',string(round(cf,3,'significant')));
    text(max(xlim)-0.025*(max(xlim)-min(xlim)),min(ylim)+0.05*(max(ylim)-min(ylim)),txt, 'Horiz','right', 'Vert','bottom','FontSize',20,'BackgroundColor','white') 

end

function [resp,om] = cwt_freq_response(filt,cent_freq,delta_om,N_om,deltaT)
    % find frequency reponse of a filter using sine-wave analysis
    %

    N = length(filt);
    N2 = (length(filt)-1)/2;

    t = (-N2:N2) * deltaT;

    N_om2 = (N_om - 1)/2;
    om = (cent_freq - delta_om*N_om2:delta_om:cent_freq+delta_om*N_om2);

    resp = zeros(length(om),1);

    for k=1:N_om
        c = cos(om(k) * t);

        % for NConvolved() data must be same dim
        resp(k) = real(NConvolved([],c,filt,ones(size(filt)),1,0.5));
    end

end


%% FUNCTIONS
function nameOut = invalidNameIn(nameIn)
    % only goes up to name beginning with '9'
    if nameIn(1)=='1'
        nameOut = strcat('One',nameIn(2:end));
    elseif nameIn(1)=='2'
        nameOut = strcat('Two',nameIn(2:end));
    elseif nameIn(1)=='3'
        nameOut = strcat('Three',nameIn(2:end));
    elseif nameIn(1)=='4'
        nameOut = strcat('Four',nameIn(2:end));
    elseif nameIn(1)=='5'
        nameOut = strcat('Five',nameIn(2:end));
    elseif nameIn(1)=='6'
        nameOut = strcat('Six',nameIn(2:end));
    elseif nameIn(1)=='7'
        nameOut = strcat('Seven',nameIn(2:end));
    elseif nameIn(1)=='8'
        nameOut = strcat('Eight',nameIn(2:end));
    elseif nameIn(1)=='9'
        nameOut = strcat('Nine',nameIn(2:end));
    else
        nameOut = nameIn;
    end
end


function plotAnnSE(txt,btm)
    % add text to upper left corner of any plot
    % xspot = max(xlim)-0.1*(max(xlim)-min(xlim));
    text(max(xlim)-0.025*(max(xlim)-min(xlim)),btm,txt, 'Horiz','right', 'Vert','bottom','FontSize',20,'BackgroundColor','white') 
end

function plotPhase(x1,y1,x2,y2,ax_in,ci,co_no,const_name,ci_txt_bottom)
    % plotting phases after adjusting mean of y2 to match mean of y1
    jump = 320;
    % a = unwrap(y2,jump) - (nanmedian(unwrap(y2,jump)) - nanmedian(unwrap(y1,jump)));
%     a = unwrap(y2,jump) - (nanmedian(y2) - nanmedian(y1));
%     y2_new = mod(a,360);
   
    y2 = y2 - (nanmedian(y2) - nanmedian(y1));
    
    
%     for i=1:(length(y2_new)-1)
%         diff = 
    if const_name=="Q1"
        y1 = y1+40;
        y2=y2+40;
    end
    plot(x2,matts_unwrap(y2,jump),'linewidth',2,'DisplayName',strcat(const_name,' UTide'))
    hold(ax_in,'on')
    plot(x1,matts_unwrap(y1,jump),'linewidth',2,'DisplayName',strcat(const_name,' CWT\_Multi'))
    % plotCI(x1,y1,ci,co_no)
    plotAnnSE(strcat('\pm',string(round(ci,3,'significant'))),ci_txt_bottom)
    hold(ax_in,'off')
    
end


function sig = matts_unwrap(sig,jump)
    nw = 1500;
    
    for i=flip([1:25:nw])
        sig = find_jumps(i,sig,jump);
    end
    

    sig = selective_mod(100,sig);

end

function sig = find_jumps(nw,sig,jump)
    for i=(nw+2):(length(sig)-nw-2)
        sig_here = sig(i-nw:i+nw);
        if max(diff(sig_here))>jump && min(diff(sig_here))< (-jump)
            [~,start_jump] = max(diff(sig_here));
            [~,end_jump] = min(diff(sig_here));
            
            if start_jump>end_jump
                if max(abs(mean(sig_here(end_jump:start_jump)))-abs(sig_here(end_jump)+jump))>0 && abs(sig_here(start_jump+1)-sig_here(end_jump))<jump
                    start_jump=start_jump; end_jump=end_jump+1;
                    sig(i-nw+end_jump-1:i-nw+start_jump-1) = sig_here(end_jump:start_jump) - (sig_here(end_jump)-sig_here(end_jump-1));
                end
            else
                if max(abs(mean(sig_here(start_jump:end_jump)))-abs(sig_here(start_jump)-jump))<1000 && (abs(sig_here(end_jump+1)-sig_here(start_jump))<jump)
                    start_jump=start_jump+1; end_jump=end_jump;
                    sig(i-nw+start_jump-1:i-nw+end_jump-1) = sig_here(start_jump:end_jump) - (sig_here(start_jump)-sig_here(start_jump-1)) ;
                end
            end
        end
    end
    

end

function sig = selective_mod(nw,sig)
    for i=(nw+1):(length(sig)-nw-1)
        sig_here = sig(i-nw:i+nw);
        if (sig_here(1)-sig_here(end))>350 && (sig_here(1)>0)
            [~,ii] = min(abs(sig_here));
            sig(i-nw+ii-1:end) = sig(i-nw+ii-1:end)+360;
        end
    end
    
    for i=(nw+1):(length(sig)-nw-1)
        sig_here = sig(i-nw:i+nw);
        if (sig_here(1)-sig_here(end))< -350 && (sig_here(1)>0)
            [~,ii] = min(abs(sig_here-360));
            sig(i-nw+ii-1:end) = sig(i-nw+ii-1:end)-360;
        end
    end
end




