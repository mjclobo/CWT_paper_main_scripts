% Comparing CWT_Multi to UTide
% clear
% [constits,species,ref,admit] = cwtMulti(dat4.dates,dat4.wl);

%% load data
addpath('../cwtMultiPackage/cwtMulti')

load('./data/Lobo_Data.mat')
load('./data/dayflow_1997to2020_NDOI_RIO_WEST_XGEO.mat')
dayFlow.dtime = datetime(dayFlow.dn,'ConvertFrom','datenum','Format','HH:mm:ss.SSS');

dat4.dates  = datetime(d.tn{4},'ConvertFrom','datenum');
dat4.wl     = d.wl{4};

dat1.dates  = datetime(d.tn{1},'ConvertFrom','datenum');
dat1.wl     = d.wl{1};

%% align datasets to common start/end
startDate = datetime(2014,06,01); % change this back to 2013,06 to avoid edge effects
endDate = datetime(2018,06,01);  % 2018,06,01

dat4.start = find(dat4.dates==startDate);
dat4.end = find(dat4.dates==endDate);

dat4.dates = dat4.dates(dat4.start:dat4.end);
dat4.wl = dat4.wl(dat4.start:dat4.end);

dat1.start = find(dat1.dates==startDate);
dat1.end = find(dat1.dates==endDate);

dat1.dates = dat1.dates(dat1.start:dat1.end);
dat1.wl = dat1.wl(dat1.start:dat1.end);

nan_ind = unique(sort([find(isnan(dat4.wl));find(isnan(dat1.wl))]));

dat4.wl(nan_ind) = NaN;
dat1.wl(nan_ind) = NaN;

% lon, lat for dat1oria
dat4.lon = d.info{4}.lon;            % W
dat4.lat = d.info{4}.lat;             % N

dat1.lon = d.info{1}.lon;            % W
dat1.lat = d.info{1}.lat;             % N

% set dates for plotting
end_date = datetime(2017,10,01);

%% UTide analysis and reconstruction
addpath('./utide_tool/')

dat1.datenums = datenum(dat1.dates);
dat4.datenums = datenum(dat4.dates);
dat1.dtime = dat1.dates;
dat4.dtime = dat4.dates;

dat1.dates.TimeZone='UTC';
dat4.dates.TimeZone='UTC';

rayleigh = 0.9;

% %% a single call to UTide for whole time period
% sf_coef_all = ut_solv(dat1.datenums, dat1.wl,[],dat1.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
% 
% fr_coef_all = ut_solv(dat4.datenums, dat4.wl,[],dat4.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
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
% writetable(tt,"SF_Jan2016_to_Jan2017_Utide.csv")
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
% writetable(tt,"FR_Jan2016_to_Jan2017_Utide.csv")


%% UTide moving window analysis and reconstruction
% Define moving parameters etc.
window = 30*1.0;         % window in days
incr = 20/24;%1/24;     % increment in days 
% rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.8;     % specify what percentage of data in window must be not NaN

maxFLength = 5*24;   % same cutoff period as low-pass filtering done in CWT\_Multi

stat_window=[datetime(2015,01,01) end_date]; stat_window.TimeZone='UTC';

[cd,OUT] = cwt_utide(dat4,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

[cdSF,OUTSF] = cwt_utide(sf_data,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1


%% change filter lengths (hoping to improve rmse and residual sprectra)
altLength = 1080;
altLengthToo = 337;
altLengthForSF = 1080;
altShortLength=337;
coFiltLength = [altLength,altLength,altLength,altLengthToo,altLength,altLengthToo,altLength,altLength,altLength,...
    altLength,altLength,altLength,altLengthToo,altLength,altLengthToo,altLength,altLength,altLength,altLengthToo,...
    altLength,altLength,altLengthToo,altLength,altLength,altShortLength,altShortLength,altShortLength,altShortLength,...
    altShortLength,altShortLength,altShortLength,altShortLength];

% constitsD1Names = ["Alp1","TwoQ1","Q1","O1","NO1","K1","J1","OO1","Ups1"];
% constitsD2Names = ["Eps2","Mu2","N2","M2","L2","S2","Eta2"];
% constitsD3Names = ["MO3","M3","MK3","SK3"];
% constitsD4PlusNames = ["MN4","M4","MS4","S4","MK5","M6","MK7","M8","MK9","M10","MK11","M12"];
% constitsNames = horzcat(constitsD1Names,constitsD2Names,constitsD3Names,constitsD4PlusNames);

%% run cwt routine on Vancouver/Astoria data

disp("Beginning CWT\_Multi analysis")

% run CWT analysis/dynamic inference on False River data
[constits,species,ref,admit] = cwtMulti(dat4.dates,dat4.wl,'dynamicInference','coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

% Run admittance on False River data, using SF Bay data as reference
% [adConstits,adSpecies,adRef,adAdmit] = cwtMulti(dat4.dates,dat4.wl,'performAdmittance',dat1.lon,dat1.lat,'refStation',dat1.wl);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);
% [~,~,tpRef,tpAdmit] = cwtMulti(dat4.dates,dat4.wl,'performAdmittance',dat1.lon,dat1.lat);

% run full CWT analysis on SF Bay data
[SFconstits,SFspecies,~,~] = cwtMulti(dat1.dates,dat1.wl,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

%% test on new SF Bay data
load('./sf_data.mat')
[SFconstits,SFspecies,~,~] = cwtMulti(sf_data.dates,sf_data.wl,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

[~,start_date_ind] = min(abs(constits.decTimesAll - datetime(2015,01,01)));

%% reconstruction comparison
% p=figure();
% plot(constits.alltimes,constits.reconHi,'DisplayName','CWT\_Multi reconstruction')
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
% plot(constits.alltimes,OUT.wl-constits.reconHi)
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
% plot(constits.alltimes,constits.dataIn_hi-constits.reconHi.')
% hold off


%% defining plotting colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];

%% Water level plots
figure(); p=tiledlayout(2,1);

ax1=nexttile;
plot(dat1.dates,dat1.wl,'color',blue,'linewidth',2,'DisplayName','water level')
hold(ax1,'on')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
grid on
legend('location','northeast')
title('SF Bay water level', 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(dat4.dates,dat4.wl,'color',blue,'linewidth',2,'DisplayName','water level')
hold(ax2,'on')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
grid on
legend('location','northeast')
title('False River water level', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
% set(gca,'xtick',[])
axis tight
xlim([datetime(2015,01,01) end_date])
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
    plot(x,adAdmit.constits.M2.amps,'linewidth',2,'DisplayName','admittance: FR/SF water level')
    hold on
    plot(x,tpAdmit.constits.M2.amps,'linewidth',2,'DisplayName','admittance: FR/SF tidal potential')
    xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
    xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
    title('M_2 admittance', 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
    hold off
    grid on
    ylim([0.15 0.7])
    ylabel('Admittance amplitude')
    xlim([datetime(2015,01,01) end_date])
    ax = gca;
    ax.FontSize = 14;
    legend('location','northeast')
    set(gcf,'Position',[100 100 1500 600])
    saveas(p,'./figs/SF_FR_M2_admittance.png')
end

%% wavelet "energy" plot
p=figure();
yyaxis right
plot(abs(constits.M2.filter),'k--','linewidth',2,'DisplayName','M_2 wavelet (abs; R)')
ylabel('Magnitude')
ylim([5.3e-7 1.5e-4])
yyaxis left
plot(real(constits.M2.filter),'xb-','linewidth',1,'DisplayName','M_2 wavelet (real; L)')
hold on
plot(imag(constits.M2.filter),'sr-','linewidth',1,'DisplayName','M_2 wavelet (imag; L)')
hold off
xline(85,'linewidth',2,'color',[0 0 0]+0.5,'HandleVisibility','off')
xline(252,'linewidth',2,'color',[0 0 0]+0.5,'DisplayName','80% energy window')
ylim([-1.75e-4 2e-4])
xlim([1 length(constits.M2.filter)])
ylabel('Amplitude')
xlabel('Index no.')
legend('location','northeast')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.FontSize = 14;
set(gcf,'Position',[100 100 1500 600])
saveas(p,'./figs/M2_wavelet.png')

%% two dynamically inferred ONLY constituents (P1 and K2) for False River
% DI AST
figure(); p=tiledlayout(2,1);

% P1 dyn inf

ax1=nexttile;
plot(constits.decTimesAll,constits.dynamicInference.P1.sixMoAmp(1:end-2),'linewidth',2,'DisplayName','P_1 CWT\_Multi (6 month filter)')
hold(ax1,'on')
plot(constits.decTimesAll,constits.dynamicInference.P1.shortAmps,'linewidth',2,'DisplayName','P_1 CWT\_Multi (dynamically inferred 2 weeks filter)')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
title('P_1 dynamic inference', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
legend('location','northeast')
grid on
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(constits.decTimesAll,constits.dynamicInference.K2.sixMoAmp(1:end-2),'linewidth',2,'DisplayName','K_2 CWT\_Multi (6 month filter)')
hold(ax2,'on')
plot(constits.decTimesAll,constits.dynamicInference.K2.shortAmps,'linewidth',2,'DisplayName','K_2 CWT\_Multi (dynamically inferred 2 weeks filter)')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
grid on
legend('location','northeast')
title('K_2 dynamic inference', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

set(gcf,'Position',[100 100 2500 1750])
saveas(p,'./figs/FR_dynamic_inference_D1_amps.png')

%% species and river flow results
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2,'color','black')
hold(ax1,'on')
qw{1}=plot(nan,'color','black','linewidth',2);
qw{2}=plot(nan,'color',[0 0.4470 0.7410],'linewidth',2); qw{3}=plot(nan,'color',[0.85 0.325 0.098],'linewidth',2);
qw{4}=plot(nan,'color',[0.929 0.694 0.125],'linewidth',2);qw{5}=plot(nan,'color',[0.494,0.184,0.556],'linewidth',2);
qw{6}=plot(nan,'r-','linewidth',2); qw{7}=plot(nan,'g-','linewidth',2);
hold(ax1,'off')
xlim([datetime(2015,01,01) end_date])
title('River flow', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylabel('NDOI discharge (m^3/s)')
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

legend([qw{:}], {'River flow','D_1','D_2','D_3','D_4','construction start','flood start'},'Location','northeast')


% xlim([constits.decTimesAll(1) constits.decTimesAll(end)])

ax2=nexttile;
plot(adSpecies.decTimesAll,adSpecies.d1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax2,'on')
plot(adSpecies.decTimesAll,adSpecies.d2.amps,'linewidth',2,'DisplayName','D_2')
plot(adSpecies.decTimesAll,adSpecies.d3.amps,'linewidth',2,'DisplayName','D_3')
plot(adSpecies.decTimesAll,adSpecies.d4.amps,'linewidth',2,'DisplayName','D_4')
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2,'color','black','DisplayName','River flow')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax2,'off')
ylabel('Amplitude (m)')
grid on
% legend('Location','northoutside') %'Position',[0.8,1.8,0.4,0.3])  %'location','northeast')
title('False River species results', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01) end_date])
ylim([0 0.45])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(adSpecies.decTimesAll,adRef.species.d1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax3,'on')
plot(species.decTimesAll,adRef.species.d2.amps,'linewidth',2,'DisplayName','D_2')
plot(species.decTimesAll,adRef.species.d3.amps,'linewidth',2,'DisplayName','D_3')
plot(species.decTimesAll,adRef.species.d4.amps,'linewidth',2,'DisplayName','D_4')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
ylabel('Amplitude (m)')
grid on
% legend('location','northeast')
title('SF Bay species results', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
set(gcf,'Position',[100 100 2400 2000])
plotAnn('(c)')

saveas(p,'./figs/SF_FR_flow_and_species.png')

%% residual spectrum
figure(); p=tiledlayout(2,1);

ax1=nexttile;
loglog(constits.input_freqs.^-1,constits.input_pwr,'r-','linewidth',1,'DisplayName','input data')
hold(ax1,'on')
plot(constits.fband.freqs.^-1,constits.fband.pwrs,'k-','linewidth',2,'DisplayName','CWT\_Multi residual')
% plot(constits.recon_freqs.^-1,constits.recon_pwr,'b-','linewidth',0.5,'DisplayName','recon data')
plot(OUT.resid_freqs.^-1,OUT.resid_pwr,'b-','linewidth',2,'DisplayName','UTide residual')
hold(ax1,'off')
grid on
xlim([0 5])
ylim auto
yma = max(constits.fband.pwrs)*1.2;
ymi = min(constits.fband.pwrs);
set(gca,'xdir','reverse')

legend('Location','northeast')
% axis tight
title('False River spectra', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
ax = gca;
ax.FontSize = 14;
plotAnnSpec('(a)')

ax2=nexttile;
loglog(constits.input_freqs.^-1,SFconstits.input_pwr,'r-','linewidth',1,'DisplayName','input data')
hold(ax2,'on')
plot(constits.fband.freqs.^-1,SFconstits.fband.pwrs,'k-','linewidth',2,'DisplayName','CWT\_Multi residual')
% plot(constits.recon_freqs.^-1,constits.recon_pwr,'b-','linewidth',0.5,'DisplayName','recon data')
plot(OUT.resid_freqs.^-1,OUTSF.resid_pwr,'b-','linewidth',2,'DisplayName','UTide residual')
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
title('SF Bay spectra', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
ax = gca;
ax.FontSize = 14;

ylabel(p,'PSD [m^2 cpd^{-1}]','FontSize',18)
set(gcf,'Position',[100 100 1800 1000])
% saveas(p,'./figs/SF_FR_spectra.png')

%% compare reconstruction & recon statistics
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(constits.alltimes,constits.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax1,'on')
plot(constits.alltimes,OUT.wlAll,'DisplayName','UTide reconstruction')
hold(ax1,'off')
grid on
legend('location','northeast')
ylabel('Water level [m]')
title('False River reconstruction', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(constits.alltimes,constits.reconAll - dat4.wl.','DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.0 1.0])
yyaxis right
plot(OUT.dtimesTrim,OUT.resid,'DisplayName','UTide error')
ylim([-0.5 1.5])
grid on
legend('location','northeast')
xlim([datetime(2015,01,01) end_date])
title('False River reconstruction error', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(constits.alltimes,SFconstits.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax3,'on')
plot(OUTSF.dtimes,OUTSF.wlAll,'DisplayName','UTide reconstruction')
hold(ax3,'off')
grid on
legend('location','northeast')
title('SF Bay reconstruction', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
ylabel('Water level [m]')
xlim([datetime(2015,01,01) end_date])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(constits.alltimes,SFconstits.reconAll-dat1.wl.','DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.0 1.0])
yyaxis right
plot(OUTSF.dtimesTrim,OUTSF.resid,'DisplayName','UTide error')
ylim([-0.5 1.5])
grid on
legend('location','northeast')
title('SF Bay reconstruction error', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
ax = gca;
ax.FontSize = 14;
xlim([datetime(2015,01,01) end_date])
plotAnn('(d)')

set(gcf,'Position',[100 100 2500 2000])
% saveas(p,'./figs/SF_FR_UTide_CWT_reconstruction.png')

%% FALSE RIVER DIURNAL AMPLITUDE
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(dayFlow.dtime,dayFlow.NDOI,'linewidth',2)
xlim([datetime(2015,01,01) end_date])
title('River flow', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylabel('NDOI discharge (m^3/s)')
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax1=nexttile;
plot(cd.datetimes,cd.Q1.amp,'color',blue,'linewidth',2,'DisplayName','Q_1 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
plotCI(constits.decTimesAll(start_date_ind:end),constits.Q1.amps(start_date_ind:end),constits.fband.ci(3),2)
hold(ax1,'off')
legend('location','northeast')
title('False River Q_1', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')


ax2=nexttile;
plot(cd.datetimes,cd.O1.amp,'color',blue,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.O1.amps(start_date_ind:end),constits.fband.ci(4),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River O_1', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

% K1
ax3=nexttile;
plot(cd.datetimes,cd.K1.amp,'color',blue,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi (original 2 weeks reponse)')
plotCI(constits.decTimesAll(start_date_ind:end),constits.K1.amps(start_date_ind:end),constits.fband.ci(6),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('False River K_1', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
legend('location','northeast')
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/FR_UTide_CWT_D1_amps.png')

%% FALSE RIVER SEMIDIURNAL AMPLITUDE

figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(cd.datetimes,cd.M2.amp,'color',blue,'linewidth',2,'DisplayName','M_2 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.M2.amps(start_date_ind:end),constits.fband.ci(13),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('False River M_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cd.N2.amp,'color',blue,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.N2.amps(start_date_ind:end),constits.fband.ci(12),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River N_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,cd.S2.amp,'color',blue,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.S2.amps(start_date_ind:end),constits.fband.ci(15),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('False River S_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,cd.L2.amp,'color',blue,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.L2.amps,'color',red,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.L2.amps(start_date_ind:end),constits.fband.ci(14),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('False River L_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
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
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.MK3.amps(start_date_ind:end),constits.fband.ci(19),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('False River MK_3', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

% plot results (M4)
ax2=nexttile;
plot(cd.datetimes,cd.M4.amp,'color',blue,'linewidth',2,'DisplayName','M_4 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.M4.amps(start_date_ind:end),constits.fband.ci(22),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River M_4', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
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
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.Q1.amps(start_date_ind:end),SFconstits.fband.ci(3),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('SF Bay Q_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cdSF.O1.amp,'color',blue,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.O1.amps(start_date_ind:end),SFconstits.fband.ci(4),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay O_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1
ax3=nexttile;
plot(cd.datetimes,cdSF.K1.amp,'color',blue,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi (original 2 weeks reponse)')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.K1.amps(start_date_ind:end),SFconstits.fband.ci(6),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('SF Bay K_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
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
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M2.amps(start_date_ind:end),SFconstits.fband.ci(13),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('SF Bay M_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,cdSF.N2.amp,'color',blue,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.N2.amps(start_date_ind:end),SFconstits.fband.ci(12),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay N_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,cdSF.S2.amp,'color',blue,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
axis tight
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.S2.amps(start_date_ind:end),SFconstits.fband.ci(15),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('SF Bay S_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,cdSF.L2.amp,'color',blue,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.L2.amps,'color',red,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.L2.amps(start_date_ind:end),SFconstits.fband.ci(14),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('SF Bay L_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
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
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.MK3.amps(start_date_ind:end),SFconstits.fband.ci(19),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('SF Bay MK_3', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

% plot results (M4)
ax2=nexttile;
plot(cd.datetimes,cdSF.M4.amp,'color',blue,'linewidth',2,'DisplayName','M_4 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M4.amps(start_date_ind:end),SFconstits.fband.ci(22),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay M_4', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
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
plot(cd.datetimes,OUTSF.Q1.phases_alt,'linewidth',2,'DisplayName','Q_1 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.Q1.phases(start_date_ind:end),SFconstits.fband.ci_ph(3),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('SF Bay Q_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,OUTSF.O1.phases_alt,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.O1.phases(start_date_ind:end),SFconstits.fband.ci_ph(4),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay O_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1
ax3=nexttile;
plot(cd.datetimes,OUTSF.K1.phases_alt,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi (original 2 weeks reponse)')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.K1.phases(start_date_ind:end),SFconstits.fband.ci_ph(6),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('SF Bay K_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
legend('location','northeast')
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/SF_UTide_CWT_D1_phases.png')

%% SF BAY SEMIDIURNAL PHASES
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(cd.datetimes,OUTSF.M2.phases_alt,'linewidth',2,'DisplayName','M_2 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.M2.phases,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.M2.phases(start_date_ind:end),SFconstits.fband.ci_ph(13),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.phases,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('SF Bay M_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 380])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,OUTSF.N2.phases_alt,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.N2.phases(start_date_ind:end),SFconstits.fband.ci_ph(12),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.phases,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('SF Bay N_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,OUTSF.S2.phases_alt,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.S2.phases(start_date_ind:end),SFconstits.fband.ci_ph(15),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('SF Bay S_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,OUTSF.L2.phases_alt,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01) end_date])
plot(SFconstits.decTimesAll,SFconstits.L2.phases,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(SFconstits.decTimesAll(start_date_ind:end),SFconstits.L2.phases(start_date_ind:end),SFconstits.fband.ci_ph(14),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('L_2 phase', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
grid on
ylim([0 440])
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/SF_UTide_CWT_D2_phases.png')



%% FALSE RIVER DIURNAL PHASES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(cd.datetimes,OUT.Q1.phases_alt,'linewidth',2,'DisplayName','Q_1 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.Q1.phases(start_date_ind:end),constits.fband.ci_ph(3),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax1,'off')
legend('location','northeast')
title('False River Q_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylim([0 420])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,OUT.O1.phases_alt,'linewidth',2,'DisplayName','O_1 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.O1.phases(start_date_ind:end),constits.fband.ci_ph(4),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.phases,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('False River O_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% K1 
ax3=nexttile;
plot(cd.datetimes,OUT.K1.phases_alt,'linewidth',2,'DisplayName','K_1 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi (original 2 weeks reponse)')
plotCI(constits.decTimesAll(start_date_ind:end),constits.K1.phases(start_date_ind:end),constits.fband.ci_ph(6),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
title('False River K_1', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
legend('location','northeast')
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/FR_UTide_CWT_D1_phases.png')


%% FALSE RIVER SEMIDIURNAL PHASES

figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(cd.datetimes,OUT.M2.phases_alt,'linewidth',2,'DisplayName','M_2 UTide')
hold(ax1,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.M2.phases,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.M2.phases(start_date_ind:end),constits.fband.ci_ph(13),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.phases,'r-','linewidth',2)
hold(ax1,'off')
legend('location','northeast')
title('False River M_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(cd.datetimes,OUT.N2.phases_alt,'linewidth',2,'DisplayName','N_2 UTide')
hold(ax2,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.N2.phases(start_date_ind:end),constits.fband.ci_ph(12),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
% plot(species.dectimes,species.d2.phases,'r-','linewidth',2)
hold(ax2,'off')
legend('location','northeast')
title('N_2 phase', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(cd.datetimes,OUT.S2.phases_alt,'linewidth',2,'DisplayName','S_2 UTide')
hold(ax3,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.S2.phases(start_date_ind:end),constits.fband.ci_ph(15),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax3,'off')
legend('location','northeast')
title('False River S_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
set(gca,'xtick',[])
grid on
ylim([0 440])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(cd.datetimes,OUT.L2.phases_alt,'linewidth',2,'DisplayName','L_2 UTide')
hold(ax4,'on')
xlim([datetime(2015,01,01) end_date])
plot(constits.decTimesAll,constits.L2.phases,'linewidth',2,'DisplayName','L_2 CWT\_Multi')
plotCI(constits.decTimesAll(start_date_ind:end),constits.L2.phases(start_date_ind:end),constits.fband.ci_ph(14),2)
xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
xline(datetime(2017,01,07),'g-','linewidth',2,'DisplayName','flood start')
hold(ax4,'off')
legend('location','northeast')
title('False River L_2', 'Units', 'normalized', 'Position', [0.5, 0.75, 0])
grid on
ylim([0 360])
ax = gca;
ax.FontSize = 14;
plotAnn('(d)')

ylabel(p,'Phase (deg)','FontSize',18)
set(gcf,'Position',[100 0 2500 1400])
saveas(p,'./figs/FR_UTide_CWT_D2_phases.png')

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
