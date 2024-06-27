% Comparing CWT\_Multi results to NS\_Tide results for Columbia River stations
% clear

[vanDIConstits,vanDISpecies,vanDIRef,~] = cwtMulti(Van.dates,Van.wl,'dynamicInference','alt_phase');


%% load data
addpath('../cwtMultiPackage/CWT_Multi')

load('./fromPascal/NS_TIDE_comparisons2/NS_TIDE_tests_LCR2.mat')

%% add timezones
Van.dates.TimeZone = 'America/Los_Angeles';
Ast.dates.TimeZone = 'America/Los_Angeles';

%% Assigning all NS\_Tide outputs to two structs
AstNS.t = t;
AstNS.dtime = datetime(t,'ConvertFrom','datenum');

VanNS.t = t;
VanNS.dtime = datetime(t,'ConvertFrom','datenum');

dt = (t(2)-t(1))*24*3600;

% align datasets to common start/end
startDate = datetime(2010,04,05);
endDate = datetime(2013,06,28);

AstNS.start = find(AstNS.dtime==startDate);
AstNS.end = find(AstNS.dtime==endDate);

AstNS.dtime = AstNS.dtime(AstNS.start:AstNS.end);
AstNS.wlOrig = wl_Ast(AstNS.start:AstNS.end);

AstNS.recon = poutAst(AstNS.start:AstNS.end);
nonNanInd = find(~isnan(AstNS.recon));
AstNS.resid = AstNS.wlOrig(nonNanInd) - AstNS.recon(nonNanInd);
AstNS.rmse = mean(sqrt((AstNS.wlOrig(nonNanInd) - AstNS.recon(nonNanInd)).^2));

VanNS.start = find(VanNS.dtime==startDate);
VanNS.end = find(VanNS.dtime==endDate);

VanNS.dtime = VanNS.dtime(VanNS.start:VanNS.end);
VanNS.wlOrig = wl_Van(VanNS.start:VanNS.end);

VanNS.recon = poutVan(VanNS.start:VanNS.end);
nonNanInd = find(~isnan(VanNS.recon));
VanNS.resid = VanNS.wlOrig(nonNanInd) - VanNS.recon(nonNanInd);
VanNS.rmse = mean(sqrt((VanNS.wlOrig(nonNanInd) - VanNS.recon(nonNanInd)).^2));

% finding names Pascal's way
d1A = strmatch('Q1  ',nameuAst) - find(fuAst<.03,1,'last');
d2A = strmatch('O1  ',nameuAst) - find(fuAst<.03,1,'last');
d3A = strmatch('P1  ',nameuAst) - find(fuAst<.03,1,'last');
d4A = strmatch('K1  ',nameuAst) - find(fuAst<.03,1,'last');

d1V = strmatch('Q1  ',nameuVan) - find(fuVan<.03,1,'last');
d2V = strmatch('O1  ',nameuVan) - find(fuVan<.03,1,'last');
d3V = strmatch('P1  ',nameuVan) - find(fuVan<.03,1,'last');
d4V = strmatch('K1  ',nameuVan) - find(fuVan<.03,1,'last');

%---Semi-Diurnals
s1A = strmatch('N2  ',nameuAst) - find(fuAst<.07,1,'last');
s2A = strmatch('M2  ',nameuAst) - find(fuAst<.07,1,'last');
s3A = strmatch('S2  ',nameuAst) - find(fuAst<.07,1,'last');

s1V = strmatch('N2  ',nameuVan) - find(fuVan<.07,1,'last');
s2V = strmatch('M2  ',nameuVan) - find(fuVan<.07,1,'last');
s3V = strmatch('S2  ',nameuVan) - find(fuVan<.07,1,'last');

%---MK3, M4
o1V = strmatch('MK3 ',nameuVan) - find(fuVan<.10,1,'last');
o2V = strmatch('M4  ',nameuVan) - find(fuVan<.14,1,'last');

o1A = strmatch('MK3 ',nameuAst) - find(fuAst<.10,1,'last');
o2A = strmatch('M4  ',nameuAst) - find(fuAst<.14,1,'last');

% amplitudes 
AstNS.Q1 = tidecon_tsAst{1,3}{2,1}(d1A,AstNS.start:AstNS.end);
AstNS.O1 = tidecon_tsAst{1,3}{2,1}(d2A,AstNS.start:AstNS.end);
AstNS.P1 = tidecon_tsAst{1,3}{2,1}(d3A,AstNS.start:AstNS.end);
AstNS.K1 = tidecon_tsAst{1,3}{2,1}(d4A,AstNS.start:AstNS.end);

VanNS.Q1 = tidecon_tsVan{1,3}{2,1}(d1V,VanNS.start:VanNS.end);
VanNS.O1 = tidecon_tsVan{1,3}{2,1}(d2V,VanNS.start:VanNS.end);
VanNS.K1 = tidecon_tsVan{1,3}{2,1}(d4V,VanNS.start:VanNS.end);

AstNS.N2 = tidecon_tsAst{1,3}{3,1}(s1A,AstNS.start:AstNS.end);
AstNS.M2 = tidecon_tsAst{1,3}{3,1}(s2A,AstNS.start:AstNS.end);
AstNS.S2 = tidecon_tsAst{1,3}{3,1}(s3A,AstNS.start:AstNS.end);

VanNS.N2 = tidecon_tsVan{1,3}{3,1}(s1V,VanNS.start:VanNS.end);
VanNS.M2 = tidecon_tsVan{1,3}{3,1}(s2V,VanNS.start:VanNS.end);
VanNS.S2 = tidecon_tsVan{1,3}{3,1}(s3V,VanNS.start:VanNS.end);

AstNS.MK3 = tidecon_tsAst{1,3}{4,1}(o1A,AstNS.start:AstNS.end);
AstNS.M4 = tidecon_tsAst{1,3}{5,1}(o2A,AstNS.start:AstNS.end);

VanNS.MK3 = tidecon_tsVan{1,3}{4,1}(o1V,VanNS.start:VanNS.end);
VanNS.M4 = tidecon_tsVan{1,3}{5,1}(o2V,VanNS.start:VanNS.end);

% phases
AstNS.Q1ph = tidecon_tsAst{1,5}{2,1}(d1A,AstNS.start:AstNS.end);
AstNS.O1ph = tidecon_tsAst{1,5}{2,1}(d2A,AstNS.start:AstNS.end);
AstNS.P1ph = tidecon_tsAst{1,5}{2,1}(d3A,AstNS.start:AstNS.end);
AstNS.K1ph = tidecon_tsAst{1,5}{2,1}(d4A,AstNS.start:AstNS.end);

VanNS.Q1ph = tidecon_tsVan{1,5}{2,1}(d1V,VanNS.start:VanNS.end);
VanNS.O1ph = tidecon_tsVan{1,5}{2,1}(d2V,VanNS.start:VanNS.end);
VanNS.K1ph = tidecon_tsVan{1,5}{2,1}(d4V,VanNS.start:VanNS.end);

AstNS.N2ph = tidecon_tsAst{1,5}{3,1}(s1A,AstNS.start:AstNS.end);
AstNS.M2ph = tidecon_tsAst{1,5}{3,1}(s2A,AstNS.start:AstNS.end);
AstNS.S2ph = tidecon_tsAst{1,5}{3,1}(s3A,AstNS.start:AstNS.end);

VanNS.N2ph = tidecon_tsVan{1,5}{3,1}(s1V,VanNS.start:VanNS.end);
VanNS.M2ph = tidecon_tsVan{1,5}{3,1}(s2V,VanNS.start:VanNS.end);
VanNS.S2ph = tidecon_tsVan{1,5}{3,1}(s3V,VanNS.start:VanNS.end);

VanNS.MK3ph = tidecon_tsVan{1,5}{4,1}(o1V,VanNS.start:VanNS.end);
VanNS.M4ph  = tidecon_tsVan{1,5}{5,1}(o2V,VanNS.start:VanNS.end);

% admittances
AdmNS.Q1 = VanNS.Q1./AstNS.Q1;
AdmNS.O1 = VanNS.O1./AstNS.O1;
AdmNS.K1 = VanNS.K1./AstNS.K1;

AdmNS.N2 = VanNS.N2./AstNS.N2;
AdmNS.M2 = VanNS.M2./AstNS.M2;
AdmNS.S2 = VanNS.S2./AstNS.S2;

% copmuting NS\_Tide residual spectra for figures below
c = VanNS.resid;

bad_vals = find(isnan(c));
for i=1:10
    c(bad_vals) = 0;
    c = c-mean(c);
end

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

nfft = 2^12;

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/dt);

VanNS.resid_pwr=pxx/3600/24;
VanNS.resid_freqs=f*3600*24;      % frequency in cpd

c = AstNS.resid;

bad_vals = find(isnan(c));
for i=1:10
    c(bad_vals) = 0;
    c = c-mean(c);
end

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

nfft = 2^12;

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/dt);

AstNS.resid_pwr=pxx/3600/24;
AstNS.resid_freqs=f*3600*24;      % frequency in cpd

%% Align dates in dataset
load('./data/Astoria/astoria_wl.mat')
load('./data/Vancouver/vancouver_wl.mat')

% align datasets to common start/end
% startDate = datetime(2010,04,01);
% endDate = datetime(2013,07,01);

Van.start = find(Van.dates==startDate);
Van.end = find(Van.dates==endDate);

Van.dates = Van.dates(Van.start:Van.end);
Van.wl = Van.wl(Van.start:Van.end);

Ast.start = find(Ast.dates==startDate);
Ast.end = find(Ast.dates==endDate);

Ast.dates = Ast.dates(Ast.start:Ast.end);
Ast.wl = Ast.wl(Ast.start:Ast.end);

nan_ind = unique(sort([find(isnan(Van.wl));find(isnan(Ast.wl))]));

Van.wl(nan_ind) = NaN;
Ast.wl(nan_ind) = NaN;

% lon, lat for Astoria
lon = 123+46.1/60;            % W
lat = 46+12.4/60;             % N

dInd = floor(length(Van.dates)/2);
vInd = floor(length(VanNS.dtime)/2);
% ref_date = datetime(0,0,0,0,0,0) + (Van.dates(vInd)-datetime(1899,12,31,12,0,0));
ref_date = datetime(1899,12,31,12,0,0);
ref_date.TimeZone = 'UTC';

%% including SO3
coFiltLength    = [1081, 1081,1081,363,1081,363,1081,1081,1081,1081,1081,363,363,1081,1081,363,1081,1081,1081,363,1081,1081,363,1081,1081,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours

constitsD1Freqs = [29.0727, 28.00621204, 26.868350, 25.81933871, 24.84120241, 23.93447213, 23.09848146, 22.30608083, 21.5782].^-1;
constitsD1Names = ["Alp1","TwoQ1","Q1","O1","NO1","K1","J1","OO1","UPS1"];

constitsD2Freqs = [13.1272668, 12.8717576, 12.65834751, 12.4206012, 12.19162085, 12, 11.75452172].^-1;
constitsD2Names = ["Eps2","Mu2","N2","M2","L2","S2","Eta2"];   % Nu2 is number 10

constitsD3Freqs = [8.3863068, 8.280400802, 0.1220639878^-1, 8.177140247, 7.9936].^-1;
constitsD3Names = ["MO3","M3","SO3","MK3","SK3"];

constitsD4PlusFreqs = [6.2691737, 6.21030, 6.103339, 6.0, 4.93088021306, 4.140200399, 3.52964079728,...
    3.10515029954, 2.748563985947, 2.48412023963, 2.25054027184, 2.07010019969].^-1;
constitsD4PlusNames = ["MN4","M4","MS4","S4","MK5","M6","MK7","M8","MK9","M10","MK11","M12"];

coFreqs = {constitsD1Freqs,constitsD2Freqs,constitsD3Freqs,constitsD4PlusFreqs};
coNames = horzcat(constitsD1Names,constitsD2Names,constitsD3Names,constitsD4PlusNames);


%% run cwt routine on Vancouver/Astoria data

[constits,species,ref,admit] = cwtMulti(Van.dates,Van.wl,'performAdmittance',lon,lat,'refStation',Ast.wl,'rel_phase',ref_date); 

% Astoria analysis with dynamic inference
[astDIConstits,astDISpecies,astDIRef,~] = cwtMulti(Ast.dates,Ast.wl,'dynamicInference','alt_phase'); %,'coFreqs',coFreqs,coFiltLength,coNames);

% Vancouver analysis with dynamic inference
[vanDIConstits,vanDISpecies,vanDIRef,~] = cwtMulti(Van.dates,Van.wl,'dynamicInference','alt_phase');

%%
% must also do admittance relative to tidal potential at GMT to get phases
% relative to GMT (NOT WORKING [WHEN COMPARED TO NS_TIDE])
% [~,~,~,admitGMTVan] = cwtMulti(Van.dates,Van.wl,'performAdmittance',0,0); % tidal potential at lat=0, lon=0
% 
% [~,~,~,admitGMTAst] = cwtMulti(Ast.dates,Ast.wl,'performAdmittance',0,0); % tidal potential at lat=0, lon=0

%% admittance
% [constits,species,ref,admit] = cwtMulti(Van.dates,Van.wl,'performAdmittance',lon,lat); 

%% Perform UTide analysis on Vancouver
addpath('./utide_tool/')

Ast.datenums = datenum(Ast.dates);
Van.datenums = datenum(Van.dates);
Ast.dtime = Ast.dates;
Van.dtime = Van.dates;

Ast.lat = lat; Ast.lon = lon;
Van.lat = lat; Van.lon = lon;

% Define moving window parameters etc.
window = 30*1.0;         % window in days
incr = 20/24;%1/24;     % increment in days 
rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.8;     % specify what percentage of data in window must be not NaN

maxFLength = 5*24;   % same cutoff period as low-pass filtering done in CWT\_Multi

[cd,OUT] = cwt_utide(Van,window,incr,good_pct,rayleigh,maxFLength,[datetime(2010,10,01) datetime(2013,01,01)]);%16P1/17K1

% Nu2 stuff at Astoria
% [nu2constits,nu2species,nu2ref,nu2admit] = cwtMulti(Ast.dates,Ast.wl,'coFreqs',coFreqs,coFiltLength,coNames);

%% defining plotting colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];

%% a dynamic inference plot
p=figure();
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.K1.shortAmps,'color',blue,'linewidth',2,'DisplayName','K_1 inferred')
hold on
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.P1.shortAmps,'color',red,'linewidth',2,'DisplayName','P_1 inferred')
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.S1.shortAmps,'color',yellow,'linewidth',2,'DisplayName','S_1 inferred')
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.K1.sixMoAmp,'--','color',blue,'linewidth',2,'DisplayName','K_1 (6 months filter)')
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.P1.sixMoAmp,'--','color',red,'linewidth',2,'DisplayName','P_1 (6 months filter)')
plot(vanDIConstits.decTimesAll,vanDIConstits.dynamicInference.S1.sixMoAmp,'--','color',yellow,'linewidth',2,'DisplayName','S_1 (6 months filter)')
hold off
grid on
ylabel('Amplitude (m)','FontSize',24)
title('Vancouver diurnal three-way dynamic inference', 'Units', 'normalized', 'Position', [0.5, 1.025, 0])
xlim([Van.dates(1) Van.dates(end)])
legend('Location','northeast')
ax = gca;
ax.FontSize = 18;
set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/Vancouver_D1_dynamic_inference.png')



%% Water level plots
figure(); p=tiledlayout(2,1);

ax1=nexttile;
plot(Ast.dates,Ast.wl,'color',blue,'linewidth',2,'DisplayName','water level (Astoria)')
grid on
legend('location','northeast')
title('Astoria water level', 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
set(gca,'xtick',[])
xlim([Ast.dates(1) Ast.dates(end)])
ylim([-1 4.25])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(Van.dates,Van.wl,'color',blue,'linewidth',2,'DisplayName','water level (Vancouver)')
grid on
legend('location','northeast')
title('Vancouver water level', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
% set(gca,'xtick',[])
xlim([Van.dates(1) Van.dates(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_Vancouver_water_levels.png')

%% ASTORIA SEMIDIURNAL AMPLITUDES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(AstNS.dtime,AstNS.M2,'color',blue,'linewidth',2,'DisplayName','M_2 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,astDIConstits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.M2.amps,astDIConstits.fband.ci(13),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Astoria M_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0]) 
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(AstNS.dtime,AstNS.N2,'color',blue,'linewidth',2,'DisplayName','N_2 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,astDIConstits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.N2.amps,astDIConstits.fband.ci(12),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Astoria N_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ylim([0.135 0.245])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')


ax3=nexttile;
plot(AstNS.dtime,AstNS.S2,'color',blue,'linewidth',2,'DisplayName','S_2 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,astDIConstits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.S2.amps,astDIConstits.fband.ci(15),2)
hold(ax3,'off')
grid on
legend('location','northeast')
title('Astoria S_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
axis tight
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_NS_CWT_D2_amps.png')


%% ASTORIA DIURNAL AMPLITUDES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(AstNS.dtime,AstNS.Q1,'color',blue,'linewidth',2,'DisplayName','Q_1 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,astDIConstits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.Q1.amps,astDIConstits.fband.ci(3),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Astoria Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(AstNS.dtime,AstNS.O1,'color',blue,'linewidth',2,'DisplayName','O_1 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,astDIConstits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.O1.amps,astDIConstits.fband.ci(4),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Astoria O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(AstNS.dtime,AstNS.K1,'color',blue,'linewidth',2,'DisplayName','K_1 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,astDIConstits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
plotCI(constits.decTimesAll,astDIConstits.K1.amps,astDIConstits.fband.ci(6),2)
% plot(constits.decTimesAll,astDIConstits.dynamicInference.K1.shortAmps,'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','K_1 CWT\_Multi (dynamically inferred 2 weeks amplitude)')
% plot(constits.decTimesAll,astDIConstits.dynamicInference.K1.sixMoAmp,':','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','K_1 CWT\_Multi (6 months filter)')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Astoria K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_NS_CWT_D1_amps.png')

%% VANCOUVER SEMIDIURNAL AMPLITUDES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(VanNS.dtime,VanNS.M2,'color',blue,'linewidth',2,'DisplayName','M_2 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,constits.M2.amps,'color',red,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
plotCI(constits.decTimesAll,constits.M2.amps,constits.fband.ci(13),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Vancouver M_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,VanNS.N2,'color',blue,'linewidth',2,'DisplayName','N_2 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,constits.N2.amps,'color',red,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
plotCI(constits.decTimesAll,constits.N2.amps,constits.fband.ci(12),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver N_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(VanNS.dtime,VanNS.S2,'color',blue,'linewidth',2,'DisplayName','S_2 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,constits.S2.amps,'color',red,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
plotCI(constits.decTimesAll,constits.S2.amps,constits.fband.ci(15),2)
% plot(constits.decTimesAll,vanDIConstits.dynamicInference.S2.shortAmps,'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi (dynamically inferred 2 weeks amplitude)')
% plot(constits.decTimesAll,vanDIConstits.dynamicInference.S2.sixMoAmp,':','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','S_2 CWT\_Multi (6 months filter)')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Vancouver S_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
plotAnn('(c)')


set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Vancouver_NS_CWT_D2_amps.png')

%% VANCOUVER DIURNAL AMPLITUDES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(VanNS.dtime,VanNS.Q1,'color',blue,'linewidth',2,'DisplayName','Q_1 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,constits.Q1.amps,'color',red,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
plotCI(constits.decTimesAll,constits.Q1.amps,constits.fband.ci(3),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Vancouver Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,VanNS.O1,'color',blue,'linewidth',2,'DisplayName','O_1 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,constits.O1.amps,'color',red,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
plotCI(constits.decTimesAll,constits.O1.amps,constits.fband.ci(4),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(VanNS.dtime,VanNS.K1,'color',blue,'linewidth',2,'DisplayName','K_1 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,constits.K1.amps,'color',red,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
plotCI(constits.decTimesAll,constits.K1.amps,constits.fband.ci(6),2)
% plot(constits.decTimesAll,vanDIConstits.dynamicInference.K1.shortAmps,'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','K_1 CWT\_Multi (dynamically inferred 2 weeks amplitude)')
% plot(constits.decTimesAll,vanDIConstits.dynamicInference.K1.sixMoAmp,':','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','K_1 CWT\_Multi (6 months filter)')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Vancouver K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Vancouver_NS_CWT_D1_amps.png')


%% ASTORIA OVERTIDE AMPLITUDES
figure(); p=tiledlayout(2,1);

ax1=nexttile;
plot(AstNS.dtime,AstNS.MK3,'color',blue,'linewidth',2,'DisplayName','MK_3 NS\_Tide')
hold(ax1,'on')
plot(astDIConstits.decTimesAll,astDIConstits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(astDIConstits.decTimesAll,astDIConstits.MK3.amps,astDIConstits.fband.ci(19),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Astoria MK_3', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([astDIConstits.decTimesAll(1) astDIConstits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,AstNS.M4,'color',blue,'linewidth',2,'DisplayName','M_4 NS\_Tide')
hold(ax2,'on')
plot(astDIConstits.decTimesAll,astDIConstits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(astDIConstits.decTimesAll,astDIConstits.M4.amps,astDIConstits.fband.ci(22),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Astoria M_4', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([astDIConstits.decTimesAll(1) astDIConstits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
plotAnn('(b)')

set(gcf,'Position',[100 100 2500 1750])
% saveas(p,'./figs/Astoria_NS_CWT_overtide_amps.png')

%% VANCOUVER OVERTIDE AMPLITUDES
figure(); p=tiledlayout(2,1);

ax1=nexttile;
plot(VanNS.dtime,VanNS.MK3,'color',blue,'linewidth',2,'DisplayName','MK_3 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,constits.MK3.amps,'color',red,'linewidth',2,'DisplayName','MK_3 CWT\_Multi')
plotCI(constits.decTimesAll,constits.MK3.amps,constits.fband.ci(19),2)
hold(ax1,'off')
grid on
legend('location','northeast')
title('Vancouver MK_3', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,VanNS.M4,'color',blue,'linewidth',2,'DisplayName','M_4 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,constits.M4.amps,'color',red,'linewidth',2,'DisplayName','M_4 CWT\_Multi')
plotCI(constits.decTimesAll,constits.M4.amps,constits.fband.ci(22),2)
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver M_4', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
plotAnn('(b)')

set(gcf,'Position',[100 100 2500 1750])
saveas(p,'./figs/Vancouver_NS_CWT_overtide_amps.png')
%% SEMMIDIURNAL ADMITTANCES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(VanNS.dtime,AdmNS.M2,'linewidth',2,'DisplayName','M_2 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,admit.constits.M2.amps,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
hold(ax1,'off')
grid on
legend('location','northeast')
title('Vancouver M_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,AdmNS.N2,'linewidth',2,'DisplayName','N_2 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,admit.constits.N2.amps,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver N_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(VanNS.dtime,AdmNS.S2,'linewidth',2,'DisplayName','S_2 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,admit.constits.S2.amps,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Vancouver S_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Admittance amplitude','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_Vancouver_NS_CWT_D2_admittances.png')


%% DIURNAL ADMITTANCES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(VanNS.dtime,AdmNS.Q1,'linewidth',2,'DisplayName','Q_1 NS\_Tide')
hold(ax1,'on')
plot(constits.decTimesAll,admit.constits.Q1.amps,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
hold(ax1,'off')
grid on
legend('location','northeast')
title('Vancouver Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(VanNS.dtime,AdmNS.O1,'linewidth',2,'DisplayName','O_1 NS\_Tide')
hold(ax2,'on')
plot(constits.decTimesAll,admit.constits.O1.amps,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(VanNS.dtime,AdmNS.K1,'linewidth',2,'DisplayName','K_1 NS\_Tide')
hold(ax3,'on')
plot(constits.decTimesAll,admit.constits.K1.amps,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Vancouver K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Admittance amplitude','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_Vancouver_NS_CWT_D1_admittances.png')

%% river flow and species results
figure(); p=tiledlayout(3,1);

ax1=nexttile;
plot(AstNS.dtime,Q_Bon(AstNS.start:AstNS.end),VanNS.dtime,Q_Wil(VanNS.start:VanNS.end),'linewidth',2);ylabel('Discharge [m^3/s]');grid on;
title('River flow', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
legend('Columbia at Bonneville Dam','Willamette at Portland');
ax = gca;
ax.FontSize = 14;
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
set(gca,'xtick',[])
plotAnn('(a)')

ax2=nexttile;
plot(species.decTimesAll,astDISpecies.D1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax2,'on')
plot(species.decTimesAll,astDISpecies.D2.amps,'linewidth',2,'DisplayName','D_2')
plot(species.decTimesAll,astDISpecies.D3.amps,'linewidth',2,'DisplayName','D_3')
plot(species.decTimesAll,astDISpecies.D4.amps,'linewidth',2,'DisplayName','D_4')
hold(ax2,'off')
ylabel('Amplitude (m)')
grid on
legend('location','northeast')
title('Astoria species results', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(species.decTimesAll,vanDISpecies.D1.amps,'linewidth',2,'DisplayName','D_1')
hold(ax3,'on')
plot(species.decTimesAll,vanDISpecies.D2.amps,'linewidth',2,'DisplayName','D_2')
plot(species.decTimesAll,vanDISpecies.D3.amps,'linewidth',2,'DisplayName','D_3')
plot(species.decTimesAll,vanDISpecies.D4.amps,'linewidth',2,'DisplayName','D_4')
hold(ax3,'off')
ylabel('Amplitude (m)')
grid on
legend('location','northeast')
title('Vancouver species results', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

set(gcf,'Position',[100 100 2400 2000])
saveas(p,'./figs/Astoria_Vancouver_flow_and_species.png')

%% compare reconstruction & recon statistics
figure(); p=tiledlayout(4,1);

ax1=nexttile;
plot(constits.alltimes,astDIConstits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax1,'on')
plot(AstNS.dtime,AstNS.recon,'DisplayName','NS\_Tide reconstruction')
hold(ax1,'off')
grid on
legend('location','northeast')
ylabel('Water level [m]')
title('Astoria reconstruction', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.alltimes(1) constits.alltimes(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
plot(constits.alltimes,Ast.wl.' - astDIConstits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.5 0.5])
yyaxis right
plot(AstNS.dtime,AstNS.wlOrig-AstNS.recon,'DisplayName','NS\_Tide error')
ylim([-1 1.5])
grid on
legend('location','northeast')
xlim([constits.alltimes(1) constits.alltimes(end)])
title('Astoria reconstruction error', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
plot(constits.alltimes,constits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent reconstruction')
hold(ax3,'on')
plot(AstNS.dtime,VanNS.recon,'DisplayName','NS\_Tide reconstruction')
hold(ax3,'off')
grid on
legend('location','northeast')
title('Vancouver reconstruction', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ylabel('Water level [m]')
xlim([constits.alltimes(1) constits.alltimes(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

ax4=nexttile;
plot(constits.alltimes,Van.wl.' - constits.reconstruction.reconAll,'DisplayName','CWT\_Multi constituent error')
ylabel('Water level error [m]')
ylim([-1.5 0.5])
yyaxis right
plot(VanNS.dtime,VanNS.wlOrig-VanNS.recon,'DisplayName','NS\_Tide error')
ylim([-1 1.5])
grid on
legend('location','northeast')
title('Vancouver reconstruction error', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
ax = gca;
ax.FontSize = 14;
xlim([constits.alltimes(1) constits.alltimes(end)])
plotAnn('(d)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_Vancouver_NS_CWT_reconstruction.png')



%% power spectrum
figure(); p=tiledlayout(2,1);
sp_names = ["D_1","D_2","D_3","D_4","D_5","D_6","D_7","D_8","D_9","D_{10}","D_{11}","D_{12}"];

ax1=nexttile;
loglog(constits.input_freqs.^-1,constits.input_pwr,'r-','linewidth',1,'DisplayName','input data')
hold(ax1,'on')
plot(VanNS.resid_freqs.^-1,VanNS.resid_pwr,'b-','linewidth',2,'DisplayName','NS\_Tide residual')
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
set(gca,'xtick',[])
legend('Location','southwest')
% axis tight
title('Vancouver spectra', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
ax = gca;
ax.FontSize = 14;
plotAnnSpec('(a)')

ax2=nexttile;
loglog(constits.input_freqs.^-1,astDIConstits.input_pwr,'r-','linewidth',1,'DisplayName','input data')
hold(ax2,'on')
plot(VanNS.resid_freqs.^-1,AstNS.resid_pwr,'b-','linewidth',2,'DisplayName','NS\_Tide residual')
plot(constits.fband.freqs.^-1,astDIConstits.fband.pwrs,'k-','linewidth',2,'DisplayName','CWT\_Multi residual')
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
legend('Location','southwest')
% axis tight
title('Astoria spectra', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
ax = gca;
ax.FontSize = 14;

ylabel(p,'PSD [m^2 cpd^{-1}]','FontSize',18)
set(gcf,'Position',[100 100 1800 1000])
saveas(p,'./figs/Ast_Van_spectra.png')

%% ASTORIA SEMIDIURNAL PHASES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
% plot(AstNS.dtime,AstNS.M2ph,'linewidth',2,'DisplayName','M_2 NS\_Tide')
% hold(ax1,'on')
% plot(constits.decTimesAll,astDIConstits.M2.phases,'linewidth',2,'DisplayName','M_2 CWT\_Multi')
% hold(ax1,'off')
plotPhase(constits.decTimesAll,astDIConstits.M2.phases,AstNS.dtime,AstNS.M2ph,ax1,astDIConstits.fband.ci_ph(12),2,'M_2',50)
grid on
legend('location','northeast')
title('Astoria M_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(AstNS.dtime,AstNS.N2ph,'linewidth',2,'DisplayName','N_2 NS\_Tide')
% hold(ax2,'on')
% plot(constits.decTimesAll,astDIConstits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% hold(ax2,'off')
plotPhase(constits.decTimesAll,astDIConstits.N2.phases,AstNS.dtime,AstNS.N2ph,ax2,astDIConstits.fband.ci_ph(13),2,'N_2',-145)
grid on
legend('location','northeast')
title('Astoria N_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
ylim([-160 400])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
% plot(AstNS.dtime,AstNS.S2ph,'linewidth',2,'DisplayName','S_2 NS\_Tide')
% hold(ax3,'on')
% plot(constits.decTimesAll,astDIConstits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% hold(ax3,'off')
plotPhase(constits.decTimesAll,astDIConstits.S2.phases,AstNS.dtime,AstNS.S2ph,ax3,astDIConstits.fband.ci_ph(15),2,'S_2',-145)
grid on
legend('location','northeast')
title('Astoria S_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Phase (deg)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_NS_CWT_D2_phases.png')

%% ASTORIA DIURNAL PHASES
%% CHANGING NW TO 1500 FIXED N2 PLOT BUT MADE K1 PLOT BAD; NOT SAVING THIS PLOT FOR NOW
figure(); p=tiledlayout(3,1);

ax1=nexttile;
% plot(AstNS.dtime,AstNS.Q1ph,'linewidth',2,'DisplayName','Q_1 NS\_Tide')
% hold(ax1,'on')
% plot(constits.decTimesAll,astDIConstits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% hold(ax1,'off')
plotPhase(constits.decTimesAll,astDIConstits.Q1.phases,AstNS.dtime,AstNS.Q1ph,ax1,astDIConstits.fband.ci_ph(3),2,'Q_1',-145)
grid on
legend('location','northeast')
title('Astoria Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(AstNS.dtime,AstNS.O1ph,'linewidth',2,'DisplayName','O_1 NS\_Tide')
% hold(ax2,'on')
% plot(constits.decTimesAll,astDIConstits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% hold(ax2,'off')
plotPhase(constits.decTimesAll,astDIConstits.O1.phases,AstNS.dtime,AstNS.O1ph,ax2,astDIConstits.fband.ci_ph(4),2,'O_1',-145)
grid on
legend('location','northeast')
title('Astoria O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
% plot(AstNS.dtime,AstNS.K1ph,'linewidth',2,'DisplayName','K_1 NS\_Tide')
% hold(ax3,'on')
% plot(constits.decTimesAll,astDIConstits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% hold(ax3,'off')
plotPhase(constits.decTimesAll,astDIConstits.K1.phases,AstNS.dtime,AstNS.K1ph,ax3,astDIConstits.fband.ci_ph(6),2,'K_1',-145)
grid on
legend('location','northeast')
title('Astoria K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Phase (deg)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Astoria_NS_CWT_D1_phases.png')



%% VANCOUVER SEMIDIURNAL PHASES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
% plot(VanNS.dtime,VanNS.M2ph,'linewidth',2,'DisplayName','M_2 NS\_Tide')
% hold(ax1,'on')
% plot(constits.decTimesAll,constits.M2.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% hold(ax1,'off')
plotPhase(constits.decTimesAll,constits.M2.phases,VanNS.dtime,VanNS.M2ph,ax1,constits.fband.ci_ph(13),2,'M_2',-145)
grid on
legend('location','northeast')
title('Vancouver M_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(VanNS.dtime,VanNS.N2ph,'linewidth',2,'DisplayName','N_2 NS\_Tide')
% hold(ax2,'on')
% plot(constits.decTimesAll,constits.N2.phases,'linewidth',2,'DisplayName','N_2 CWT\_Multi')
% hold(ax2,'off')
plotPhase(constits.decTimesAll,constits.N2.phases,VanNS.dtime,VanNS.N2ph,ax2,constits.fband.ci_ph(12),2,'N_2',-145)
grid on
legend('location','northeast')
title('Vancouver N_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')


ax3=nexttile;
% plot(VanNS.dtime,VanNS.S2ph,'linewidth',2,'DisplayName','S_2 NS\_Tide')
% hold(ax3,'on')
% plot(constits.decTimesAll,constits.S2.phases,'linewidth',2,'DisplayName','S_2 CWT\_Multi')
% hold(ax3,'off')
plotPhase(constits.decTimesAll,constits.S2.phases,VanNS.dtime,VanNS.S2ph,ax3,constits.fband.ci_ph(15),2,'S_2',-145)
grid on
legend('location','northeast')
title('Vancouver S_2', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel(p,'Phase (deg)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
saveas(p,'./figs/Vancouver_NS_CWT_D2_phases.png')




%% VANCOUVER DIURNAL PHASES
figure(); p=tiledlayout(3,1);

ax1=nexttile;
% plot(VanNS.dtime,VanNS.Q1ph,'linewidth',2,'DisplayName','Q_1 NS\_Tide')
% hold(ax1,'on')
% plot(constits.decTimesAll,constits.Q1.phases,'linewidth',2,'DisplayName','Q_1 CWT\_Multi')
% hold(ax1,'off')
plotPhase(constits.decTimesAll,constits.Q1.phases,VanNS.dtime,VanNS.Q1ph,ax1,constits.fband.ci_ph(3),2,'Q_1',-145)
grid on
legend('location','northeast')
title('Vancouver Q_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2=nexttile;
% plot(VanNS.dtime,VanNS.O1ph,'linewidth',2,'DisplayName','O_1 NS\_Tide')
% hold(ax2,'on')
% plot(constits.decTimesAll,constits.O1.phases,'linewidth',2,'DisplayName','O_1 CWT\_Multi')
% hold(ax2,'off')
plotPhase(constits.decTimesAll,constits.O1.phases,VanNS.dtime,VanNS.O1ph,ax2,constits.fband.ci_ph(4),2,'O_1',-145)
grid on
legend('location','northeast')
title('Vancouver O_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ylim([-160 400])
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

ax3=nexttile;
% plot(VanNS.dtime,VanNS.K1ph,'linewidth',2,'DisplayName','K_1 NS\_Tide')
% hold(ax3,'on')
% plot(constits.decTimesAll,constits.K1.phases,'linewidth',2,'DisplayName','K_1 CWT\_Multi')
% hold(ax3,'off')
plotPhase(constits.decTimesAll,constits.K1.phases,VanNS.dtime,VanNS.K1ph,ax3,constits.fband.ci_ph(6),2,'K_1',-145)
grid on
legend('location','northeast')
title('Vancouver K_1', 'Units', 'normalized', 'Position', [0.5, 1.05,0])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylim([-160 400])
ylabel(p,'Phase (deg)','FontSize',18)
plotAnn('(c)')

set(gcf,'Position',[100 100 2500 2500])
% saveas(p,'./figs/Vancouver_NS_CWT_D1_phases.png')

%%
load('./data/Artificial/artifS.mat')

dat.dtime.TimeZone = 'UTC';

% [constits,species,ref,admit] = cwtMulti(dat.dtime,dat.wlnoise);
% 
% [astDIConstits,astDISpecies,astDIRef,astAdmit] = cwtMulti(Ast.dates,Ast.wl,'performAdmittance',lon,lat);

%%
sp_names = ["D_1","D_2","D_3","D_4","D_5","D_6","D_7","D_8","D_9","D_{10}","D_{11}","D_{12}"];


deltaT = seconds(dat.dtime(2) - dat.dtime(1));

c = dat.wlnoise;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

dat.input_pwr=pxx/3600/24;
dat.input_freqs=f*3600*24;      % frequency in cpd

c = constits.reconstruction.reconAll; %

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

dat.recon_pwr=pxx/3600/24;
dat.recon_freqs=f*3600*24;      % frequency in cpd

deltaT = 3600.;

c = astDIRef.ref_wl;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

input.pwr=pxx/3600/24;
input.freqs=f*3600*24;      % frequency in cpd

astDIRef.input_pwr = input.pwr;
astDIRef.input_freqs = input.freqs;

c = Ast.wl;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

input.pwr=pxx/3600/24;
input.freqs=f*3600*24;      % frequency in cpd

Ast.input_pwr = input.pwr;
Ast.input_freqs = input.freqs;



figure(); p=tiledlayout(2,1);

ax1=nexttile;
loglog(astDIRef.input_freqs.^-1,astDIRef.input_pwr,'b-','linewidth',2,'DisplayName','Astoria tidal potential')
hold(ax1,'on')
% plot(dat.input_freqs.^-1,dat.input_pwr,'k-','linewidth',2,'DisplayName','Artificial water level')
plot(Ast.input_freqs.^-1,Ast.input_pwr,'r-','linewidth',1,'DisplayName','Astoria water level')
for i=1:length(species.input.omegas(3:end))
    x = (species.input.omegas(i+2)/2/pi)^-1/3600/24;
    xline(x,'k--','linewidth',2,'alpha',0.5,'HandleVisibility','off')
    text(x, min(ylim)*10^-3,sp_names(i), 'Horiz','left', 'Vert','bottom','FontSize',16) 
end
hold(ax1,'off')
grid on
xlim([0 35])
yma = max(astDIRef.input_pwr)*100;
ymi = min(astDIRef.input_pwr)*10^-3;
ylim([ymi yma])
set(gca,'xdir','reverse')
plotAnnSpec('(a)')
% xlabel('Period [days]')
legend('Location','southwest')
% axis tight
title('Astoria spectra', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;

ylabel(p,'PSD [m^2 cpd^{-1}]','FontSize',18)

ax2=nexttile;
loglog(dat.input_freqs.^-1,dat.input_pwr,'b-','linewidth',2,'DisplayName','Artificial water level')
hold(ax2,'on')
plot(dat.recon_freqs.^-1,dat.recon_pwr,'r-','linewidth',1,'DisplayName','Artificial reconstruction')
for i=1:length(species.input.omegas(3:end))
    x = (species.input.omegas(i+2)/2/pi)^-1/3600/24;
    xline(x,'k--','linewidth',2,'alpha',0.5,'HandleVisibility','off')
    text(x, min(ylim),sp_names(i), 'Horiz','left', 'Vert','bottom','FontSize',16) 
end
hold(ax2,'off')
grid on
xlim([0 35])
ylim auto
yma = max(astDIConstits.fband.pwrs)*1.2;
ymi = min(astDIConstits.fband.pwrs)-10^-8;
set(gca,'xdir','reverse')
plotAnnSpec('(b)')
xlabel('Period [days]')
legend('Location','southwest')
% axis tight
title('Artificial spectra', 'Units', 'normalized', 'Position', [0.5, 1.05, 0])
ax = gca;
ax.FontSize = 14;

ylabel(p,'PSD [m^2 cpd^{-1}]','FontSize',18)

set(gcf,'Position',[100 100 1800 800])
saveas(p,'./figs/Ast_art_TPot_spectra.png')

%% Filter with lobes
nfft = 2^12;
c = constits.M2.filter;
c = c - mean(c);
dt = 360;

% pxx = (fft(c,nfft)/nfft);
% f = linspace(-1/dt/2,1/dt/2 ,nfft);    % - 1/dt/nfft  

[pxx,f] = periodogram(real(c),[],[],1/dt);

pwr=pxx/3600/24;
freqs=f*3600*24;   

semilogy(freqs,pow2db(pwr))
% set(gca,'xdir','reverse')


%%
figure(); p=tiledlayout(1,1);
plot(constits.decTimesAll,matts_unwrap(constits.S2.phases,300))
hold on
plot(constits.decTimesAll,constits.S2.phases)
plot(VanNS.dtime,matts_unwrap(VanNS.S2ph,300))
hold off
% ylim([-1000 1000])

% a = selective_mod(2000,matts_unwrap(constits.N2.phases,300),300);


%% %%%%%%%%%%%%%% PLOTTING FUNCTION(S)
function plotAnn(txt)
    % add text to upper left corner of any plot
    text(min(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end

function plotAnnSpec(txt)
    % add text to upper left corner of any plot
    text(max(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
end

function plotAnnSE(txt,btm)
    % add text to upper left corner of any plot
    % xspot = max(xlim)-0.1*(max(xlim)-min(xlim));
    text(max(xlim)-0.025*(max(xlim)-min(xlim)),btm,txt, 'Horiz','right', 'Vert','bottom','FontSize',20,'BackgroundColor','white') 
end

function plotPhase(x1,y1,x2,y2,ax_in,ci,co_no,const_name,ci_txt_bottom)
    % plotting phases after adjusting mean of y2 to match mean of y1
    jump = 300;
    % a = unwrap(y2,jump) - (nanmedian(unwrap(y2,jump)) - nanmedian(unwrap(y1,jump)));
%     a = unwrap(y2,jump) - (nanmedian(y2) - nanmedian(y1));
%     y2_new = mod(a,360);
   
    y1 = y1 - (nanmedian(y1) - nanmedian(y2));
    
    
%     for i=1:(length(y2_new)-1)
%         diff = 
    
    plot(x2,matts_unwrap(y2,jump),'linewidth',2,'DisplayName',strcat(const_name,' NS\_Tide'))
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
                if max(abs(mean(sig_here(end_jump:start_jump)))-abs(sig_here(end_jump)+jump))>100 && abs(sig_here(start_jump+1)-sig_here(end_jump))<jump
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
        if (sig_here(1)-sig_here(end))>360 && (sig_here(1)>0)
            [~,ii] = min(abs(sig_here));
            sig(i-nw+ii-1:end) = sig(i-nw+ii-1:end)+360;
        end
    end
    
    for i=(nw+1):(length(sig)-nw-1)
        sig_here = sig(i-nw:i+nw);
        if (sig_here(1)-sig_here(end))< -360 && (sig_here(1)>0)
            [~,ii] = min(abs(sig_here-360));
            sig(i-nw+ii-1:end) = sig(i-nw+ii-1:end)-360;
        end
    end
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

    % Alternately, just add the text
    txt = strcat('\pm',string(round(cf,3,'significant')));
    text(max(xlim)-0.025*(max(xlim)-min(xlim)),min(ylim)+0.05*(max(ylim)-min(ylim)),txt, 'Horiz','right', 'Vert','bottom','FontSize',20,'BackgroundColor','white') 
end


