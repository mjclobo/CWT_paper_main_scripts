%% example of overfitting using UTide versus CWT_Multi increased resolution

% loading routines
addpath('../cwtMultiPackage/CWT_Multi')
addpath('./utide_tool/')

% loading data
load('./data/Artificial/artifS_1yr_noise.mat')

% 15-day analysis window
startDate = datetime(1970,01,01);
endDate   = datetime(1970,12,31);

dat.dtime.TimeZone = 'UTC';
startDate.TimeZone='UTC'; endDate.TimeZone='UTC';

startInd = find(dat.dtime==startDate);
endInd = find(dat.dtime==endDate);

% define constituents to use in analysis
co = {'Q1  ','O1  ','K1  ','N2  ','M2  ','S2  '};

dat.wl = dat.wlnoise10;

% analysis
coef = ut_solv(datenum(dat.dtime(startInd:endInd)), dat.wl(startInd:endInd),[],46+12.4/60,co,'RunTimeDisp',...
            'nnn','White','OLS','OrderCnstit','frq','NodsatNone','GwchNone');

coef.mean = 0.0;
coef.slope = 0.0;

% synthesis
ut_recon = ut_reconstr(datenum(dat.dtime),coef);

%% CWT analysis and synthesis
coNames = ["Q1", "O1", "K1", "N2", "M2", "S2"];
coFiltLength = (length(dat.wl(startInd:endInd))-1)*ones(1,6);


spFL = (length(dat.wl(startInd:endInd))-1)*ones(1,14);         % filter lengths in hours for 4d, 2d, d1, d2, d3, d4, d6, d8, ...

[constits,species,~,~] = cwtMulti(dat.dtime(startInd:endInd), dat.wl(startInd:endInd),'spFiltLength',spFL,...
    'coFreqs',{[1/26.868350, 1/25.81933871, 1/23.93447213],[1/12.65834751, 1/12.4206012, 1/12],[],[]},coFiltLength,coNames);
% 
% [~,species,~,~] = cwtMulti(dat.dtime(startInd:endInd), dat.wl(startInd:endInd),'spFiltLength',spFL,...
%     'coFreqs',{[(0.0417807462 * 1.0)],[(0.0805114007 * 1.0)],[],[]},...
%      (length(dat.wl(startInd:endInd))-1)*ones(1,2),["D1", "D2"]);

CWT_coef = coef;
CWT_coef.A(1) = constits.Q1.amps(1); CWT_coef.g(1) = constits.Q1.phases(1);
CWT_coef.A(2) = constits.O1.amps(1); CWT_coef.g(2) = constits.O1.phases(1);
CWT_coef.A(3) = constits.K1.amps(1); CWT_coef.g(3) = constits.K1.phases(1);
CWT_coef.A(4) = constits.N2.amps(1); CWT_coef.g(4) = constits.N2.phases(1);
CWT_coef.A(5) = constits.M2.amps(1); CWT_coef.g(5) = constits.M2.phases(1);
CWT_coef.A(6) = constits.S2.amps(1); CWT_coef.g(6) = constits.S2.phases(1);

CWT_coef_sp = coef;
CWT_coef_sp.name = {'D1SP';'D2SP'};
CWT_coef_sp.A = [species.D1.amps(1); species.D2.amps(1)];
CWT_coef_sp.g = [species.D1.phases(1); species.D2.phases(1)];
CWT_coef_sp.aux.frq = [(0.0417807462 * 0.96882); (0.0805114007 * 1.00764)];

% hacky but works
CWT_coef_sp.A_ci = [0.001,0.001]'; % CWT_coef.A_ci(1:2);
CWT_coef_sp.g_ci = [0.005,0.005]'; % CWT_coef.g_ci(1:2);

CWT_coef_sp.diagn.SNR = [100.,100.]';

% CWT_coef_sp.aux.lind = [coef.aux.lind(3); coef.aux.lind(5)];
CWT_coef_sp.aux.lind = [147; 148];


% do recon with UTide routine
CWT_recon = ut_reconstr(datenum(dat.dtime),CWT_coef,'MinSNR',0.0001);
CWT_recon_sp = ut_reconstr(datenum(dat.dtime),CWT_coef_sp,'MinSNR',0.0001);


%% FFT stuff
% params
deltaT = 3600.;
nfft = 2^12;

% UTide resid spectrum
c = dat.wl - ut_recon;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

UT_resid_pwr=pxx/3600/24;
UT_resid_freqs=f*3600*24;      % frequency in cpd

% CWT resid spectrum
c = dat.wl - CWT_recon;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

CWT_resid_pwr=pxx/3600/24;
CWT_resid_freqs=f*3600*24;      % frequency in cpd

% input spectrum
c = dat.wl;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

input_pwr=pxx/3600/24;
input_freqs=f*3600*24;      % frequency in cpd

% UTide reconstruction spectrum
c = ut_recon;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

UT_recon_pwr=pxx/3600/24;
UT_recon_freqs=f*3600*24;      % frequency in cpd

% CWT reconstruction spectrum
c = CWT_recon;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

CWT_recon_pwr=pxx/3600/24;
CWT_recon_freqs=f*3600*24;      % frequency in cpd

% CWT reconstruction spectrum
c = CWT_recon_sp;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

CWT_recon_pwr_sp=pxx/3600/24;
CWT_recon_freqs_sp=f*3600*24;      % frequency in cpd
%% plot
figure()
[~,end_ind] = min(abs(UT_resid_freqs - 6.2));

semilogy(UT_resid_freqs,UT_resid_pwr,'linewidth',1.5,'DisplayName','UTide residual')
hold('on')
% plot(CWT_resid_freqs,CWT_resid_pwr,'linewidth',1.5,'DisplayName','CWT species residual')
plot(input_freqs,input_pwr,'linewidth',1.5,'DisplayName','input wl')
% plot(UT_recon_freqs,UT_recon_pwr,'linewidth',1.5,'DisplayName','UTide recon wl')
% plot(CWT_recon_freqs,CWT_recon_pwr,'linewidth',1.5,'DisplayName','CWT recon wl')
hold('off')
legend('FontSize',20)
% title('Rmin = 0.3')
xlabel('Frequency [cpd]','FontSize',18)
title("UTide residual power spectrum (3 days)",'FontSize',20.)
ylim([min(vertcat(UT_resid_pwr(1:end_ind),input_pwr(1:end_ind))) max(vertcat(UT_resid_pwr,input_pwr))])
xlim([0.1 6.2])
ylabel('PSD [m^2 cpd^{-1}]','FontSize',18)


%% plot
figure()
[~,end_ind] = min(abs(UT_resid_freqs - 6.2));

semilogy(CWT_resid_freqs,CWT_resid_pwr,'linewidth',1.5,'DisplayName','CWT species residual')
hold('on')
plot(input_freqs,input_pwr,'linewidth',1.5,'DisplayName','input wl')
% plot(UT_recon_freqs,UT_recon_pwr,'linewidth',1.5,'DisplayName','UTide recon wl')
% plot(CWT_recon_freqs,CWT_recon_pwr,'linewidth',1.5,'DisplayName','CWT recon wl')
hold('off')
legend('FontSize',20)
% title('Rmin = 0.3')
xlabel('Frequency [cpd]','FontSize',18)
title("CWT\_Multi residual power spectrum (3 days)",'FontSize',20.)
ylim([min(vertcat(input_pwr(1:end_ind),CWT_resid_pwr(1:end_ind))) max(vertcat(input_pwr,CWT_resid_pwr))])
xlim([0.1 6.2])
ylabel('PSD [m^2 cpd^{-1}]','FontSize',18)


%% more plots
p = figure();

plot(dat.dtime,dat.wl,'DisplayName','input wl')
hold on
plot(dat.dtime,CWT_recon_sp,'DisplayName','CWT\_Multi species reconstruction')
% plot(dat.dtime,dat.wl - CWT_recon_sp,'DisplayName','CWT species residual')
hold off
legend('FontSize',16,'Location','northwest')
ylim([-3.0 6.0])
ylabel('Amplitude [m]','FontSize',18)
title("CWT\_Multi species results (2 day analysis)",'FontSize',20.)

axes('Position',[.65 .7 .2 .2])
box on
plot(dat.dtime(startInd:endInd), dat.wl(startInd:endInd))
hold on
plot(dat.dtime(startInd:endInd),CWT_recon_sp(startInd:endInd))
hold off
box off

set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/overfitting_2d_CWT_sp.png')


%
p = figure();

dec = 20;

plot(dat.dtime,dat.wl,'DisplayName','input wl')
hold on
plot(dat.dtime(1:dec:end),CWT_recon(1:dec:end),'DisplayName','CWT\_Multi contits reconstruction','linewidth',0.05)
% plot(dat.dtime(1:dec:end),dat.wl(1:dec:end) - CWT_recon(1:dec:end),'DisplayName','CWT contits residual','linewidth',0.05)
hold off
legend('FontSize',16,'Location','northwest')
ylim([-35 55.0])
ylabel('Amplitude [m]','FontSize',18)
title("CWT\_Multi constituents results (2 day analysis)",'FontSize',20.)

axes('Position',[.65 .7 .2 .2])
box on
plot(dat.dtime(startInd:endInd), dat.wl(startInd:endInd))
hold on
plot(dat.dtime(startInd:endInd),CWT_recon(startInd:endInd))
hold off
box off

set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/overfitting_2d_CWT_co.png')

%
p = figure();

plot(dat.dtime,dat.wl,'DisplayName','input wl')
hold on
plot(dat.dtime(1:dec:end),ut_recon(1:dec:end),'DisplayName','UTide reconstruction','linewidth',0.05)
% plot(dat.dtime(1:dec:end),dat.wl(1:dec:end) - ut_recon(1:dec:end),'DisplayName','UTide residual','linewidth',0.05)
hold off
legend('FontSize',16,'Location','northwest')
ylim([-35 55.0])
ylabel('Amplitude [m]','FontSize',18)
title("UTide results (2 day analysis)",'FontSize',20.)

axes('Position',[.65 .7 .2 .2])
box on
plot(dat.dtime(startInd:endInd), dat.wl(startInd:endInd))
hold on
plot(dat.dtime(startInd:endInd),ut_recon(startInd:endInd))
hold off
box off

set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/overfitting_2d_UT.png')







