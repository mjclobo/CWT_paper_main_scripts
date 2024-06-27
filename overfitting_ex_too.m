%% last attempt to illustrate overfitting in CWT_Multi versus HA
% loading routines
addpath('../cwtMultiPackage/CWT_Multi')
addpath('./utidec_tool/')

% loading data
load('./data/Artificial/artifS_1yr_noise.mat')

%%
lon = 123+46.1/60;            % W
lat = 46+12.4/60;             % N

%% UTide moving window analysis and reconstruction
% Define moving parameters etc.
window = 30*1.0;         % window in days
incr = 20/24;%1/24;     % increment in days 
rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.8;     % specify what percentage of data in window must be not NaN

maxFLength = 5*24;   % same cutoff period as low-pass filtering done in CWT\_Multi

dat.datenums = datenum(dat.dtime);
dat.dates = dat.dtime;
dat.dates = dat.dates.' ; dat.wl = dat.wl.' ; dat.datenums = dat.datenums.'; dat.dtime = dat.dtime.';
dat.lat = lat; dat.lon = lon;

stat_window=[dat.dtime(1) dat.dtime(end)];

dat.wl0 = dat.wl;

[cd0,OUT0] = cwt_utide(dat,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

dat.wl = dat.wlnoise5.';
[cd5,OUT5] = cwt_utide(dat,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

dat.wl = dat.wlnoise10.';
[cd10,OUT10] = cwt_utide(dat,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1

dat.wl = dat.wlnoise20.';
[cd20,OUT20] = cwt_utide(dat,window,incr,good_pct,rayleigh,maxFLength,stat_window);%16P1/17K1



%% run CWT analysis/dynamic inference on False River data
[constits0,species0,~,~] = cwtMulti(dat.dtime,dat.wl0,'dynamicInference');%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);
[constits5,species5,~,~] = cwtMulti(dat.dtime,dat.wlnoise5,'dynamicInference');%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);
[constits10,species10,~,~] = cwtMulti(dat.dtime,dat.wlnoise10,'dynamicInference');%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);
[constits20,species20,~,~] = cwtMulti(dat.dtime,dat.wlnoise20,'dynamicInference');%,'coFiltLength',coFiltLength);  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);


%% input power spectrum
% input spectrum
c = dat.wl;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

deltaT = 3600.;
nfft = 2^12;

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

input_pwr=pxx/3600/24;
input_freqs=f*3600*24;      % frequency in cpd


%% plotting
[~,end_ind] = min(abs(OUT0.resid_freqs - 6.2));


p=figure();
semilogy(input_freqs,input_pwr,'color','red','linewidth',1.0,'DisplayName','Input wl')
hold on
plot(OUT0.resid_freqs,OUT0.resid_pwr,'linewidth',1.5,'DisplayName','residual (0%)')
plot(OUT5.resid_freqs,OUT5.resid_pwr,'linewidth',1.5,'DisplayName','residual (5%)')
plot(OUT10.resid_freqs,OUT10.resid_pwr,'linewidth',1.5,'DisplayName','residual (10%)')
plot(OUT20.resid_freqs,OUT20.resid_pwr,'linewidth',1.5,'DisplayName','residual (20%)')
hold off
title('UTide residual spectra','Fontsize',20)
xlim([0.25 6.3])
ylim([min(OUT0.resid_pwr(1:end_ind)) max(input_pwr(1:end_ind))])
legend('FontSize',16,'Location','northeast')

xlabel('Frequency [cpd]','FontSize',18)
ylabel('PSD [m^2 cpd^{-1}]','FontSize',18)

set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/overfitting_UTide_noise.png')


%
p=figure();
semilogy(input_freqs,input_pwr,'color','red','linewidth',1.0,'DisplayName','Input wl')
hold on
plot(constits0.fband.freqs,constits0.fband.pwrs,'linewidth',1.5,'DisplayName','residual (0%)')
plot(constits5.fband.freqs,constits5.fband.pwrs,'linewidth',1.5,'DisplayName','residual (5%)')
plot(constits10.fband.freqs,constits10.fband.pwrs,'linewidth',1.5,'DisplayName','residual (10%)')
plot(constits20.fband.freqs,constits20.fband.pwrs,'linewidth',1.5,'DisplayName','residual (20%)')
hold off
title('CWT\_Multi residual spectra','Fontsize',20)
xlim([0.25 6.3])
ylim([min(constits0.fband.pwrs(1:end_ind)) max(input_pwr(1:end_ind))])
legend('FontSize',16,'Location','northeast')

xlabel('Frequency [cpd]','FontSize',18)
ylabel('PSD [m^2 cpd^{-1}]','FontSize',18)

set(gcf,'Position',[100 100 2000 500])
saveas(p,'./figs/overfitting_CWT_Multi_noise.png')






