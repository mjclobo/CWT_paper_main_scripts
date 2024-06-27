% Comparing CWT_Multi to UTide
clear

%% load data
addpath('../cwtMultiPackage/cwtMulti')

load('./data/Artificial/artifNS.mat')

% load('./data/Lobo_Data.mat')
% 
% 
% dat4.dates  = datetime(d.tn{4},'ConvertFrom','datenum');
% dat4.wl     = d.wl{4};
% 
% dat1.dates  = datetime(d.tn{1},'ConvertFrom','datenum');
% dat1.wl     = d.wl{1};
% 
% %% align datasets to common start/end
% startDate = datetime(2013,06,01);
% endDate = datetime(2018,06,01);
% 
% dat4.start = find(dat4.dates==startDate);
% dat4.end = find(dat4.dates==endDate);
% 
% dat4.dates = dat4.dates(dat4.start:dat4.end);
% dat4.wl = dat4.wl(dat4.start:dat4.end);
% 
% dat1.start = find(dat1.dates==startDate);
% dat1.end = find(dat1.dates==endDate);
% 
% dat1.dates = dat1.dates(dat1.start:dat1.end);
% dat1.wl = dat1.wl(dat1.start:dat1.end);
% 
% nan_ind = unique(sort([find(isnan(dat4.wl));find(isnan(dat1.wl))]));
% 
% dat4.wl(nan_ind) = NaN;
% dat1.wl(nan_ind) = NaN;
% 
% % lon, lat for dat1oria
% dat4.lon = d.info{4}.lon;            % W
% dat4.lat = d.info{4}.lat;             % N

% lon, lat for Astoria
dat.lon = 123+46.1/60;            % W
dat.lat = 46+12.4/60;             % N

%% filtering data so that we look at tides only
% defining filter criteria


% dt = seconds(nanmedian(diff(dat1.dates)));
% 
% stopHrs = 5*24;
% 
% co_stop = (stopHrs*3600)^-1;
% co_pass = ((stopHrs+48)*3600)^-1;     % i.e., pass band of 2 cyc/day
% 
% dpass_pct_lo = 0.1; % 5 % passband ripple
% dstop_db_lo = 0.005; % 47db stopband attenuation
% 
% % defining optimal filters
% [co_n,wn,beta,ftype] = kaiserord([co_pass,co_stop],[1,0],[dpass_pct_lo,dstop_db_lo],dt^-1);
% co_filt_lo = fir1(co_n,wn,ftype,kaiser(co_n+1,beta),'noscale');
% 
% % applying filters to dat1
% dataLo = conv(dat1.wl,co_filt_lo,'same');
% dataHi = dat1.wl - dataLo;
% 
% dat1.wl = dataHi;
% 
% % applying filters to dat4
% dataLo = conv(dat4.wl,co_filt_lo,'same');
% dataHi = dat4.wl - dataLo;
% 
% dat4.wl = dataHi;

%% UTide analysis
% addpath('/home/matt/Desktop/PSU/S21/CE489/project/')
addpath('./utide_tool')

dat.datenums = datenum(dat.dtime);
dat.dates = dat.dtime;
% dat4.datenums = datenum(dat4.dates);

% Define moving parameters etc.
window = 30*1.5;         % window in days
incr = 7;%1/24;     % increment in days MUST BE SMALL for good reconstruction (i.e. ~<0.1) depends on your dt value below; trying to change this
plot_flag = 0;      % do basic plotting (TBD)
rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.9;     % specify what percentage of data in window must be not NaN


maxFLength = 5*24;   % should match max constit filter length being used below (in hours)

[cd,OUT] = cwt_utide(dat,window,incr,plot_flag,good_pct,rayleigh,maxFLength);%16P1/17K1

%% run cwt routine on dat4couver/dat1oria data
coFiltLength    = [1081,1080,1080,1081,1080,1081,1081,1081,1081,1080,425,1081,1080,1081,1080,1080,1080,1080,1080,1080,1080,1080,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours

disp('Beginning cwtMulti analysis')

[constits,species,ref,admit] = cwtMulti(dat.dates,dat.wl,'pfResp');%,'coFiltLength',coFiltLength); % ,'dynamicInference');  %,'performAdmittance',dat4.lon,dat4.lat,'refStation',dat1.wl);

%% plot results (M2 S2)
p=figure();
plot(cd.datetimes,cd.M2.amp,'b-','linewidth',2,'DisplayName','M_2 UTide')
hold on
plot(constits.decTimesAll,constits.M2.amps,'k-','linewidth',2,'DisplayName','M_2 cwtMulti')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off
legend('location','northeast')
title('M_2 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

p=figure();
plot(cd.datetimes,cd.S2.amp,'b-','linewidth',2,'DisplayName','S_2 UTide')
hold on
plot(constits.decTimesAll,constits.S2.amps,'k-','linewidth',2,'DisplayName','S_2 cwtMulti')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
hold off
legend('location','northeast')
title('S_2 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])


%% plot results (N2 L2)
p=figure();
plot(cd.datetimes,cd.N2.amp,'b-','linewidth',2,'DisplayName','N_2 UTide')
hold on
plot(constits.decTimesAll,constits.N2.amps,'k-','linewidth',2,'DisplayName','N_2 cwtMulti')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off
legend('location','northeast')
title('N_2 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

p=figure();
plot(cd.datetimes,cd.L2.amp,'b-','linewidth',2,'DisplayName','L_2 UTide')
hold on
plot(constits.decTimesAll,constits.L2.amps,'k-','linewidth',2,'DisplayName','L_2 cwtMulti')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
hold off
legend('location','northeast')
title('L_2 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])


%% plot results (O1 K1)
p=figure();
plot(cd.datetimes,cd.O1.amp,'b-','linewidth',2,'DisplayName','O_1 UTide')
hold on
plot(constits.decTimesAll,constits.O1.amps,'k-','linewidth',2,'DisplayName','O_1 cwtMulti')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off
legend('location','northeast')
title('O_1 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

%% K1 dyn inf
p=figure();
plot(cd.datetimes,cd.K1.amp,'b-','linewidth',2,'DisplayName','K_1 UTide')
hold on
plot(constits.decTimesAll,constits.K1.amps,'k-','linewidth',2,'DisplayName','K_1 cwtMulti (original 2 weeks reponse)')
plot(constits.decTimesAll,constits.dynamicInference.K1.sixMoAmp(1:end-2),'r-','linewidth',2,'DisplayName','K_1 cwtMulti (6 month filter)')
plot(constits.decTimesAll,constits.dynamicInference.K1.shortAmps,'c-','linewidth',2,'DisplayName','K_1 cwtMulti (dynamically inferred 2 weeks filter)')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
hold off
title('K_1 amplitude comparison')
legend('location','northeast')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

%% P1 dyn inf
% plot(cd.datetimes,cd.P1.amp,'b-','linewidth',2,'DisplayName','P_1 UTide')
% plot(constits.decTimesAll,constits.P1.amps,'k-','linewidth',2,'DisplayName','P_1 cwtMulti (original 2 weeks reponse)')
plot(constits.decTimesAll,constits.dynamicInference.P1.sixMoAmp(1:end-2),'r-','linewidth',2,'DisplayName','P_1 cwtMulti (6 month filter)')
hold on
plot(constits.decTimesAll,constits.dynamicInference.P1.shortAmps,'c-','linewidth',2,'DisplayName','P_1 cwtMulti (dynamically inferred 2 weeks filter)')
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
hold off
title('P_1 amplitude comparison')
legend('location','northeast')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

%% plot results (MK3)
p=figure();
plot(constits.decTimesAll,constits.MK3.amps,'b-','linewidth',2,'DisplayName','MK3 cwtMulti')
hold on
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off
legend('location','northeast')
title('MK3 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])


%% plot results (M4)
p=figure();
plot(constits.decTimesAll,constits.M4.amps,'b-','linewidth',2,'DisplayName','M4 cwtMulti')
hold on
% xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off
legend('location','northeast')
title('M4 amplitude comparison')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])

%% dynamic inference results
% p=figure();
% plot(cd.datetimes,cd.K1.amp,'b-','linewidth',2,'DisplayName','K_1 cwtMulti')
% hold on
% plot(constits.decTimesAll,constits.K1.amps,'k-','linewidth',2,'DisplayName','K_1 dynamic inference')
% % xline(datetime(2015,05,08),'r-','linewidth',2,'DisplayName','construction start')
% hold off
% title('K_1 amplitude comparison')
% legend('location','northeast')
% grid on
% % xlim([datetime(2014,01,01) datetime(2017,12,01)])

%%
p=figure();
plot(OUT.dtimes,OUT.wl,'DisplayName','UTide reconstruction')
hold on
plot(constits.alltimes,constits.reconHi,'DisplayName','cwtMulti reconstruction')
plot(OUT.dtimes,OUT.wl,'DisplayName','UTide reconstruction')

hold off
legend('location','northeast')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])


p=figure();
plot(OUT.dtimes,OUT.wl-constits.reconHi.')
grid on
% xlim([datetime(2014,01,01) datetime(2017,12,01)])
title('UTide reconstruction - cwtMulti reconstruction')

%% plot admittance
% x = constits.decTimesAll;
% plot(x,admit.constits.M2.amps,'linewidth',2)
% hold on
% title('M_2 admittance (False River / San Francisco)')
% hold off
% grid on
% % xlim([datetime(2014,01,01) datetime(2017,12,01)])




