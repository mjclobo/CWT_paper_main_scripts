% Comparing CWT_Multi to UTide
clear

%% load data
addpath('../cwtMultiPackage/cwtMulti')

load('./data/Lobo_Data.mat')

%% align datasets to common start/end
% startDate = datetime(2011,01,01);
% endDate = datetime(2011,07,01);

startDate = datetime(2006,01,01);
endDate = datetime(2014,12,30);

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

Ast.lat=lat;
Van.lat=lat;

%% filtering data so that we look at tides only
% defining filter criteria
dt = seconds(nanmedian(diff(Ast.dates)));

stopHrs = 3*24;

co_stop = (stopHrs*3600)^-1;
co_pass = ((stopHrs+24)*3600)^-1;     % i.e., pass band of 1 cyc/day

dpass_pct_lo = 0.1; % 5 % passband ripple
dstop_db_lo = 0.005; % 47db stopband attenuation

% defining optimal filters
[co_n,wn,beta,ftype] = kaiserord([co_pass,co_stop],[1,0],[dpass_pct_lo,dstop_db_lo],dt^-1);
co_filt_lo = fir1(co_n,wn,ftype,kaiser(co_n+1,beta),'noscale');

% applying filters to Ast
dataLo = conv(Ast.wl,co_filt_lo,'same');
dataHi = Ast.wl - dataLo;

Ast.wl = dataHi;

% applying filters to Van
dataLo = conv(Van.wl,co_filt_lo,'same');
dataHi = Van.wl - dataLo;

Van.wl = dataHi;

%% UTide analysis
addpath('/home/matt/Desktop/PSU/S21/CE489/project/')

Ast.datenums = datenum(Ast.dates);
Van.datenums = datenum(Van.dates);

% Define moving parameters etc.
window = 30*6.3;         % window in days
incr = 31;%1/24;        % increment in days MUST BE SMALL (i.e. ~<0.1) depends on your dt value below
plot_flag = 0;      % do basic plotting (TBD)
rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.9;     % specify what percentage of data in window must be not NaN


maxFLength = 732;   % should match max constit filter length being used below (in hours)

[cd,OUT] = cwt_utide(Ast,window,incr,plot_flag,good_pct,rayleigh,maxFLength);%16P1/17K1

x=r
%% add custom filter lengths for the constituent analysis

coFiltLength = [1080,425,425,1080,425,1080,1080,1080,424,425,425,1080,425,1080,425,425,425,425,425,425,425,425];% 355 semidiurnals, 328 diurnals
coNames = ["TwoQ1","Q1","O1","NO1","K1","J1","OO1","Eps2","Mu2","N2","M2","L2","S2","Eta2","MO3","M3","MK3","SK3","MN4","M4","MS4","S4"];


%% run cwt routine on Vancouver/Astoria data
[constits,species,ref,admit] = cwtMulti(Ast.dates,Ast.wl,lon,lat,'refStation',Van.wl,'coFiltLength',coFiltLength,'dynamicInference',...
    'pfAmps','pfResid','pfResidSpectra');


%% plot results
p=figure();
plot(cd.datetimes,cd.M2.amp,'b-','linewidth',2)
hold on
plot(constits.dectimes,constits.M2.amps,'k-','linewidth',2)
% plot(species.dectimes,species.d2.amps,'r-','linewidth',2)
hold off

p=figure();
plot(cd.datetimes,cd.S2.amp,'b-','linewidth',2)
hold on
plot(constits.dectimes,constits.S2.amps,'k-','linewidth',2)
hold off

%%
p=figure();
plot(OUT.dtimes,OUT.wl)
hold on
plot(constits.alltimes,constits.recon)
hold off

p=figure();
plot(OUT.dtimes,OUT.wl-constits.recon.')
%% example of using more semidiurnal and diurnal frequencies
% constits_omega = cellfun(@(x) ((2*pi)/3600) * x, { ...
%      [0.0372185026, ...                     % Q1
%       0.0387306544, ...                     % O1 
%       0.0417807462, ...                     % K1
%       1/23.09848146], ...                   % added diurnal frequency: J1
%      [0.0789992488, ...                     % N2 
%       0.0805114007, ...                     % M2
%       0.0833333333, ...                     % S2
%       1/12.19162085], ...                   % added semidiurnal frequency: L2
%      [(0.0387306544 + 0.0805114007), ...    % MO3
%        (0.0417807462 + 0.0805114007)], ...  % MK3
%      [(0.0789992488 + 0.0805114007), ...    % MN4   
%        (0.0805114007*2), ...                % M4
%        (0.0805114007 + 0.0833333333)]},...  % MS4
%        'UniformOutput',false);
% 
% constits_names = ["Q1","O1","K1","J1","N2","M2","S2","L2","MO3","MK3","MN4","M4","MS4"]; 
% 
% constits_flength = 363*ones(1,13);                       % filter lengths in hours  for 4d, 2d, d1, d2...
% 
% % run cwt routine again
% [constits,species,poten,admit] = cwt_multi(Van.dates,Van.wl,lon,lat,'ref_station',Ast.wl,'constits_flength',constits_flength,...
%     'constits_freqs',constits_omega,constits_flength,constits_names,'plot_flag_lod2','plot_flag_amps','plot_flag_resid','plot_flag_resid_spectra');
% 



