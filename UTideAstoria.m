% Comparing CWT_Multi to UTide

%% load data
addpath('../cwtMultiPackage/cwtMulti')

load('./data/Astoria/astoria_wl.mat')
load('./data/Vancouver/vancouver_wl.mat')

% align datasets to common start/end
startDate = datetime(2010,04,01);
endDate = datetime(2013,07,01);

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
Ast.lon = 123+46.1/60;            % W
Ast.lat = 46+12.4/60;             % N

Van.lon = 123+46.1/60;            % W
Van.lat = 46+12.4/60;             % N


%% UTide analysis
addpath('/home/matt/Desktop/PSU/S21/CE489/project/')

Ast.datenums = datenum(Ast.dates);
Van.datenums = datenum(Van.dates);

% Define moving parameters etc.
window = 30*12;         % window in days
incr = 7;%1/24;     % increment in days MUST BE SMALL for good reconstruction (i.e. ~<0.1) depends on your dt value below; trying to change this
rayleigh = 0.9;     % specify Rayleigh criterion, implies that we use 'auto' cnstit arg in ut_solv
good_pct = 0.8;     % specify what percentage of data in window must be not NaN

maxFLength = 5*24;   % same cutoff period as low-pass filtering done in cwtMulti

[cd,OUT] = cwt_utide(Van,window,incr,good_pct,rayleigh,maxFLength);%16P1/17K1

