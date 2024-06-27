% outputting UTide SNR avg. for all stations

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

rayleigh = 0.9;

% %% a single call to UTide for whole time period
% sf_coef_all = ut_solv(dat1.datenums, dat1.wl,[],dat1.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
% 
% fr_coef_all = ut_solv(dat4.datenums, dat4.wl,[],dat4.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
%         
% % export to csv for Davido
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

[cdFR,outFR] = cwt_utide(dat4,window,incr,good_pct,rayleigh,maxFLength,[datetime(2015,01,01) end_date]);%16P1/17K1

[cdSF,outSF] = cwt_utide(dat1,window,incr,good_pct,rayleigh,maxFLength,[datetime(2015,01,01) end_date]);%16P1/17K1

%% CR stuff
%% Align dates in dataset
load('./data/Astoria/astoria_wl.mat')
load('./data/Vancouver/vancouver_wl.mat')

% align datasets to common start/end
startDate = datetime(2010,04,05);
endDate = datetime(2013,06,28);

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
% ref_date = datetime(0,0,0,0,0,0) + (Van.dates(vInd)-datetime(1899,12,31,12,0,0));
ref_date = datetime(1899,12,31,12,0,0);
ref_date.TimeZone = 'UTC';

%% Perform UTide analysis on Astoria/Vancouver
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

[cdVan,outVan] = cwt_utide(Van,window,incr,good_pct,rayleigh,maxFLength,[datetime(2010,10,01) datetime(2013,01,01)]);%16P1/17K1

[cdAst,outAst] = cwt_utide(Ast,window,incr,good_pct,rayleigh,maxFLength,[datetime(2010,10,01) datetime(2013,01,01)]);%16P1/17K1

%% Plot M2 SNR at FR
plot(cdFR.datetimes,cdFR.M2.snr,'linewidth',2,'DisplayName','FR')
hold on
plot(cdSF.datetimes,cdSF.M2.snr,'linewidth',2,'DisplayName','SF')
hold off
legend()
title('M2')

figure();
plot(cdFR.datetimes,cdFR.Q1.snr,'linewidth',2,'DisplayName','FR')
hold on
plot(cdSF.datetimes,cdSF.Q1.snr,'linewidth',2,'DisplayName','SF')
hold off
legend()
title('Q1')

%% output to a table
% Not sure that there exists a function to do this well

tt = table(zeros(4,1),'RowNames',[{'False River'},{'SF Bay'},{'Astoria'},{'Vancouver'}]);
n=1;
for k=1:length(outFR.names)
    nombre = invalidNameIn(outFR.names{k});
    tt = addvars(tt,[cdFR.(nombre).avgSNR...
    ,cdSF.(nombre).avgSNR...
    ,cdAst.(nombre).avgSNR...
    ,cdVan.(nombre).avgSNR]','NewVariableNames',string(outFR.names{k}));%...
    %[string(outFR.names{k}),string(outSF.names{k}),string(outAst.names{k}),string(outVan.names{k})]');
    n = n+1;
end

tt = removevars(tt,'Var1');

writetable(tt,"UTide_avg_SNR_all.csv",'WriteRowNames',true)

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
