% Comparing CWT\_Multi results to NS\_Tide results for Columbia River stations
clear

%% load data
addpath('../cwtMultiPackage/CWT_Multi')

load('./fromPascal/NS_Tide_comp/NS_TIDE_tests_LCR.mat')

%% Assigning all NS\_Tide outputs to struct
AstNS.t = t;
AstNS.dtime = datetime(t,'ConvertFrom','datenum');

VanNS.t = t;
VanNS.dtime = datetime(t,'ConvertFrom','datenum');

AstNS.wlOrig = wl_Ast;
AstNS.recon = poutAst;
nonNanInd = find(~isnan(AstNS.recon));
AstNS.rmse = mean(sqrt((AstNS.wlOrig(nonNanInd) - AstNS.recon(nonNanInd)).^2));

VanNS.wlOrig = wl_Van;
VanNS.recon = poutVan;
nonNanInd = find(~isnan(VanNS.recon));
VanNS.rmse = mean(sqrt((VanNS.wlOrig(nonNanInd) - VanNS.recon(nonNanInd)).^2));

AstNS.Q1 = tidecon_tsAst{1,3}{2,1}(3,:);
AstNS.O1 = tidecon_tsAst{1,3}{2,1}(5,:);
AstNS.P1 = tidecon_tsAst{1,3}{2,1}(9,:);
AstNS.K1 = tidecon_tsAst{1,3}{2,1}(10,:);

VanNS.Q1 = tidecon_tsVan{1,3}{2,1}(2,:);
VanNS.O1 = tidecon_tsVan{1,3}{2,1}(3,:);
VanNS.K1 = tidecon_tsVan{1,3}{2,1}(5,:);

AstNS.N2 = tidecon_tsAst{1,3}{3,1}(3,:);
AstNS.M2 = tidecon_tsAst{1,3}{3,1}(4,:);
AstNS.S2 = tidecon_tsAst{1,3}{3,1}(6,:);

VanNS.N2 = tidecon_tsVan{1,3}{3,1}(3,:);
VanNS.M2 = tidecon_tsVan{1,3}{3,1}(4,:);
VanNS.S2 = tidecon_tsVan{1,3}{3,1}(6,:);

VanNS.MK3 = tidecon_tsVan{1,3}{4,1}(2,:);
VanNS.M4 = tidecon_tsVan{1,3}{5,1}(2,:);

AdmNS.Q1 = VanNS.Q1./AstNS.Q1;
AdmNS.O1 = VanNS.O1./AstNS.O1;
AdmNS.K1 = VanNS.K1./AstNS.K1;

AdmNS.N2 = VanNS.N2./AstNS.N2;
AdmNS.M2 = VanNS.M2./AstNS.M2;
AdmNS.S2 = VanNS.S2./AstNS.S2;


%% run CWT\_Multi analysis
load('./data/Astoria/astoria_wl.mat')
load('./data/Vancouver/vancouver_wl.mat')

% align datasets to common start/end
startDate = datetime(2010,04,01);
endDate = datetime(2011,07,01);

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

%% Define dynamic inference parameters


%% run cwt routine on Vancouver data with dynamic inference
[constits,species,ref,admit] = cwtMulti(Ast.dates,Ast.wl,'dynamicInference'); %,'dynamicInferenceCustom',["K1","P1","S1"],[1/23.93447213,1/24.06588766,1/24],[4383,4383,4383]); 


%%
p=figure();
plot(constits.decTimesAll,constits.dynamicInferenceCust.K1.shortAmps,'linewidth',2,'DisplayName','K1 inferred')
hold on
plot(constits.decTimesAll,constits.dynamicInferenceCust.P1.shortAmps,'linewidth',2,'DisplayName','P1 inferred')
plot(constits.decTimesAll,constits.dynamicInferenceCust.S1.shortAmps,'linewidth',2,'DisplayName','S1 inferred')
plot(constits.decTimesAll,constits.dynamicInferenceCust.K1.sixMoAmp,':','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','K1 (6 months filter)')
plot(constits.decTimesAll,constits.dynamicInferenceCust.P1.sixMoAmp,':','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','P1 (6 months filter)')
plot(constits.decTimesAll,constits.dynamicInferenceCust.S1.sixMoAmp,':','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','S1 (6 months filter)')
hold off
grid on
legend('Location','Northwest')
title('Three-way dynamic inference diurnal results at Astoria')
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel('Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 1800 600])
% saveas(p,'./figs/AstDynInfCustDiurn.png')


p=figure();
plot(constits.decTimesAll,constits.K1.amps,'linewidth',2,'DisplayName','K1 (2 weeks filters)')
hold on
plot(constits.decTimesAll,abs(constits.dynamicInferenceCust.P1toGroup.ratio),'linewidth',2,'DisplayName','P1/G')
plot(constits.decTimesAll,abs(constits.dynamicInferenceCust.S1toGroup.ratio),'linewidth',2,'DisplayName','S1/G')
hold off
grid on
legend('Location','Northwest')
title('Further dynamic inference results')
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel('Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 1800 600])
% saveas(p,'./figs/AstDynInfCustDiurn2.png')
