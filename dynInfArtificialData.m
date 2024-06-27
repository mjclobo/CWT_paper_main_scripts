% Comparing CWT\_Multi results to NS\_Tide results for Columbia River stations
clear

%% load data
addpath('../cwtMultiPackage/cwtMulti')

load('./data/Artificial/artifS.mat')

%% run cwt routine on Vancouver data with dynamic inference
[constits,species,ref,admit] = cwtMulti(dat.dtime,dat.wl,'dynamicInferenceCustom',["K1","P1","S1"],[1/23.93447213,1/24.06588766,1/24],[4383,4383,4383]); 


%%
base = ones(size(constits.decTimesAll));
p=figure();
plot(constits.decTimesAll,constits.dynamicInferenceCust.K1.shortAmps,'linewidth',2,'DisplayName','K1 inferred')
hold on
plot(constits.decTimesAll,constits.dynamicInferenceCust.P1.shortAmps,'linewidth',2,'DisplayName','P1 inferred')
plot(constits.decTimesAll,constits.dynamicInferenceCust.S1.shortAmps,'linewidth',2,'DisplayName','S1 inferred')
plot(constits.decTimesAll,constits.dynamicInferenceCust.K1.sixMoAmp,':','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','K1 (6 months filter)')
plot(constits.decTimesAll,constits.dynamicInferenceCust.P1.sixMoAmp,':','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','P1 (6 months filter)')
plot(constits.decTimesAll,constits.dynamicInferenceCust.S1.sixMoAmp,':','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','S1 (6 months filter)')
plot(constits.decTimesAll,base*(53011/90809),'--','color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','K1 true')
plot(constits.decTimesAll,base*(17543/90809),'--','color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','P1 true')
plot(constits.decTimesAll,base*(416/90809),'--','color',[0.9290 0.6940 0.1250],'linewidth',2,'DisplayName','S1 true')
hold off
grid on
legend('Location','Northwest')
title('Three-way dynamic inference diurnal results (Artificial)')
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;
% xlabel(p,'Time','FontSize',18)
ylabel('Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 1800 600])
saveas(p,'./figs/ArtifDynInfCustDiurn.png')

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
saveas(p,'./figs/ArtifDynInfCustDiurn2.png')
