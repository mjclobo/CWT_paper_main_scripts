% script for the response function plot
% very messy

%% run cwt routine on artficial data
addpath('../cwtMultiPackage/CWT_Multi')

load('./data/Artificial/artifS.mat')

dat.dtime.TimeZone = 'UTC';

[constits,species,ref,admit] = cwtMulti(dat.dtime,dat.wlnoise);

%% plotting freq response functions
    
N = length(constits.names);

N_om = 501;                             % must be odd

responses = zeros(N,N_om);
freqs = zeros(N,N_om);

for k=1:N
    [responses(k,:),freqs(k,:)] = cwt_freq_response(constits.(constits.names{k}).filter,constits.omegas(k),abs(constits.omegas(k)/500),N_om,3600.0);
end

figure(); p=tiledlayout(3,1);

ax1 = nexttile;
hold(ax1,'on')
for k=1:N
    plot(((freqs(k,:)*24*3600)/(2*pi)),responses(k,:),'linewidth',2)
end

hold(ax1,'off')
grid on
ylim([0 1])
xlim([0.725 1.2])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(a)')

ax2 = nexttile;
hold(ax2,'on')
for k=1:N
    plot(((freqs(k,:)*24*3600)/(2*pi)),responses(k,:),'linewidth',2)
end

hold(ax2,'off')
ylim([0 1])
xlim([1.75 2.15])
grid on
ax = gca;
ax.FontSize = 14;
plotAnn('(b)')

% add species
N = length(species.names);

N_om = 501;                             % must be odd

responses = zeros(N,N_om);
freqs = zeros(N,N_om);

for k=1:N
    [responses(k,:),freqs(k,:)] = cwt_freq_response(species.(species.names{k}).filter,species.omegas(k),abs(species.omegas(k)/500),N_om,3600.0);
end

ax3=nexttile;
hold(ax3,'on')
for k=1:N
    plot(((freqs(k,:)*24*3600)/(2*pi)),responses(k,:),'linewidth',2)
end

hold(ax3,'off')
grid on
ylim([0 1])
xlim([-1 12])
ax = gca;
ax.FontSize = 14;
plotAnn('(c)')

xlabel(p,'Freq [cy day^{-1}]','FontSize',18)
ylabel(p,'Response','FontSize',18)

% set(gcf,'Position',[100 100 1000 1300])
set(gcf,'Position',[100 100 1400 900])
saveas(p,'./figs/filter_responses.png')
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS
function plotAnn(txt)
    % add text to upper left corner of any plot
    text(min(xlim), max(ylim),txt, 'Horiz','left', 'Vert','top','FontSize',30) 
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
