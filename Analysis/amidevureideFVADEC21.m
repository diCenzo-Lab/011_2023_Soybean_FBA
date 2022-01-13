function amidevureideFVADEC21
% Analysis without constrained amide model 
clear
close all
%changeCobraSolver ('glpk'); % solver 
model_ureide=readCbModel('OGnewestmodelSEPT21.mat');

%% add reversibility rules for every reaction to model
model_ureide.rev=[];
for i=1:length(model_ureide.rxns)
    if model_ureide.lb(i) < 0
        model_ureide.rev(i,1)=1;
    else
        model_ureide.rev(i,1)=0;
    end
end

%% Solving ureide model
sol_ureide = optimizeCbModel(model_ureide);
ratio_cn_ureide = -1 * (12 * sol_ureide.x(findRxnIDs(model_ureide, 'Nodule_CO2_tx'))) / (14 * sol_ureide.x(findRxnIDs(model_ureide, 'Bacteroid_convFixed')))
[minFlux, maxFlux] = fluxVariability(model_ureide, 99, 'max', {'Nodule_CO2_tx'});
ratio_cn_ureide_min = -1 * (12 * minFlux) / (14 * sol_ureide.x(findRxnIDs(model_ureide, 'Bacteroid_convFixed')));
ratio_cn_ureide_max = -1 * (12 * maxFlux) / (14 * sol_ureide.x(findRxnIDs(model_ureide, 'Bacteroid_convFixed')));
minFluxUREIDE=abs(ratio_cn_ureide_min);
maxFluxUREIDE=abs(ratio_cn_ureide_max);
Nit=find(contains(model_ureide.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N_ureide=sol_ureide.v(Nit);
Car=find(contains(model_ureide.rxns,'Nodule_CO2_tx'));
C_ureide=sol_ureide.v(Car)
%% Generating amide model 
model_amide=model_ureide;
model_amide = addReaction(model_amide, 'TNR_GLN', ...
    'reactionFormula', '0.25 Nodule_ATP[c] + Nodule_MNXM37[c] + 0.25 Root_ATP[c]  -> 0.25 Nodule_ADP[c] + 0.25 Nodule_MNXM1[c] + 0.25 Nodule_MNXM9[c] + 0.25 Root_ADP[c] + Root_GLN[c] + 0.25 Root_PROTON[c] + 0.25 Root_Pi[c]', ...
    'reversible', false, 'upperBound', 1000000);
model_amide = addReaction(model_amide, 'TNR_ASN', ...
    'reactionFormula', '0.25 Nodule_ATP[c] + Nodule_MNXM147[c] + 0.25 Root_ATP[c]  -> 0.25 Nodule_ADP[c] + 0.25 Nodule_MNXM1[c] + 0.25 Nodule_MNXM9[c] + 0.25 Root_ADP[c] + Root_ASN[c] + 0.25 Root_PROTON[c] + 0.25 Root_Pi[c]', ...
    'reversible', false, 'upperBound', 1000000);
model_amide = changeRxnBounds(model_amide, 'TRN_ASN', 0, 'b');
model_amide = changeRxnBounds(model_amide, 'TRN_GLT', 0, 'b');
model_amide = changeRxnBounds(model_amide, 'TRN_GLN', 0, 'b');
model_amide = changeRxnBounds(model_amide, 'EXCT_rev_MNXM131[e]', 0, 'b');
sol_amide = optimizeCbModel(model_amide);
ratio_cn_amide = -1 * (12 * sol_amide.x(findRxnIDs(model_amide, 'Nodule_CO2_tx'))) / (14 * sol_amide.x(findRxnIDs(model_amide, 'Bacteroid_convFixed')))
[minFlux, maxFlux] = fluxVariability(model_amide, 99, 'max', {'Nodule_CO2_tx'});
ratio_cn_amide_min = -1 * (12 * minFlux) / (14 * sol_amide.x(findRxnIDs(model_amide, 'Bacteroid_convFixed')));
ratio_cn_amide_max = -1 * (12 * maxFlux) / (14 * sol_amide.x(findRxnIDs(model_amide, 'Bacteroid_convFixed')));
minFluxAMIDE=abs(ratio_cn_amide_min);
maxFluxAMIDE=abs(ratio_cn_amide_max);
Nit=find(contains(model_amide.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N_amide=sol_amide.v(Nit);
Car=find(contains(model_amide.rxns,'Nodule_CO2_tx'));
C_amide=sol_amide.v(Car)
%% Generating allantoate model
model_allantoate=model_ureide;
model_allantoate = changeRxnBounds(model_allantoate, 'TNR_S-ALLANTOIN', 0, 'b');
sol_allantoate = optimizeCbModel(model_allantoate);
ratio_cn_allantoate = -1 * (12 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Nodule_CO2_tx'))) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')))
[minFlux, maxFlux] = fluxVariability(model_allantoate, 99, 'max', {'Nodule_CO2_tx'});
ratio_cn_allantoate_min = -1 * (12 * minFlux) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')));
ratio_cn_allantoate_max = -1 * (12 * maxFlux) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')));
minFluxALLAN=abs(ratio_cn_allantoate_min);
maxFluxALLAN=abs(ratio_cn_allantoate_max);



% 
% %%
% figure(3) 
% % x=[maxFluxAMIDE maxFluxUREIDE];
% % b=bar(x)
% % b.FaceColor = [0 0.4470 0.7410];
% hold on
% y=[minFluxAMIDE minFluxUREIDE];
% k=bar(y)
% hold on
% set(gca,'LineWidth',2,'FontSize',50)
%  xlabel('Nodule reactions','FontSize',50)
%   ylabel('Flux (\mu molg^{-1}DWh^{-1}) ','FontSize',50)
%    x_width=19.65 ;y_width=14.99;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%  %print('FVAnodule','-depsc','-loose');
%  
 
figure(3) 
x=categorical(["Amide","Ureide","Allantoate"]);
x = reordercats(x,{'Amide','Ureide','Allantoate'});
 y=[maxFluxAMIDE (minFluxAMIDE-maxFluxAMIDE); maxFluxUREIDE (minFluxUREIDE-maxFluxUREIDE); maxFluxALLAN (minFluxALLAN-maxFluxALLAN)];
ba=bar(x,y,'stacked','FaceColor','flat','EdgeColor',[1 1 1])
ba(1).CData=[1 1 1];
ba(2).CData=[0 0.4470 0.7410];
ylim([2 3.5])
set(gca,'LineWidth',2,'FontSize',50)
  ylabel('Carbon cost of N fixation (g/g/d) ','FontSize',50)
x_width=15.65 ;y_width=14.99;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% print('amideureideFVAagainDec21','-depsc','-loose');

end