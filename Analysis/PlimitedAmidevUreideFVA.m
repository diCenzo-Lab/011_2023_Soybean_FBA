clear
close all
%model=readCbModel('NEWcombinedModel.mat')
model=readCbModel('combined_model_with_databases_subs.mat')

changeCobraSolver ('glpk')
%changeCobraSolver ('ibm_cplex');
model = changeRxnBounds(model,'Nodule_Pi_tx',4,'u');
model = changeRxnBounds(model,'Root_Pi_tx',1000,'u');



%% add reversibility rules for every reaction to model
model.rev=[];
for i=1:length(model.rxns)
    if model.lb(i) < 0
        model.rev(i,1)=1;
    else
        model.rev(i,1)=0;
    end
end
%% solving 
sol_all = optimizeCbModel(model);
all_rgr=sol_all.f*24/1000
ratio_cn_all = -1 * (12 * sol_all.x(findRxnIDs(model, 'Nodule_CO2_tx'))) / (14 * sol_all.x(findRxnIDs(model, 'Bacteroid_convFixed')))
[minFlux, maxFlux] = fluxVariability(model, 99, 'max', {'Nodule_CO2_tx'});
ratio_cn_all_min = -1 * (12 * minFlux) / (14 * sol_all.x(findRxnIDs(model, 'Bacteroid_convFixed')));
ratio_cn_all_max = -1 * (12 * maxFlux) / (14 * sol_all.x(findRxnIDs(model, 'Bacteroid_convFixed')));
minFluxAMIDE=abs(ratio_cn_all_min);
maxFluxAMIDE=abs(ratio_cn_all_max);
Nit=find(contains(model.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N_all=sol_all.v(Nit);
Car=find(contains(model.rxns,'Nodule_CO2_tx'));
C_all=sol_all.v(Car)


%% Generating ureide only model 
model_ureide=model;
model_ureide = changeRxnBounds(model_ureide, 'TNR_ASN', 0, 'b');
model_ureide = changeRxnBounds(model_ureide, 'TNR_GLN', 0, 'b');
sol_ureide = optimizeCbModel(model_ureide);
ureide_rgr=sol_ureide.f*24/1000
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
%% Allantoate model
model_allantoate=model_ureide;
model_allantoate = changeRxnBounds(model_allantoate, 'TNR_S-ALLANTOIN', 0, 'b');
sol_allantoate = optimizeCbModel(model_allantoate);
ratio_cn_allantoate = -1 * (12 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Nodule_CO2_tx'))) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')))
[minFlux, maxFlux] = fluxVariability(model_allantoate, 99, 'max', {'Nodule_CO2_tx'});
ratio_cn_allantoate_min = -1 * (12 * minFlux) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')));
ratio_cn_allantoate_max = -1 * (12 * maxFlux) / (14 * sol_allantoate.x(findRxnIDs(model_allantoate, 'Bacteroid_convFixed')));
minFluxALLAN=abs(ratio_cn_allantoate_min);
maxFluxALLAN=abs(ratio_cn_allantoate_max);

amide=transpose([maxFluxAMIDE (minFluxAMIDE-maxFluxAMIDE)]);
udeide=transpose([maxFluxUREIDE (minFluxUREIDE-maxFluxUREIDE)]);
allan=transpose([maxFluxALLAN (minFluxALLAN-maxFluxALLAN)]);
jeff = cell2table(horzcat(num2cell(amide),num2cell(udeide),num2cell(allan)));
  writetable(jeff, 'Bacteroid.txt', 'Delimiter', '\t')
figure(3) 
x=categorical(["Amide","Ureide","Allantoate"]);
x = reordercats(x,{'Amide','Ureide','Allantoate'});
 y=[maxFluxAMIDE (minFluxAMIDE-maxFluxAMIDE); maxFluxUREIDE (minFluxUREIDE-maxFluxUREIDE); maxFluxALLAN (minFluxALLAN-maxFluxALLAN)];
ba=bar(x,y,'stacked','FaceColor','flat','EdgeColor',[1 1 1])
ba(1).CData=[1 1 1];
ba(2).CData=[0 0.4470 0.7410];
%ylim([2 4])
set(gca,'LineWidth',2,'FontSize',50)
  ylabel('Carbon cost of N fixation (g/g/d) ','FontSize',50)
x_width=15.65 ;y_width=14.99;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 %print('amideureideFVAHIGHRGR','-depsc','-loose');
