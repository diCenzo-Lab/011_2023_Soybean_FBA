%% code for adding energy requirements for transport reactions


%model=readCbModel('finalNodulatedPlant.mat');
rxnFormula800='-> Root_AMMONIUM[c]';
model=addReaction(model,'Root_R_NH4_tx',rxnFormula800, [], 0, 0, 0, 0);

rxnFormula8='-> Nodule_MNXM9[c]';
model=addReaction(model,'Nodule_R_Pi_tx',rxnFormula8, [], 0, 0, 1000000, 0);

rxnFormula10='Root_AMMONIUM[c] -> Leave_AMMONIUM[c]';
model=addReaction(model,'TRS_AMMONIUM',rxnFormula10, [], 0, 0, 1000000, 0);

%% Removing artifacts of old soybean model
biomass=model.rxns(find(contains(model.rxns,'biomass')));
vac=model.mets(find(contains(model.mets,'[v]')));
vacrxns =findRxnsFromMets(model,vac{1});
for n=2:length(vac)
    vacrxns = vertcat(vacrxns,findRxnsFromMets(model,vac{n}));
end
biomass=union(biomass,vacrxns)
  model = tncore_remove_reactions(model, biomass)
  %model = tncore_remove_reactions(model,{'Nodule_R_SULFATE_vc'})

% optimizeCbModel(model)
% cle
%  biomass=setdiff(biomass,'Nodule_R_Pi_biomass')
% lemon=[];
% for n=1:length(biomass)   
% model = tncore_remove_reactions(model, biomass{n})
% rip=optimizeCbModel(model);
% lemon=[lemon,rip.f];
% end
%    jeff = cell2table(horzcat(biomass, num2cell(transpose(lemon))));

pumped=model.mets(find(contains(model.mets,'Pumped-PROTON')));
pumpedrxns =findRxnsFromMets(model,pumped{1});
for n=2:length(pumped)
    pumpedrxns = vertcat(pumpedrxns,findRxnsFromMets(model,pumped{n}));
end
mnxm = findRxnsFromMets(model, 'Leave_PROTON[c]');
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Leave_PROTON[m]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Leave_PROTON[p]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Root_PROTON[c]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Root_PROTON[m]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Root_PROTON[p]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Nodule_MNXM1[m]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Nodule_MNXM1[p]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Nodule_MNXM1[c]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Bacteroid_MNXM1[c]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Bacteroid_MNXM1[m]'));
mnxm = vertcat(mnxm,findRxnsFromMets(model, 'Bacteroid_MNXM1[p]'));
mnxm = intersect(mnxm, pumpedrxns);
pumpedrxns = setdiff(pumpedrxns,mnxm);
% atp=mnxm(find(contains(mnxm,'ATP')));
% mnxm=setdiff(mnxm,atp);
% atp=vertcat(pumpedrxns,atp);
% m=mnxm(find(contains(mnxm,'_m')));
% mc=mnxm(find(contains(mnxm,'_mc')));
% m=setdiff(m,mc)
% p=mnxm(find(contains(mnxm,'_p')));
% mnxm=setdiff(mnxm,p);
% mnxm=setdiff(mnxm,m);

%%
rxnFormula877='Nodule_ATP[c] + Nodule_MNXM2[c] -> Nodule_ADP[c] + Nodule_MNXM9[c] + Nodule_MNXM1[c]';
model=addReaction(model,'Nodule_R_H+_ATPase_c',rxnFormula877, [], 0, 0, 1000000, 0);

rxnFormula1='Root_ATP[c] + Root_WATER[c] -> Root_ADP[c] + Root_Pi[c] + Root_PROTON[c]';
model=addReaction(model,'Root_R_H+_ATPase_c',rxnFormula1, [], 0, 0, 1000000, 0);

%%
rxnFormula6='Leave_PYRUVATE[m]  <=> Leave_PYRUVATE[c]'
model=addReaction(model,'Leave_R_PYRUVATE_H+_mc',rxnFormula6, [], 0, -1000000, 1000000, 0);

rxnFormula7='Leave_Pi[m] <=> Leave_Pi[c]';
model=addReaction(model,'Leave_R_Pi_H+_mc',rxnFormula7, [], 0, -1000000, 1000000, 0);

rxnFormula808='Nodule_MNXM9[m] <=> Nodule_MNXM9[c]';
model=addReaction(model,'Nodule_R_Pi_H+_mc',rxnFormula808, [], 0, -1000000, 1000000, 0);

rxnFormula9=' -> Nodule_MNXM58[c]';
model=addReaction(model,'Nodule_R_SO4_tx',rxnFormula9, [], 0, 0, 1000000, 0);

rxnFormula101='-> Root_MG+2[c] '
model=addReaction(model,'Root_R_Mg_tx',rxnFormula101, [], 0, 0, 1000000, 0);

rxnFormula11='Root_PYRUVATE[m]  <=> Root_PYRUVATE[c]'
model=addReaction(model,'Root_R_PYRUVATE_H+_mc',rxnFormula11, [], 0, -1000000, 1000000, 0);

rxnFormula12='Root_Pi[m] <=> Root_Pi[c]'
model=addReaction(model,'Root_R_Pi_H+_mc',rxnFormula12, [], 0, -1000000, 1000000, 0);

rxnFormula13=' -> Root_Pi[c]'
model=addReaction(model,'Root_R_Pi_tx',rxnFormula13, [], 0, 0, 1000000, 0);

rxnFormula14=' -> Root_SULFATE[c] '
model=addReaction(model,'Root_R_SO4_tx',rxnFormula14, [], 0, 0, 1000000, 0);

optimizeCbModel(model)


