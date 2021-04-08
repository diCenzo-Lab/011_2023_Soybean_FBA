%% code for adding energy requirements for transport reactions
%changeCobraSolver ('ibm_cplex');
%model=readCbModel('finalNodulatedPlant.mat');
% rxnFormula800='-> Root_MNXM15[c]';
% model=addReaction(model,'Root_NH4_tx',rxnFormula800, [], 0, 0, 1000000, 0); 
% rxnFormula8='-> Nodule_MNXM9[c]';
% model=addReaction(model,'Nodule_Pi_tx',rxnFormula8, [], 0, 0, 1000000, 0); 
% rxnFormula10='Root_MNXM15[c] -> Leave_MNXM15[c]';
% model=addReaction(model,'TRS_AMMONIUM',rxnFormula10, [], 0, 0, 1000000, 0); 

co2=model.rxns(find(contains(model.rxns,'TRN_CO2')));
model = tncore_remove_reactions(model, co2);
%%
tfer = model.rxns(strmatch('TRS_', model.rxns));
tfer = vertcat(tfer,model.rxns(strmatch('TSR_', model.rxns)));
tfer = vertcat(tfer,model.rxns(strmatch('TNR_', model.rxns)));
TRANS = vertcat(tfer,model.rxns(strmatch('TRN_', model.rxns)));
tx=TRANS(find(contains(TRANS,'_tx')));
model = tncore_remove_reactions(model,tx);
TRANS=setdiff(TRANS,tx')
reactions = model.rxns(find(contains(model.rxns,'_pc')));
reactions = vertcat(reactions, model.rxns(find(contains(model.rxns,'_mc')))); 
reactions = vertcat(reactions, model.rxns(find(contains(model.rxns,'_xc'))));
reactions = vertcat(reactions, model.rxns(find(contains(model.rxns,'_tx'))));
reactions = unique(reactions);


reactionsToExclude = {'Leave_R_CO2_mc';'Leave_R_CO2_pc';'Leave_R_CO2_tx';'Leave_R_CO2_xc';'Nodule_R_CO2_pc';'Nodule_R_CO2_tx';...
    'Nodule_R_CO2_xc';'Leave_R_H2O_pc';'Leave_R_H2O_xc';'Leave_R_O2_mc';'Leave_R_O2_pc';'Leave_R_O2_tx';'Leave_R_O2_xc';...
    'Nodule_R_H2O_mc';'Nodule_R_H2O_pc';'Nodule_R_H2O_tx';'Nodule_R_H2O_xc';'Nodule_R_N2_tx';'Nodule_R_O2_mc';'Nodule_R_O2_tx';'Nodule_R_O2_xc';...
    'Root_R_CO2_mc';'Root_R_CO2_pc';'Root_R_CO2_tx';'Root_R_CO2_xc';'Root_R_H2O_tx';'Root_R_H2O_xc';'Root_R_O2_mc';'Root_R_O2_pc';'Root_R_O2_tx';'Root_R_O2_xc';...
    'TRN_WATER';'TRS_WATER';'TSR_WATER';'Leave_R_Photon_tx'};
    reactions = setdiff(reactions, reactionsToExclude);
    TRANS= setdiff(TRANS, reactionsToExclude);

% Identify reactions with one compound
model.S = full(model.S);
reactionsTemp = {};
for n = 1:length(reactions)
    pos = findRxnIDs(model, reactions{n});
    if sum(abs(model.S(:,pos))) < 3
        reactionsTemp = vertcat(reactionsTemp, reactions{n});
    end
end
reactions = reactionsTemp;
% separate TRANS rxns with a proton
TRANSproton = findRxnsFromMets(model, 'Leave_PROTON[c]');
TRANSproton = vertcat(TRANSproton,findRxnsFromMets(model, 'Root_PROTON[c]'));
TRANSproton = vertcat(TRANSproton,findRxnsFromMets(model, 'Nodule_MNXM1[c]'));
TRANSproton = vertcat(TRANSproton,findRxnsFromMets(model, 'Bacteroid_MNXM1[e]'));
TRANSproton = intersect(TRANSproton, TRANS);
TRANS = setdiff(TRANS, TRANSproton);
% separate TRANS rxns with a phosphate
TRANSphosphate = findRxnsFromMets(model, 'Leave_Pi[c]');
TRANSphosphate = vertcat(TRANSphosphate,findRxnsFromMets(model, 'Root_Pi[c]'));
TRANSphosphate = vertcat(TRANSphosphate,findRxnsFromMets(model, 'Nodule_MNXM9[c]'));
TRANSphosphate = vertcat(TRANSphosphate,findRxnsFromMets(model, 'Bacteroid_MNXM9[e]'));
TRANSphosphate = intersect(TRANSphosphate, TRANS);
TRANS = setdiff(TRANS, TRANSphosphate);

% Separate reactions with a proton
protonRxns = findRxnsFromMets(model, 'Leave_PROTON[c]');
protonRxns = vertcat(protonRxns,findRxnsFromMets(model, 'Root_PROTON[c]'));
protonRxns = vertcat(protonRxns,findRxnsFromMets(model, 'Nodule_MNXM1[c]'));
protonRxns = vertcat(protonRxns,findRxnsFromMets(model, 'Bacteroid_MNXM1[e]'));
protonRxns = intersect(protonRxns, reactions);
reactions = setdiff(reactions, protonRxns);

% Separate reactions with a phosphate molecule
phosphateRxns = findRxnsFromMets(model, 'Leave_Pi[c]');
phosphateRxns = vertcat(phosphateRxns,findRxnsFromMets(model, 'Root_Pi[c]'));
phosphateRxns = vertcat(phosphateRxns,findRxnsFromMets(model, 'Nodule_MNXM9[c]'));
phosphateRxns = vertcat(phosphateRxns,findRxnsFromMets(model, 'Bacteroid_MNXM9[e]'));
phosphateRxns = intersect(phosphateRxns, reactions);
reactions = setdiff(reactions, phosphateRxns);

% Separate reactions with a ATP
atpRxns = findRxnsFromMets(model, 'Leave_ATP[c]');
atpRxns = vertcat(atpRxns,findRxnsFromMets(model, 'Root_ATP[c]'));
atpRxns = vertcat(atpRxns,findRxnsFromMets(model, 'Nodule_ATP[c]'));
atpRxns = vertcat(atpRxns,findRxnsFromMets(model, 'Bacteroid_MNXM3[e]'));
atpRxns = intersect(atpRxns, reactions);
reactions = setdiff(reactions, atpRxns);

% Separate reactions with a ADP
adpRxns = findRxnsFromMets(model, 'Leave_ADP[c]');
adpRxns = vertcat(adpRxns,findRxnsFromMets(model, 'Root_ADP[c]'));
adpRxns = vertcat(adpRxns,findRxnsFromMets(model, 'Nodule_ADP[c]'));
adpRxns = vertcat(adpRxns,findRxnsFromMets(model, 'Bacteroid_MNXM7[e]'));
adpRxns = intersect(adpRxns, reactions);
reactions = setdiff(reactions, adpRxns);

% Separate transfer reactions (TSR,TRS,TNR,TRN)

%% Add ATP to reactions with one compound exchanged 
for n = 1:length(reactions)
    if strmatch('Leave_', reactions{n})
ID_atp = findMetIDs(model, 'Leave_ATP[c]');
ID_adp = findMetIDs(model, 'Leave_ADP[c]');
ID_proton = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Leave_Pi[c]');
elseif strmatch('Root_', reactions{n})
ID_atp = findMetIDs(model, 'Root_ATP[c]');
ID_adp = findMetIDs(model, 'Root_ADP[c]');
ID_proton = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Root_Pi[c]');
elseif strmatch('Nodule_', reactions{n})
ID_atp = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate = findMetIDs(model, 'Nodule_MNXM9[c]');
    else
 ID_atp = findMetIDs(model, 'Bacteroid_MNXM3[e]');
ID_adp = findMetIDs(model, 'Bacteroid_MNXM7[e]');
ID_proton = findMetIDs(model, 'Bacteroid_MNXM1[e]');
ID_phosphate = findMetIDs(model, 'Bacteroid_MNXM9[e]');
    end
    
    formula = printRxnFormula(model, reactions{n}, false);
    if contains(formula, '<=>')
        pos = findRxnIDs(model, reactions{n});
        mets = model.mets(logical(abs(model.S(:,pos))));  
        comp1 = reactions{n}(end);
        comp2 = reactions{n}(end-1);
        rxn1 = strcat(model.rxns{pos}(1:end-2),comp1,comp2);
        rxn2 = strcat( model.rxns{pos}(1:end-2),comp2,comp1);
        rxnName1 =rxn1;
        rxnName2 =rxn2;
        model = tncore_remove_reactions(model, reactions{n});
        if length(mets) < 2
              model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}}, [1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{1,1}}, [-1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        else
        model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        end
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
        
    else     
        pos = findRxnIDs(model, reactions{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
       
    end
end

%% Deal with proton containing reactions
for n = 1:length(protonRxns)
       if strmatch('Leave_', protonRxns{n})
ID_atp = findMetIDs(model, 'Leave_ATP[c]');
ID_adp = findMetIDs(model, 'Leave_ADP[c]');
ID_proton = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Leave_Pi[c]');
elseif strmatch('Root_', protonRxns{n})
ID_atp = findMetIDs(model, 'Root_ATP[c]');
ID_adp = findMetIDs(model, 'Root_ADP[c]');
ID_proton = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Root_Pi[c]');
elseif strmatch('Nodule_', protonRxns{n})
ID_atp = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate = findMetIDs(model, 'Nodule_MNXM9[c]');
       else 
 ID_atp = findMetIDs(model, 'Bacteroid_MNXM3[e]');
ID_adp = findMetIDs(model, 'Bacteroid_MNXM7[e]');
ID_proton = findMetIDs(model, 'Bacteroid_MNXM1[e]');
ID_phosphate = findMetIDs(model, 'Bacteroid_MNXM9[e]');
    end
 formula = printRxnFormula(model, protonRxns{n}, false);
    if contains(formula, '<=>')
        pos = findRxnIDs(model, protonRxns{n});
        mets = model.mets(logical(abs(model.S(:,pos))));
        comp1 = protonRxns{n}(end);
        comp2 = protonRxns{n}(end-1);
        rxn1 = strcat(model.rxns{pos}(1:end-2),comp1,comp2);
        rxn2 = strcat( model.rxns{pos}(1:end-2),comp2,comp1);
        rxnName1 = rxn1;
        rxnName2 = rxn2;
        if length(mets) < 2
            
            model = tncore_remove_reactions(model, protonRxns{n});
             model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}}, [1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{1,1}}, [-1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        else
            model = tncore_remove_reactions(model, protonRxns{n});
            model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
                false, 0, [], [], [], [], [], [], false, 0);
            model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
                false, 0, [], [], [], [], [], [], false, 0);
        end
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = -0.75;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 1.25;
    else
        pos = findRxnIDs(model, protonRxns{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 1.25;
    end
end

%% Deal with phosphate containing reactions 
for n = 1:length(phosphateRxns)
      if strmatch('Leave_', phosphateRxns{n})
ID_atp = findMetIDs(model, 'Leave_ATP[c]');
ID_adp = findMetIDs(model, 'Leave_ADP[c]');
ID_proton = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Leave_Pi[c]');
elseif strmatch('Root_', phosphateRxns{n})
ID_atp = findMetIDs(model, 'Root_ATP[c]');
ID_adp = findMetIDs(model, 'Root_ADP[c]');
ID_proton = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Root_Pi[c]');
elseif strmatch('Nodule_', phosphateRxns{n})
ID_atp = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate = findMetIDs(model, 'Nodule_MNXM9[c]');
       else 
 ID_atp = findMetIDs(model, 'Bacteroid_MNXM3[e]');
ID_adp = findMetIDs(model, 'Bacteroid_MNXM7[e]');
ID_proton = findMetIDs(model, 'Bacteroid_MNXM1[e]');
ID_phosphate = findMetIDs(model, 'Bacteroid_MNXM9[e]');
    end
   formula = printRxnFormula(model, phosphateRxns{n}, false);
    if contains(formula, '<=>')
               pos = findRxnIDs(model, phosphateRxns{n});
        mets = model.mets(logical(abs(model.S(:,pos))));
         comp1 = phosphateRxns{n}(end);
        comp2 = phosphateRxns{n}(end-1);
        rxn1 = strcat(model.rxns{pos}(1:end-2),comp1,comp2);
        rxn2 = strcat( model.rxns{pos}(1:end-2),comp2,comp1);
        rxnName1 = rxn1;
        rxnName2 = rxn2;
        if length(mets) < 2
             model = tncore_remove_reactions(model, protonRxns{n});
             model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}}, [1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{1,1}}, [-1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        else
        model = tncore_remove_reactions(model, phosphateRxns{n});
        model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        end
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = -0.75;
        model.S(ID_proton,pos) = 0.25;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 1.25;
        model.S(ID_proton,pos) = 0.25;
    
    else
        pos = findRxnIDs(model, phosphateRxns{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 1.25;
        model.S(ID_proton,pos) = 0.25;
    end
end

%% Deal with ATP containing reactions
for n = 1:length(atpRxns)
      if strmatch('Leave_', atpRxns{n})
ID_atp = findMetIDs(model, 'Leave_ATP[c]');
ID_adp = findMetIDs(model, 'Leave_ADP[c]');
ID_proton = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Leave_Pi[c]');
elseif strmatch('Root_', atpRxns{n})
ID_atp = findMetIDs(model, 'Root_ATP[c]');
ID_adp = findMetIDs(model, 'Root_ADP[c]');
ID_proton = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Root_Pi[c]');
elseif strmatch('Nodule_', atpRxns{n})
ID_atp = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate = findMetIDs(model, 'Nodule_MNXM9[c]');
       else 
 ID_atp = findMetIDs(model, 'Bacteroid_MNXM3[e]');
ID_adp = findMetIDs(model, 'Bacteroid_MNXM7[e]');
ID_proton = findMetIDs(model, 'Bacteroid_MNXM1[e]');
ID_phosphate = findMetIDs(model, 'Bacteroid_MNXM9[e]');
    end
    pos = findRxnIDs(model, atpRxns{n});
    mets = model.mets(logical(abs(model.S(:,pos))));
     comp1 = atpRxns{n}(end);
        comp2 = atpRxns{n}(end-1);
        rxn1 = strcat(model.rxns{pos}(1:end-2),comp1,comp2);
        rxn2 = strcat( model.rxns{pos}(1:end-2),comp2,comp1);
        rxnName1 = rxn1;
        rxnName2 = rxn2;
    model = tncore_remove_reactions(model, atpRxns{n});
    model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    pos = findRxnIDs(model, rxn1);
    model.S(ID_atp,pos) = -1.25;
    model.S(ID_adp,pos) = 0.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
    pos = findRxnIDs(model, rxn2);
    model.S(ID_atp,pos) = 0.75;
    model.S(ID_adp,pos) = 0.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
end

%% Deal with ADP containing reactions
for n = 1:length(adpRxns)
      if strmatch('Leave_', adpRxns{n})
ID_atp = findMetIDs(model, 'Leave_ATP[c]');
ID_adp = findMetIDs(model, 'Leave_ADP[c]');
ID_proton = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Leave_Pi[c]');
elseif strmatch('Root_', adpRxns{n})
ID_atp = findMetIDs(model, 'Root_ATP[c]');
ID_adp = findMetIDs(model, 'Root_ADP[c]');
ID_proton = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate = findMetIDs(model, 'Root_Pi[c]');
elseif strmatch('Nodule_', adpRxns{n})
ID_atp = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate = findMetIDs(model, 'Nodule_MNXM9[c]');
       else 
 ID_atp = findMetIDs(model, 'Bacteroid_MNXM3[e]');
ID_adp = findMetIDs(model, 'Bacteroid_MNXM7[e]');
ID_proton = findMetIDs(model, 'Bacteroid_MNXM1[e]');
ID_phosphate = findMetIDs(model, 'Bacteroid_MNXM9[e]');
    end
    pos = findRxnIDs(model, adpRxns{n});
    mets = model.mets(logical(abs(model.S(:,pos))));
     comp1 = adpRxns{n}(end);
        comp2 = adpRxns{n}(end-1);
        rxn1 = strcat(model.rxns{pos}(1:end-2),comp1,comp2);
        rxn2 = strcat( model.rxns{pos}(1:end-2),comp2,comp1);
        rxnName1 = rxn1;
        rxnName2 = rxn2;
    model = tncore_remove_reactions(model, adpRxns{n});
    model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    pos = findRxnIDs(model, rxn1);
    model.S(ID_atp,pos) = -0.25;
    model.S(ID_adp,pos) = -0.75;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
    pos = findRxnIDs(model, rxn2);
    model.S(ID_atp,pos) = -0.25;
    model.S(ID_adp,pos) = 1.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
end
%% Deal with transfer reactions
% Get IDs of metabolites to add

for n=1:length(TRANS)
    if  strmatch('TNR_', TRANS{n}) 
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANS{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = 0.25;
        model.S(ID_proton_nod,pos) = 0.25;
    elseif strmatch('TRN_', TRANS{n})
        ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANS{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = 0.25;
        model.S(ID_proton_nod,pos) = 0.25;
    else
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_shoot = findMetIDs(model, 'Leave_ATP[c]');
ID_adp_shoot = findMetIDs(model, 'Leave_ADP[c]');
ID_proton_shoot = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate_shoot = findMetIDs(model, 'Leave_Pi[c]'); 
      pos = findRxnIDs(model,TRANS{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_shoot,pos) = -0.25;
        model.S(ID_adp_shoot,pos) = 0.25;
        model.S(ID_phosphate_shoot,pos) = 0.25;
        model.S(ID_proton_shoot,pos) = 0.25;
    end
     
   
end
%% deal with proton transfer reactions
for n=1:length(TRANSproton)
    if  strmatch('TNR_', TRANSproton{n}) 
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANSproton{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = 1.25;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = 0.25;
        model.S(ID_proton_nod,pos) = -0.75;
    elseif strmatch('TRN_', TRANSproton{n})
        ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANSproton{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = -0.75;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = 0.25;
        model.S(ID_proton_nod,pos) = 1.25;
    elseif strmatch('TRS_', TRANSproton{n})
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_shoot = findMetIDs(model, 'Leave_ATP[c]');
ID_adp_shoot = findMetIDs(model, 'Leave_ADP[c]');
ID_proton_shoot = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate_shoot = findMetIDs(model, 'Leave_Pi[c]'); 
      pos = findRxnIDs(model,TRANSproton{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = -0.75;
        model.S(ID_atp_shoot,pos) = -0.25;
        model.S(ID_adp_shoot,pos) = 0.25;
        model.S(ID_phosphate_shoot,pos) = 0.25;
        model.S(ID_proton_shoot,pos) = 1.25;
    else 
        ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_shoot = findMetIDs(model, 'Leave_ATP[c]');
ID_adp_shoot = findMetIDs(model, 'Leave_ADP[c]');
ID_proton_shoot = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate_shoot = findMetIDs(model, 'Leave_Pi[c]'); 
      pos = findRxnIDs(model,TRANSproton{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 0.25;
        model.S(ID_proton_root,pos) = 1.25;
        model.S(ID_atp_shoot,pos) = -0.25;
        model.S(ID_adp_shoot,pos) = 0.25;
        model.S(ID_phosphate_shoot,pos) = 0.25;
        model.S(ID_proton_shoot,pos) = -0.75;
    end
     
   
end
%% deal with proton transfer reactions
for n=1:length(TRANSphosphate)
    if  strmatch('TNR_', TRANSphosphate{n}) 
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANSphosphate{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 1.25;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = -0.75;
        model.S(ID_proton_nod,pos) = 0.25;
    elseif strmatch('TRN_', TRANSphosphate{n})
        ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_nod = findMetIDs(model, 'Nodule_ATP[c]');
ID_adp_nod = findMetIDs(model, 'Nodule_ADP[c]');
ID_proton_nod = findMetIDs(model, 'Nodule_MNXM1[c]');
ID_phosphate_nod = findMetIDs(model, 'Nodule_MNXM9[c]');
pos = findRxnIDs(model,TRANSphosphate{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = -0.75;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_nod,pos) = -0.25;
        model.S(ID_adp_nod,pos) = 0.25;
        model.S(ID_phosphate_nod,pos) = 1.25;
        model.S(ID_proton_nod,pos) = 0.25;
    elseif strmatch('TRS_', TRANSphosphate{n})
ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_shoot = findMetIDs(model, 'Leave_ATP[c]');
ID_adp_shoot = findMetIDs(model, 'Leave_ADP[c]');
ID_proton_shoot = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate_shoot = findMetIDs(model, 'Leave_Pi[c]'); 
      pos = findRxnIDs(model,TRANSphosphate{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = -0.75;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_shoot,pos) = -0.25;
        model.S(ID_adp_shoot,pos) = 0.25;
        model.S(ID_phosphate_shoot,pos) = 1.25;
        model.S(ID_proton_shoot,pos) = 0.25;
    else 
        ID_atp_root = findMetIDs(model, 'Root_ATP[c]');
ID_adp_root = findMetIDs(model, 'Root_ADP[c]');
ID_proton_root = findMetIDs(model, 'Root_PROTON[c]');
ID_phosphate_root = findMetIDs(model, 'Root_Pi[c]');
ID_atp_shoot = findMetIDs(model, 'Leave_ATP[c]');
ID_adp_shoot = findMetIDs(model, 'Leave_ADP[c]');
ID_proton_shoot = findMetIDs(model, 'Leave_PROTON[c]');
ID_phosphate_shoot = findMetIDs(model, 'Leave_Pi[c]'); 
      pos = findRxnIDs(model,TRANSphosphate{n});
        model.S(ID_atp_root,pos) = -0.25;
        model.S(ID_adp_root,pos) = 0.25;
        model.S(ID_phosphate_root,pos) = 1.25;
        model.S(ID_proton_root,pos) = 0.25;
        model.S(ID_atp_shoot,pos) = -0.25;
        model.S(ID_adp_shoot,pos) = 0.25;
        model.S(ID_phosphate_shoot,pos) = -0.75;
        model.S(ID_proton_shoot,pos) = 0.25;
    end
     
   
end

        nodulatedPlant=model; 
