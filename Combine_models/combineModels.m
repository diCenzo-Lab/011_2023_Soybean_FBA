%% Load the models

% Load ureides model
load('finalNodulatedPlant.mat');
model_ureides = finalNodulatedPlant;

% Load amides model
load('finalNodulatedPlant_amide.mat');
model_amides = finalNodulatedPlant;

clearvars finalNodulatedPlant;

%% Modify the amide model for integration

% Find reactions specific to the amide model
rxn_unique = setdiff(model_amides.rxns, model_ureides.rxns);

% Get optimized growth for amide model
sol_amide = optimizeCbModel(model_amides);

% Test which reactions are required
rxns_to_remove = {};
model_tmp = model_amides;
for n = 1:length(rxn_unique)
    model = tncore_remove_reactions(model_tmp, rxn_unique{n});
    sol_test = optimizeCbModel(model);
    if sol_test.f > 0.99 * sol_amide.f
        rxns_to_remove = vertcat(rxns_to_remove, rxn_unique{n});
        model_tmp = model;
    end
end
rxns_to_remove = setdiff(rxns_to_remove, {'TNR_GLN'; 'TNR_ASN'});

% Remove these reactions from the amide model
model_amides = tncore_remove_reactions(model_amides, rxns_to_remove);

% Check growth
optimizeCbModel(model_amides)

%% Combine the models

% Add rev field to the ureides model, if it does not exist
if isfield(model_ureides, 'rev') == 0
    model_ureides.rev = [];
    for n = 1:length(model_ureides.rxns)
        if model.lb(n,1) < 0
            model_ureides.rev(n,1) = 1;
        else
            model_ureides.rev(n,1) = 0;
        end
    end
end

% Add rev field to the amides model, if it does not exist
if isfield(model_amides, 'rev') == 0
    model_amides.rev = [];
    for n = 1:length(model_amides.rxns)
        if model.lb(n,1) < 0
            model_amides.rev(n,1) = 1;
        else
            model_amides.rev(n,1) = 0;
        end
    end
end

% Find reactions specific to the updated amide model
rxn_unique = setdiff(model_amides.rxns, model_ureides.rxns);

% Remove common reactions from amide model
common_rxns = setdiff(model_amides.rxns, rxn_unique);
model_amides = tncore_remove_reactions(model_amides, common_rxns);
model_amides = tncore_remove(model_amides);

% Get reaction formulas 
model_amides_formulas = printRxnFormula(model_amides, model_amides.rxns, false);
model_ureides_formulas = printRxnFormula(model_ureides, model_ureides.rxns, false);

% Get the relevant information from the models
rxnNameList = vertcat(model_amides.rxnNames, model_ureides.rxnNames);                        
rxnAbrList = vertcat(model_amides.rxns, model_ureides.rxns);                        
rxnList = vertcat(model_amides_formulas, model_ureides_formulas);
revFlagList = vertcat(model_amides.rev, model_ureides.rev);
lowerBoundList = vertcat(model_amides.lb, model_ureides.lb);
upperBoundList= vertcat(model_amides.ub, model_ureides.ub);
subSystemList= vertcat(model_amides.subSystems, model_ureides.subSystems);
grRuleList = vertcat(model_amides.grRules, model_ureides.grRules);
geneNameList = vertcat(model_amides.genes, model_ureides.genes);

% Combine the models
combinedModel = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, ...
    lowerBoundList, upperBoundList, subSystemList, grRuleList);

% Fix a few fields, if necessary
if size(combinedModel.lb,1) == 1
    combinedModel.lb = transpose(combinedModel.lb);
end
if size(combinedModel.ub,1) == 1
    combinedModel.ub = transpose(combinedModel.ub);
end
if size(combinedModel.c,1) == 1
    combinedModel.c = transpose(combinedModel.c);
end

% Set the objective
combinedModel = changeObjective(combinedModel, 'R_Biomass');

% Solve the model
sol = optimizeCbModel(combinedModel)

%% Save
save('combinedModel.mat', 'combinedModel');
save('allWorkspace.mat');
clear
