function PRODUCESOY
addpath /Users/bethanyh/Documents/FASTCORE_1

%% Make tissue specific model

% Clear workspace
clear

% Make the model
plantModel=BUILDTISSUE();

% Fix grRules
plantModel = rmfield(plantModel, 'grRules');
plantModel = tncore_fix(plantModel);

%% Make tissue specific gene names

% Make models for each tissue
rootModel = plantModel;
shootModel = plantModel;

% Rename genes in each model
for n = 1:length(plantModel.genes)
    geneID = findGeneIDs(rootModel, plantModel.genes{n});
    rootModel.genes{geneID} = strcat('Root_', rootModel.genes{geneID});
    shootModel.genes{geneID} = strcat('Leave_', shootModel.genes{geneID});
end
for n = 1:length(plantModel.grRules)
    for m = 1:length(plantModel.genes)
        rootModel.grRules{n} = strrep(rootModel.grRules{n}, plantModel.genes{m}, rootModel.genes{m});
        shootModel.grRules{n} = strrep(shootModel.grRules{n}, plantModel.genes{m}, shootModel.genes{m});
    end
end

% Get list or root and shoot reactions
rootRxns = plantModel.rxns(strmatch('Root_', plantModel.rxns));
shootRxns = plantModel.rxns(strmatch('Leave_', plantModel.rxns));

% Get necessary info to add root to the shoot
rootRxnNames = plantModel.rxnNames(strmatch('Root_', plantModel.rxns));
plantFormulas = printRxnFormula(plantModel);
rootFormulas = plantFormulas(strmatch('Root_', plantModel.rxns));
rootRev = plantModel.rev(strmatch('Root_', plantModel.rxns));
rootLb = plantModel.lb(strmatch('Root_', plantModel.rxns));
rootUb = plantModel.ub(strmatch('Root_', plantModel.rxns));
rootC = plantModel.c(strmatch('Root_', plantModel.rxns));
rootSubsystem = plantModel.subSystems(strmatch('Root_', plantModel.rxns));
rootGrrules = rootModel.grRules(strmatch('Root_', plantModel.rxns));

% Delete root reactions from the shoot models
shootModel = tncore_remove_reactions(shootModel, rootRxns);

% Add root reactions back to the shoot model
rebuiltModel = shootModel;
for n = 1:length(rootRxns)
    rebuiltModel = addReaction(rebuiltModel, {rootRxns{n}, rootRxnNames{n}}, rootFormulas{n}, ...
        [], rootRev(n), rootLb(n), rootUb(n), rootC(n), rootSubsystem{n}, ...
        rootGrrules{n}, [], [], false);
end
rebuiltModel.metChEBIID = cell(1712,1);

% Remove any unused genes (if they exist)
plantModel = removeUnusedGenes(rebuiltModel);

% Set objective to biomass
plantModel = changeObjective(plantModel, 'R_Biomass');

% Test growth
final_sol = optimizeCbModel(plantModel)

%% Save and clean workspace

save('allWorkspace.mat');
save('plantModel.mat', 'plantModel');
clear;
end