addpath /Users/bethanyh/Documents/TIGER
addpath /Users/bethanyh/Documents/
addpath /Users/bethanyh/Documents/Tn-Core-master/Most-Recent-Version
load('noduleModel.mat');
load('plantModel.mat');
%load('medicagoModel'); % are these needed?
%load('melilotiModel');
changeCobraSolver ('ibm_cplex');

%% Make a nodule model

noduleModel.rev=[];
for i=1:length(noduleModel.rxns)
    if noduleModel.lb(i) < 0
        noduleModel.rev(i,1)=1;
    else
        noduleModel.rev(i,1)=0;
    end
end
plantModel.rev=[];
for i=1:length(plantModel.rxns)
    if plantModel.lb(i) < 0
        plantModel.rev(i,1)=1;
    else
        plantModel.rev(i,1)=0;
    end
end

%% Delete all exchange reactions from the nodule zones

% Get exchange reactions
exchangeRxns = noduleModel.rxns(find(contains(noduleModel.rxns,'_tx')));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(find(contains( noduleModel.rxns,'EX_'))));
exchangeRxns = setdiff(exchangeRxns, {'R_PROTON_tx'; 'Nodule_R_MNXM111_tx'; 'Nodule_R_MNXM1026_tx'});

% Delete the exchange reactions
noduleModel_noEx = removeRxns(noduleModel, exchangeRxns);

%% Add the nodule to the plant

% Get reaction formulas for all models
noduleFormulas = printRxnFormula(noduleModel_noEx, noduleModel_noEx.rxns, false);
plantFormulas = printRxnFormula(plantModel, plantModel.rxns, false);

% Get the relevant information from all models
rxnNameList = vertcat(plantModel.rxnNames, noduleModel_noEx.rxnNames);                        
rxnAbrList = vertcat(plantModel.rxns, noduleModel_noEx.rxns);                        
rxnList = vertcat(plantFormulas, noduleFormulas);
revFlagList = vertcat(plantModel.rev, noduleModel_noEx.rev);
lowerBoundList = vertcat(plantModel.lb, noduleModel_noEx.lb);
upperBoundList= vertcat(plantModel.ub, noduleModel_noEx.ub);
subSystemList= vertcat(plantModel.subSystems, noduleModel_noEx.subSystems);
grRuleList = vertcat(plantModel.grRules, noduleModel_noEx.grRules);
geneNameList = vertcat(plantModel.genes, noduleModel_noEx.genes);

% Combine the models
nodulatedPlantDisconnected = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, ...
    lowerBoundList, upperBoundList, subSystemList, grRuleList);

nodulatedPlantDisconnected = changeRxnBounds(nodulatedPlantDisconnected, 'Bacteroid_ATPMR', 0, 'l');
nodulatedPlantDisconnected = changeObjective(nodulatedPlantDisconnected, 'R_Biomass');

%% Connect nodule metabolism to plant metabolism


% Update symport reactions
transportRxns0 = findRxnsFromMets(nodulatedPlantDisconnected, 'Bacteroid_MNXM1[e]');
transportRxns1=findRxnsFromMets(nodulatedPlantDisconnected, 'Bacteroid_h[c]');
transportRxns=intersect(transportRxns0,transportRxns1);
ext=find(contains(nodulatedPlantDisconnected.rxns,'EXCT'));
transportRxns2=nodulatedPlantDisconnected.rxns(ext);
transportRxns=setdiff(transportRxns,transportRxns2);
rxnsToKeep = {'Bacteroid_ATPS4r';'Bacteroid_NADH11';'Bacteroid_NO3RxM';'Bacteroid_NO3RxU';...
    'Bacteroid_MPH';'Bacteroid_COBDM';'Bacteroid_COBDU';'Bacteroid_COBOU';'Bacteroid_UPDH';...
    'Bacteroid_UPH';'Bacteroid_THD2'};
transportRxns = setdiff(transportRxns,rxnsToKeep);
transportRxnIDs = findRxnIDs(nodulatedPlantDisconnected, transportRxns);
nodulatedPlantDisconnected.rev(transportRxnIDs) = 0;
nodulatedPlantDisconnected.lb(transportRxnIDs) = 0;

% Add transfer reactions
connectNoduleAPR;

%% Set the objective function

% Make reaction to combine plant biomass
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'R_PlantBiomass', ...
    {'Leave_BiomassShoot[c]', 'Root_BiomassRoot[c]', 'BiomassPlant[c]'}, [-0.33333 -0.66667 1], 0, 0, 1000, 0);
makegrowingnodule_gd2

% Add biomass exchange reaction
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_Biomass', ...
    {'BiomassTotal[c]'}, [-1], 0, 0, 1000, 0);
%nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_Biomass', ...
%    {'Biomass[c]'}, [-1], 0, 0, 1000, 0); % When you uncomment the 'makegrowingnodule_gd' line, comment out this biomass reaction and uncomment the biomass reaction on the previous line
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_NoduleBiomass', ...
    {'BiomassNodule[c]'}, [-1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_PlantBiomass', ...
    {'BiomassPlant[c]'}, [-1], 0, 0, 1000, 0);

% Fix reversibility field
nodulatedPlantConnected.rev=[];
for i=1:length(nodulatedPlantConnected.rxns)
    if nodulatedPlantConnected.lb(i) < 0
        nodulatedPlantConnected.rev(i,1)=1;
    else
        nodulatedPlantConnected.rev(i,1)=0;
    end
end

% Set the objective function
nodulatedPlantConnected.c(:) = 0;
nodulatedPlantConnected = changeObjective(nodulatedPlantConnected, 'R_Biomass'); 

% Adjust nitrogen-related reactions
nodulatedPlantConnected = changeRxnMets(nodulatedPlantConnected, {'Leave_AMMONIA[c]'}, {'Leave_AMMONIUM[c]'}, 'TRS_AMMONIUM');
nodulatedPlantConnected = changeRxnMets(nodulatedPlantConnected, {'Root_AMMONIA[c]'}, {'Root_AMMONIUM[c]'}, 'TRS_AMMONIUM');
nodulatedPlantConnected = changeRxnMets(nodulatedPlantConnected, {'Root_AMMONIA[c]'}, {'Root_AMMONIUM[c]'}, 'TRN_AMMONIUM');
nodulatedPlantConnected = changeRxnMets(nodulatedPlantConnected, {'Root_AMMONIA[c]'}, {'Root_AMMONIUM[c]'}, 'R_NoduleBiomass');
nodulatedPlantConnected = changeRxnMets(nodulatedPlantConnected, {'Root_AMMONIA[c]'}, {'Root_AMMONIUM[c]'}, 'Root_R_NH4_tx');

nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Root_R_NH4_tx', 0, 'b');
nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Nodule_R_NH4_tx', 0, 'b');
nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Root_R_NO3_tx', 0, 'b');
nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Nodule_R_NO3_tx', 0, 'b');
nodulatedPlantConnected = removeRxns(nodulatedPlantConnected,{'Leave_R_RXN-11811_c','Nodule_R_RXN-11811_c'});
% nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Nodule_R_ALLANTOINASE-RXN_c', 1000, 'b');
% nodulatedPlantConnected = changeRxnBounds(nodulatedPlantConnected,'Nodule_R_ALLANTOATE-DEIMINASE-RXN_c', 1000, 'b');
% nodulatedPlantConnected.rev(1851)=1;
% nodulatedPlantConnected.rev(1853)=1;
% Adjust ASN in biomass
% Temporary - to be removed once biomass reaction is properly updated
RootBiomass = find(ismember(nodulatedPlantConnected.rxns,'Root_R_BiomassRoot'));
LeaveBiomass = find(ismember(nodulatedPlantConnected.rxns,'Leave_R_BiomassShoot'));
AsparagineLeave = find(ismember(nodulatedPlantConnected.mets,{'Leave_ASN[c]'}));
AsparagineRoot = find(ismember(nodulatedPlantConnected.mets,{'Root_ASN[c]'}));
ASNMolWeight = 132.12;
ASNLeaveAmounts = nodulatedPlantConnected.S(AsparagineLeave,LeaveBiomass);
AsnLeaveWeights = abs(ASNMolWeight/1e6 * ASNLeaveAmounts);
AsnLeaveWeightChange = AsnLeaveWeights - AsnLeaveWeights * 0.05;
nodulatedPlantConnected.S(AsparagineLeave,LeaveBiomass) = 0.05 * nodulatedPlantConnected.S(AsparagineLeave,LeaveBiomass);
nodulatedPlantConnected.S(:,LeaveBiomass) = nodulatedPlantConnected.S(:,LeaveBiomass) / (1-AsnLeaveWeightChange);
ASNRootAmounts = nodulatedPlantConnected.S(AsparagineRoot,RootBiomass);
AsnRootWeights = abs(ASNMolWeight/1e6 * ASNRootAmounts);
AsnRootWeightChange = AsnRootWeights - AsnRootWeights * 0.05;
nodulatedPlantConnected.S(AsparagineRoot,RootBiomass) = 0.05 * nodulatedPlantConnected.S(AsparagineRoot,RootBiomass);
nodulatedPlantConnected.S(:,RootBiomass) = nodulatedPlantConnected.S(:,RootBiomass) / (1-AsnRootWeightChange);

% Check solution
optimizeCbModel(nodulatedPlantConnected)





%% Add maintenance costs

% Root and shoot maintenance based on GrowthRateComparison (Thomas Pfau)
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000;
ATPMaintenanceShoot = find(ismember(nodulatedPlantConnected.rxns,'Leave_R_ADENOSINETRIPHOSPHATASE-RXN_c'));
ATPMaintenanceRoot = find(ismember(nodulatedPlantConnected.rxns,'Root_R_H+_ATPase_c'));
nodulatedPlantConnected.lb(ATPMaintenanceShoot) = MaintenanceRequirement * 2/3;
nodulatedPlantConnected.lb(ATPMaintenanceRoot) = MaintenanceRequirement * 1/3;

% Add nodule maintenance costs
MaintenanceRequirementNodule = MaintenanceRequirement * 0.02;
ATPMaintenanceZI = findRxnIDs(nodulatedPlantConnected, 'Nodule_R_H+_ATPase_c');
nodulatedPlantConnected.lb(ATPMaintenanceZI) = MaintenanceRequirementNodule * 0.05;

% Add bacteroid maintenance costs
MaintenanceRequirementBacteroid = 8400 * 0.02 * 0.3;
ATPMaintenanceZIId = findRxnIDs(nodulatedPlantConnected, 'Bacteroid_ATPM');
ATPMaintenanceZII = findRxnIDs(nodulatedPlantConnected, 'Bacteroid_ATPMR');
nodulatedPlantConnected.lb(ATPMaintenanceZIId) = MaintenanceRequirementBacteroid * 0.45 * 0.25;
nodulatedPlantConnected.lb(ATPMaintenanceZII) = MaintenanceRequirementBacteroid * 0.45 * 0.25;

% Switch name
nodulatedPlant = nodulatedPlantConnected;

%% Set boundaries to reasonable limit

% Change boundaries
for n = 1:length(nodulatedPlant.rxns)
    if nodulatedPlant.ub(n) >= 1000
        nodulatedPlant.ub(n) = 1000;
    end
    if nodulatedPlant.lb(n) <= -1000
        nodulatedPlant.lb(n) = -1000;
    end
end

%% Fix the rules and grRules
% Just in case

% Change model name
model = nodulatedPlant;

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Fix the grRules
for n = 1:length(model.grRules)
    if ~isempty(model.grRules{n})
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), ' or', ''),'AAA', ' or ');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), 'or ', ''),'AAA', ' or ');
    end
end

% Fix the rules
for n = 1:length(model.rules)
    if ~isempty(model.rules{n})
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(strrep(strrep(model.rules{n}, ' | ', 'AAA'), '| ', ''),'AAA', ' | ');
    end
end

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Change model name

nodulatedPlant = model;

%% Remove deadends

model = nodulatedPlant;

model.lb = model.lb * 1000;
model.ub = model.ub * 1000;
optimizeCbModel(model)
deadEndMetabolites = detectDeadEnds(model);
if ~isempty(deadEndMetabolites)
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(model);
        if isempty(deadEndMetabolites)
            test = 0;
        else
            deadEndMetabolites3 = cell(length(deadEndMetabolites),1);
            for n = 1:length(deadEndMetabolites)
                deadEndMetabolites2 = num2cell(deadEndMetabolites);
                deadEndMetabolites3{n,1} = ...
                    model.mets{deadEndMetabolites2{n,1},1};
            end
            [deadEndReactions] = ...
                findRxnsFromMets(model,deadEndMetabolites3);
            model = tncore_remove_reactions(model,deadEndReactions);
            clear deadEndMetabolites
            clear deadEndMetabolites2
        end
    end
    model = tncore_remove(model);
end
optimizeCbModel(model)

nodulatedPlantOriginal = nodulatedPlant;
nodulatedPlant = model;

%% Set up plant transport boundaries
% This code is mostly taken from GrowthRateComparison of Thomas Pfau

Reversibility
nodulatedPlant = model;

%% Save and clean workspace

save('allWorkspace.mat');
save('nodulatedPlant.mat', 'nodulatedPlant');
clear;
