%% Load data
%addpath /Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox
%addpath /Users/bethanyh/Documents/

addpath /Users/bethanyh/Documents/tigerNEW
addpath(genpath('/Users/bethanyh/Documents/tigerNEW'))

addpath(genpath('/Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox'))
addpath /Users/bethanyh/Documents/Tn-Core-master/Most-Recent-Version
addpath /Users/bethanyh/Documents/FASTCORE_1
changeCobraSolver ('ibm_cplex');

%% adding ATP requirements for transfer reactions

%%
nodulatedPlant=readCbModel('nodulatedPlant.mat');
%nodulatedPlant=readCbModel('WorkingSoy.mat');
biomassrxns
ureides
pumpedprotonfix 
transportATP

%rna_sino=readcell('bradyexp.xlsx');
rna_sino = table2cell(readtable('bradyexp.txt'));

load('NoduleEXP.mat');
rna_medic=reactionNames;
% rna_medic = table2cell(readtable('Medicago_TPM_values_average.txt')); % assuming thats the mean of expression data for each gene
% rna_sino = table2cell(readtable('Sinorhizobium_TPM_values_average.txt'));
% rna_medic_full = table2cell(readtable('Medicago_TPM_values.txt'));
% rna_sino_full = table2cell(readtable('Sinorhizobium_TPM_values.txt'));
% rna_medic_original = rna_medic;
% rna_sino_original = rna_sino;


%% Increase flux

% nodulatedPlant.ub = nodulatedPlant.ub * 1000;
% nodulatedPlant.lb = nodulatedPlant.lb * 1000;
optimizeCbModel(nodulatedPlant)

%% Set model constraints
% 
% % Determine O2 consumption limit of nodule zone III   WHY???
% nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TRN_OXYGEN-MOLECULE')) = 1000 * 649 * 0.02;
% sol = optimizeCbModel(nodulatedPlant);
% nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TNIII_OXYGEN-MOLECULE')) = ...
%     sol.x(findRxnIDs(nodulatedPlant, 'TNIII_OXYGEN-MOLECULE'));
% nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TRN_OXYGEN-MOLECULE')) = 1000000;
% optimizeCbModel(nodulatedPlant)

% Remove unwanted NoduleIII_EXCT reactions
toRemove = {'EXCT_for_MNXM167';'EXCT_for_MNXM621';...
    'EXCT_for_MNXM1503';'EXCT_for_MNXM198';...
    'EXCT_for_MNXM165';'EXCT_for_MNXM468';...
    'EXCT_for_MNXM615'};
nodulatedPlant = tncore_remove_reactions(nodulatedPlant, toRemove);
sol = optimizeCbModel(nodulatedPlant);
reactions = nodulatedPlant.rxns(strmatch('EXCT_for', nodulatedPlant.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(nodulatedPlant, reactions{n});
    if abs(sol.x(pos)) < 0.00001
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'EXCT_for_MNXM25';'EXCT_for_MNXM98';'EXCT_for_MNXM93'});
nodulatedPlant = tncore_remove_reactions(nodulatedPlant, toRemove);
optimizeCbModel(nodulatedPlant)

%% Set growth rate threshold

sol = optimizeCbModel(nodulatedPlant);
growthThresh = 0.99 * sol.f;

%% Prepare the RNA-seq data 

processRNA;

%% Constrain the nodule 

% Turn on tiger
start_tiger('cplex');
% additional code to get cobra_to_tiger to work

nodulatedPlant.grRules{3126} = '( Bacteroid_blr2485 and Bacteroid_blr2486 and ( Bacteroid_bll0282 or Bacteroid_blr3729 ) and ( Bacteroid_bll0283 or Bacteroid_blr3728 ))';
nodulatedPlant.grRules{3127} = '( Bacteroid_blr2485 and Bacteroid_blr2486 and ( Bacteroid_bll0282 or Bacteroid_blr3729 ) and ( Bacteroid_bll0283 or Bacteroid_blr3728 ))';
nodulatedPlant.grRules{2712} = '( Bacteroid_blr2485 and Bacteroid_blr2486 and ( Bacteroid_bll0282 or Bacteroid_blr3729 ) and ( Bacteroid_bll0283 or Bacteroid_blr3728 ))';
nodulatedPlant.grRules{2713} = '( Bacteroid_blr2485 and Bacteroid_blr2486 and ( Bacteroid_bll0282 or Bacteroid_blr3729 ) and ( Bacteroid_bll0283 or Bacteroid_blr3728 ))';
nodulatedPlant.rules{2712} = '( x(679) & x(680) & ( x(22) | x(809) ) & ( x(23) | x(808) ))';
nodulatedPlant.rules{2713} = '( x(679) & x(680) & ( x(22) | x(809) ) & ( x(23) | x(808) ))';
nodulatedPlant.rules{3126} = '( x(679) & x(680) & ( x(22) | x(809) ) & ( x(23) | x(808) ))';
nodulatedPlant.rules{3127} = '( x(679) & x(680) & ( x(22) | x(809) ) & ( x(23) | x(808) ))';
 
% Produce tiger model
nodulatedPlant_tiger = cobra_to_tiger(nodulatedPlant);

% Split the rnaseq files
medicago_rnaseq = cell2mat(medicago_rnaseq_data_model(:,2));
meliloti_rnaseq = cell2mat(meliloti_rnaseq_data_model(:,2));
medicago_genes = medicago_rnaseq_data_model(:,1);
meliloti_genes = meliloti_rnaseq_data_model(:,1);

% % Force FixNOQP to be on
meliloti_rnaseq(strmatch('Bacteroid_sma0767', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('Bacteroid_sma0767', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('Bacteroid_sma1214', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('Bacteroid_sma1214', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('Bacteroid_sma0767', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('Bacteroid_sma0767', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('Bacteroid_sma1214', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('Bacteroid_sma1214', meliloti_genes)) * 10;

% Prepare the input for tncore_multiGIMME
fields = {'Medicago'; 'Meliloti'};
expressStruct = struct();
genesStruct = struct();
threshStruct = struct();
expressStruct.Medicago = medicago_rnaseq;
genesStruct.Medicago = medicago_genes;
threshStruct.Medicago = medic_thresh;
expressStruct.Meliloti = meliloti_rnaseq;
genesStruct.Meliloti = meliloti_genes;
threshStruct.Meliloti = sino_thresh;

save('temp_beforeGIMME.mat');

% Run the modified version of gimme
[geneStates, genesOut, sol, tiger, weightStruct] = tncore_multi_gimme(nodulatedPlant_tiger, ...
    fields, expressStruct, genesStruct, threshStruct, 0.99);

save('temp_afterGIMME.mat');

%% Build context specific model in COBRA based on active rxns %% UP TO HERE

buildContextModelNEW;

%% Fix the rules and grRules
% Just in case, may not be necessary

% Change model name
model = finalNodulatedPlant;

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
constrainedNodule = model;

%% Save and clean workspace

save('allWorkspace.mat');
save('constrainedNodule.mat', 'constrainedNodule');
clear;

