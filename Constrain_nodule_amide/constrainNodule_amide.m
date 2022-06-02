%% Load data
%addpath /Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox
%addpath /Users/bethanyh/Documents/

% addpath /Users/bethanyh/Documents/tigerNEW
% addpath(genpath('/Users/bethanyh/Documents/tigerNEW'))
% 
% addpath(genpath('/Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox'))
% addpath /Users/bethanyh/Documents/Tn-Core-master/Most-Recent-Version
% addpath /Users/bethanyh/Documents/FASTCORE_1
% changeCobraSolver ('ibm_cplex');

%% 

% Load previous files
load('temp_beforeGIMME.mat');

% Add amide export
nodulatedPlant = addReaction(nodulatedPlant,'TNR_GLN', ...
    '0.25 Nodule_ATP[c] + Nodule_MNXM37[c] + 0.25 Root_ATP[c] -> 0.25 Nodule_ADP[c] + 0.25 Nodule_MNXM1[c] + 0.25 Nodule_MNXM9[c] + 0.25 Root_ADP[c] + Root_GLN[c] + 0.25 Root_PROTON[c] + 0.25 Root_Pi[c]', ...
    [], 0, 0, 1000000, 0);
nodulatedPlant = addReaction(nodulatedPlant,'TNR_ASN', ...
    '0.25 Nodule_ATP[c] + Nodule_MNXM147[c] + 0.25 Root_ATP[c] -> 0.25 Nodule_ADP[c] + 0.25 Nodule_MNXM1[c] + 0.25 Nodule_MNXM9[c] + 0.25 Root_ADP[c] + Root_ASN[c] + 0.25 Root_PROTON[c] + 0.25 Root_Pi[c]', ...
    [], 0, 0, 1000000, 0);

% Remove ureide export
nodulatedPlant = removeRxns(nodulatedPlant, 'TNR_S-ALLANTOIN');
nodulatedPlant = removeRxns(nodulatedPlant, 'TNR_ALLANTOATE');

% Turn on tiger
start_tiger('cplex');

% Produce updated tiger model
nodulatedPlant_tiger = cobra_to_tiger(nodulatedPlant);

% Run the modified version of gimme
[geneStates, genesOut, sol, tiger, weightStruct] = tncore_multi_gimme(nodulatedPlant_tiger, ...
    fields, expressStruct, genesStruct, threshStruct, 0.99);

save('temp_afterGIMME_amide.mat');

%% Build context specific model in COBRA based on active rxns %% UP TO HERE

buildContextModelNEW;

%% Fix the rules and grRules
% Just in case, may not be necessary

% Change model name
nodulatedPlant = finalNodulatedPlant;

% Fix the grRules string formatting
for n = 1:length(nodulatedPlant.grRules)
    if isstring(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = str2mat(nodulatedPlant.grRules{n});
    else
        nodulatedPlant.grRules{n} = nodulatedPlant.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(nodulatedPlant.rules)
    if isstring(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = str2mat(nodulatedPlant.rules{n});
    else
        nodulatedPlant.rules{n} = nodulatedPlant.rules{n};
    end
end

% Fix the grRules
for n = 1:length(nodulatedPlant.grRules)
    if ~isempty(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(strrep(strrep(nodulatedPlant.grRules{n}, ' or ', 'AAA'), ' or', ''),'AAA', ' or ');
        nodulatedPlant.grRules{n} = strrep(strrep(strrep(nodulatedPlant.grRules{n}, ' or ', 'AAA'), 'or ', ''),'AAA', ' or ');
    end
end

% Fix the rules
for n = 1:length(nodulatedPlant.rules)
    if ~isempty(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(strrep(strrep(nodulatedPlant.rules{n}, ' | ', 'AAA'), '| ', ''),'AAA', ' | ');
    end
end

% Fix the grRules string formatting
for n = 1:length(nodulatedPlant.grRules)
    if isstring(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = str2mat(nodulatedPlant.grRules{n});
    else
        nodulatedPlant.grRules{n} = nodulatedPlant.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(nodulatedPlant.rules)
    if isstring(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = str2mat(nodulatedPlant.rules{n});
    else
        nodulatedPlant.rules{n} = nodulatedPlant.rules{n};
    end
end

% Change model name
constrainedNodule = nodulatedPlant;

%% Save and clean workspace

save('allWorkspace_amide.mat');
save('constrainedNodule_amide.mat', 'constrainedNodule');
clear;

