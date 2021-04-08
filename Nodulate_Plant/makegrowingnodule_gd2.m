%% Identify information from ViNE
changeCobraSolver ('ibm_cplex');
addpath /Users/bethanyh/Documents/FASTCORE_1
% Load ViNE
m = readCbModel('vin.mat');

% Get solution
sol = optimizeCbModel(m);

% Get metabolites and initial stoichiometries for nodule-transferred metabolites
tni_rxns = m.rxns(strmatch('TNI_', m.rxns));
tniid_rxns = m.rxns(strmatch('TNIId_', m.rxns));
tniip_rxns = m.rxns(strmatch('TNIIp_', m.rxns));
tniz_rxns = m.rxns(strmatch('TNIZ_', m.rxns));
all_rxns = vertcat(tni_rxns, tniid_rxns, tniip_rxns, tniz_rxns);
metsToTransfer = {};
stoichToTransfer = {};
for n = 1:length(all_rxns)
    if strfind(all_rxns{n}, '_PROTON')
        posFinal = strmatch('Nodule_MNXM1[C]', metsToTransfer, 'exact');
        if isempty(posFinal)
            metsToTransfer{end+1,1} = 'Nodule_MNXM1[C]';
            stoichToTransfer{end+1,1} = sol.x(findRxnIDs(m, all_rxns{n}));
        else
            stoichToTransfer{posFinal,1} = stoichToTransfer{posFinal,1} + sol.x(findRxnIDs(m, all_rxns{n}));
        end
    else
        pos = findRxnIDs(m, all_rxns{n});
        met = m.mets(logical(m.S(:,pos)));
        metPos = strmatch('Nodule_', met);
        met = met(metPos);
        posFinal = strmatch(met, metsToTransfer, 'exact');
        if isempty(posFinal)
            metsToTransfer{end+1,1} = met{1,1};
            stoichToTransfer{end+1,1} = sol.x(findRxnIDs(m, all_rxns{n}));
        else
            stoichToTransfer{posFinal,1} = stoichToTransfer{posFinal,1} + sol.x(findRxnIDs(m, all_rxns{n}));
        end
    end
end

% Normalize per gram of growing nodule tissue
for n = 1:length(stoichToTransfer)
    stoichToTransfer{n,1} = -1 * stoichToTransfer{n,1} / sol.x(findRxnIDs(m, 'NoduleBiomass'));
end

%% Add "fake" growing nodule
% Replace ViNE with the soybean model in this section

for n=1:length(metsToTransfer)
    metsToTransfer{n}=strrep(metsToTransfer{n},'Nodule_','');
    metsToTransfer{n}=strrep(metsToTransfer{n},'[C]','');
end
metacyc = table2cell(readtable('chem_xref_metacyc.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
for n = 1:length(metsToTransfer)
    met =  metsToTransfer{n};
    %met = strsplit(met, '[');
    %met = met{1};
    pos = strmatch(met, metacyc(:,2), 'exact');
    if n==3
        pos=pos(3);
    elseif n==4
        pos=pos(2);
    elseif n==5
        pos=pos(5);
    else
    pos=pos(1);
    end
    nim= strsplit(metacyc{pos}, 'metacyc:');
    nim=nim{2};
    %meta = [meta, nim]
    metsToTransfer{n}=strrep(metsToTransfer{n},met,nim)
     
end
metsToTransfer=strcat('Root_',metsToTransfer);
metsToTransfer=strcat(metsToTransfer,'[c]');
metsToTransfer(31)=[];
metsToTransfer(30)=[];
metsToTransfer(29)=[];
metsToTransfer(28)=[];
metsToTransfer(27)=[];
stoichToTransfer(31)=[];
stoichToTransfer(30)=[];
stoichToTransfer(29)=[];
stoichToTransfer(28)=[];
stoichToTransfer(27)=[];
 stoichToTransfer{26}=0; stoichToTransfer{24}=0;

% for n=1:length(metsToTransfer)
%     metsToTransfer{n}=strrep(metsToTransfer{n},'[C]','[c]');
%     metsToTransfer{n}=strrep(metsToTransfer{n},'Nodule_','Root_');
% end


% Add "fake" nodule biomass reaction
metsToTransfer{end+1,1} = 'BiomassNodule[c]';
stoichToTransfer{end+1,1} = 1;
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'R_NoduleBiomass', 'metaboliteList', metsToTransfer, ...
    'stoichCoeffList', cell2mat(stoichToTransfer), 'reversible', false, 'lowerBound', 0, ...
    'upperBound', 1000);

% Add nodule biomass to overall biomass reaction
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'R_Biomass', 'reactionFormula', ...
    ' 0.02 BiomassNodule[c] + 0.98 BiomassPlant[c] -> BiomassTotal[c]', ...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000);
