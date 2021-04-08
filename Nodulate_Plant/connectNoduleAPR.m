%% Find reactions to connect nodule to root metabolism and the soil

% Identify reactions for import and export between root and soil
importLogical = cell(length(plantModel.rxns),1);
exportLogical = cell(length(plantModel.rxns),1);
for n = 1:length(plantModel.rxns)
    if strmatch('Root_', plantModel.rxns{n})
        if strfind(plantModel.rxns{n}, '_tx')
            importLogical{n,1} = 1;
            exportLogical{n,1} = 1;
            
        else
            importLogical{n,1} = 0;
            exportLogical{n,1} = 0;
            
        end
    else
        importLogical{n,1} = 0;
        exportLogical{n,1} = 0;
    end
end
importLogical = logical(cell2mat(importLogical));
importRxns = plantModel.rxns(importLogical);
importIDs = findRxnIDs(plantModel, importRxns);
exportLogical = logical(cell2mat(exportLogical));
exportRxns = plantModel.rxns(exportLogical);
exportIDs = findRxnIDs(plantModel, exportRxns);

% Identify reactions for exchange between root and shoot
TRS_reactions = plantModel.rxns(strmatch('TRS_', plantModel.rxns));
TSR_reactions = plantModel.rxns(strmatch('TSR_', plantModel.rxns));
transferRxns = vertcat(TRS_reactions, TSR_reactions);
transferIDs = findRxnIDs(plantModel, transferRxns);

%% Add export reactions to the nodule

% Set new model name
nodulatedPlantDisconnected_export = nodulatedPlantDisconnected;

% Add export reactions for zone III
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_R_', 'Nodule_R_');
    rxnFormula = printRxnFormula(noduleModel, rxnAbr);
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

%% Allow compounds imported by root to be imported by the nodule 
% Set new model name
nodulatedPlantDisconnected_import = nodulatedPlantDisconnected_export;

% Add import to general nodule from soil
for n = 1:length(importRxns)
    rxnAbr = strrep(importRxns{n}, 'Root_R_', 'TRN_');
    importRxns{n} = strrep(importRxns{n}, 'Root_R_', 'Nodule_R_');
    rxnFormula = printRxnFormula(noduleModel, importRxns{n});
    [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
    if strmatch('Nodule_R_Pi_tx', importRxns{n}, 'exact')
        metToAdd = {'Nodule_MNXM9[c]'};
    elseif strmatch('Nodule_R_PROTON_tx', importRxns{n}, 'exact')
        metToAdd = {'Nodule_MNXM1[c]'};
    else
        metToAdd = setdiff(metaboliteList, ...
            {'Nodule_MNXM9[c]', 'Nodule_ATP[c]', 'Nodule_ADP[c]', 'Nodule_MNXM1[c]','Nodule_Pumped-PROTON[c]','Nodule_h[c]'});
    end
    if strmatch('Nodule_MNXM1[c]', metToAdd, 'exact')
    else
        nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
            rxnAbr, 'metaboliteList', metToAdd, 'stoichCoeffList', [1], ...
            'reversible', 0, 'lowerBound', 0, 'upperBound', 1000, 'objectiveCoef', 0);
    end
end

% Prevent proton import in the root
nodulatedPlantDisconnected_import = changeRxnBounds(nodulatedPlantDisconnected_import, ...
    {'Root_R_PROTON_tx'}, 0, 'u');
%GOT UP TO HERE-------------------------


% % Add transfer from general nodule to zone III
% for n = 1:length(importRxns)
%     if strmatch('Root_R_Pi_tx', importRxns{n}, 'exact')
%         rxnAbr = strrep(importRxns{n}, 'Root_R', 'TNB_');
%         rxnFormula = 'Nodule_Pi[c] + 0.25 Bacteroid_atp[c] -> 0.25 Bacteroid_adp[c] + 0.25 Bacteroid_h[c] + 1.25 Bacteroid_pi[c]';
%     else
%         rxnAbr = strrep(importRxns{n}, 'Root_R_', 'TNB_');
%         rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
%         [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
%         metToAdd = setdiff(metaboliteList, {'Root_Pi[c]', 'Root_ATP[c]', 'Root_ADP[c]', 'Root_PROTON[c]','Root_Pumped-PROTON[c]','Root_h[c]'});
%         metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
%         rxnFormula = strrep(rxnFormula, 'Root_', 'Bacteroid_');
%         rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
%         rxnFormula = strrep(rxnFormula, '_____', ' ');
%     end
%     nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
%         rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
% end

%% Allow compounds transferred between root and shoot to be transfered to the nodule

% Set new model name
nodulatedPlantDisconnected_transfer = nodulatedPlantDisconnected_import;

% Load name conversions
metNames = table2cell(readtable('MappedCompoundsTrunca', 'ReadVariableNames', false));

% Add transfer to general nodule from root
for n = 1:length(transferRxns)
    if strmatch('TSR_PROTON', transferRxns{n})
        rxnAbr = 'TRN_PROTON';
        rxnFormula = 'Root_PROTON[c] -> Nodule_MNXM1[c]';
        if strmatch('TRN_PROTON', nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = 'TRN_2_PROTON';
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    elseif strmatch('TRS_PROTON', transferRxns{n})
        rxnAbr = 'TRN_PROTON';
        rxnFormula = 'Root_PROTON[c] -> Nodule_MNXM1[c]';
        if strmatch('TRN_PROTON', nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = 'TRN_2_PROTON';
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    elseif strmatch('TRS_Pi', transferRxns{n})
        rxnAbr = 'TRN_Pi';
        rxnFormula = 'Root_Pi[c] -> Nodule_MNXM9[c]';
        if strmatch('TRN_PROTON', nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = 'TRN_2_PROTON';
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    else
        unwantedMets = {'Leave_ATP[c]'; 'Leave_Pi[c]'; 'Leave_PROTON[c]';'Leave_h[c]'; 'Leave_ADP[c]'; ...
            'Root_ATP[c]'; 'Root_Pi[c]'; 'Root_PROTON[c]'; 'Root_ADP[c]';'Root_h[c]'};
        transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
            logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
            nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
        transferredMet = transferredMet(strmatch('Root_', transferredMet));
        rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TRN_');
        rxnAbr = strrep(rxnAbr, 'TRS_', 'TRN_');
        receivedMet = strrep(transferredMet, 'Root_', 'Nodule_');
        metToFind = strrep(receivedMet, 'Nodule_', '');
        metToFind = strrep(metToFind, '[c]', '_c');
        metPos = strmatch(metToFind, metNames(:,2));
        metToFind = unique(metNames(metPos,1));
        receivedMet = strcat('Nodule_', metToFind, '[c]');      
        rxnFormula = [cell2mat(transferredMet) ' -> ' cell2mat(receivedMet)];
        if strmatch(rxnAbr, nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = strrep(rxnAbr, 'TRN_', 'TRN_2_');
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    end
end


% Add transfer from general nodule to zone III
% for n = 1:length(transferRxns)
%     rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNB_');
%     rxnAbr = strrep(rxnAbr, 'TRS_', 'TNB_');
%     if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
%     else
%         if strmatch('TSR_PROTON', transferRxns{n})
%             rxnFormula = '0.25 Bacteroid_ATP[c] + Bacteroid_PROTON[c] -> 0.25 Bacteroid_ADP[c] + 0.25 Bacteroid_Pi[c] + 1.25 Bacteroid_PROTON[c]';
%             nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
%                 rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
%         else
%             unwantedMets = {'Leave_ATP[c]'; 'Leave_Pi[c]'; 'Leave_PROTON[c]'; 'Leave_h[c]';'Leave_ADP[c]'; ...
%                 'Root_ATP[c]'; 'Root_Pi[c]'; 'Root_PROTON[c]'; 'Root_ADP[c]';'Root_h[c]'};
%             transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
%                 logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
%                 nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
%             transferredMet = transferredMet(strmatch('Root_', transferredMet));
%             donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
%             receivedMet = strrep(transferredMet, 'Root_', 'Bacteroid_');
%             rxnFormula = ['0.25 Bacteroid_ATP[c] +_____' cell2mat(donorMet) ...
%                 ' -> 0.25 Bacteroid_ADP[C] + 0.25 Bacteroid_Pi[c] + 0.25 Bacteroid_PROTON[c] +_____' ...
%                 cell2mat(receivedMet)];
%             rxnFormula = strrep(rxnFormula, '+_____', '+ ');
%             nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
%                 rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
%         end
%     end
% end

%% Add transfer of ureides from nodule to root
 rxnFormula1 = 'Nodule_S-ALLANTOIN[c] -> Root_S-ALLANTOIN[c]';
 nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'TNR_S-ALLANTOIN',rxnFormula1, [], 0, 0, 1000, 0);
 rxnFormula2='Nodule_ALLANTOATE[c] -> Root_ALLANTOATE[c]' ;
 nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'TNR_ALLANTOATE',rxnFormula2, [], 0, 0, 1000, 0);
% rxnformula1='  0.25 Bacteroid_PROTON[c] + Bacteroid_AMMONIA[c] + 0.25 Bacteroid_Pi[c] + 0.25 Bacteroid_ADP[C] -> Nodule_AMMONIA[c] + 0.25 Bacteroid_ATP[c]' ;
% nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'TBN_AMMONIUM',rxnformula1, [], 0, 0, 1000, 0);
% rxnformula1=' Bacteroid_nh3[c] -> Root_AMMONIUM[c]' ;
% nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'TBR_AMMONIUM',rxnformula1, [], 0, 0, 1000, 0);
rxnformula3=' -> Nodule_MNXM724[c]' ;
 nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'Nodule_R_N2_tx',rxnformula3, [], 0, 0, 1000, 0);
rxnformula4='Bacteroid_n2[e] -> Nodule_MNXM724[c]' ;
 nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'EXCT_rev_MNXM724[e]',rxnformula4, [], 0, 0, 1000, 0);
rxnformula5='Nodule_MNXM724[c] -> Bacteroid_n2[e]' ;
 nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'EXCT_for_MNXM724[e]',rxnformula5, [], 0, 0, 1000, 0);
 %rxnformula8='Nodule_ALLANTOATE[c] + 2 Nodule_MNXM1[c] + Nodule_MNXM2[c] <=> Nodule_CPD0-2298[c] + Nodule_MNXM13[c] + Nodule_MNXM15[c]';  
  %nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'Nodule_R_ALLANTOATE-DEIMINASE-RXN_c',rxnformula8, [], 0, -1000, 1000, 0);
rxnformula9='Nodule_MNXM2[c] + Nodule_S-ALLANTOIN[c] <=> Nodule_ALLANTOATE[c] + Nodule_MNXM1[c]  ';
  nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'Nodule_R_ALLANTOINASE-RXN_c',rxnformula9, [], 0, -1000, 1000, 0);

 %rxnFormula6 = 'Nodule_MNXM15[c] -> Root_AMMONIA[c]';
 %nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'TNR_AMMONIUM',rxnFormula6, [], 0, 0, 1000, 0);
% rxnformula6='Bacteroid_MNXM169[e] -> Nodule_MNXM169[c]' ;Nodule_MNXM15[c]
% nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'EXCT_rev_MNXM169[e]',rxnformula6, [], 0, 0, 1000, 0);
% rxnformula7='Nodule_MNXM169[c] -> Bacteroid_MNXM169[e]' ;
% nodulatedPlantDisconnected_transfer=addReaction(nodulatedPlantDisconnected_transfer,'EXCT_for_MNXM169[e]',rxnformula7, [], 0, 0, 1000, 0);


% Rename the model
nodulatedPlantConnected = nodulatedPlantDisconnected_transfer;

%% Deal with carbon dioxide

% Add carbon dioxide entry to nodule
%nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TRN_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]'}, [1], 0, 0, 1000, 0);
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNI_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]', 'NoduleI_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIId_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIId_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIp_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIIp_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIZ_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIZ_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIII_CARBON-DIOXIDE', ...
%     {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIII_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);

%% Transfer asparagine and glutamine from zone III to the root
% 
% % Add asparagine transfer
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIIR_ASN', ...
%     {'NoduleIII_ASN[C]', 'NoduleIII_ATP[C]', 'NoduleIII_ADP[C]', ...
%     'NoduleIII_Pi[C]', 'NoduleIII_PROTON[C]', 'Root_ASN[C]'}, ...
%     [-1 -0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
% 
% % Add glutamine transfer
% nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIIR_GLN', ...
%     {'NoduleIII_GLN[C]', 'NoduleIII_ATP[C]', 'NoduleIII_ADP[C]', ...
%     'NoduleIII_Pi[C]', 'NoduleIII_PROTON[C]', 'Root_GLN[C]'}, ...
%     [-1 -0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);

%% Remove reaction adding ammonium to the nodule from the soil

nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'TRN_NH4_tx');
nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'TRN_NO3_tx');
nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'Nodule_R_BiomassShoot');
nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'Nodule_R_BiomassRoot');


