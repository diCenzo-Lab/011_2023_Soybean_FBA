% reduced from 3h to 1
model = addReaction(model, 'Leave_MNXR144933', 'reactionFormula', ...
    'Leave_MNXM89795[c] + Leave_MNXM11[c] <=> Leave_MNXM336[c] +  Leave_MNXM1[c] + Leave_MNXM121[c]',...
    'reversible', true, 'lowerBound', -1000000, 'upperBound', 1000000);
% no HIS on nodule side
model = addReaction(model, 'TRN_HIS', 'reactionFormula', ...
    ['0.25 Nodule_MNXM3[c] + 0.25 Root_MNXM3[c] + Root_MNXM134[c] -> 0.25 Nodule_MNXM7[c] + 0.25 Nodule_MNXM1[c] ' ...
    '  + 0.25 Nodule_MNXM9[c] + 0.25 Root_MNXM7[c] + 0.25 Root_MNXM1[c] + 0.25 Root_MNXM9[c] + Nodule_MNXM134[c] '],...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
% missing nodule GLY
model = addReaction(model, 'TRN_GLY', 'reactionFormula', ...
    ['0.25 Nodule_MNXM3[c] + 0.25 Root_MNXM3[c] + Root_MNXM29[c] -> 0.25 Nodule_MNXM7[c] + 0.25 Nodule_MNXM1[c] ' ...
    '+ 0.25 Nodule_MNXM9[c] + 0.25 Root_MNXM7[c] + 0.25 Root_MNXM1[c] + 0.25 Root_MNXM9[c] + Nodule_MNXM29[c]'],...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
% missing Bacteroid H
model = addReaction(model, 'Bacteroid_rxnGFTPcpd00919', 'reactionFormula', ...
    ['Bacteroid_MNXM722779[e] + Bacteroid_MNXM3[c] + Bacteroid_MNXM2[c] -> Bacteroid_MNXM7[c] + 2 Bacteroid_MNXM1[c] ' ...
    '  + Bacteroid_MNXM5782[c] + Bacteroid_MNXM9[c] '],'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
% extra h
rxn=find(contains(model.rxns,'Bacteroid_MNXR146573'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0
% missing h20
rxn=find(contains(model.rxns,'Bacteroid_MNXR104650'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=-3
% extra h
rxn=find(contains(model.rxns,'Bacteroid_MNXR104868'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_MNXR146556'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_PPCOAOR'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;
% missing Bacteroid H
model = addReaction(model, 'Bacteroid_MNXR145285', 'reactionFormula', ...
    ['Bacteroid_MNXM2[c] + Bacteroid_MNXM11[c] -> 2 Bacteroid_MNXM9[c] + Bacteroid_MNXM1[c]'],'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
% missing Bacteroid H
model = addReaction(model, 'Bacteroid_MNXR102507', 'reactionFormula', ...
    ['Bacteroid_MNXM3[c] + Bacteroid_MNXM89621[c] -> Bacteroid_MNXM7[c] + Bacteroid_MNXM90006[c] + Bacteroid_MNXM1[c]'],'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
% extra h
rxn=find(contains(model.rxns,'Bacteroid_OSUCCCOL'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=-1;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_MNXR101660'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_MNXR101657'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_HSCH'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=-1;
% extra h
rxn=find(contains(model.rxns,'Bacteroid_HC3E'));
pos = find(ismember(model.mets,'Bacteroid_MNXM1[c]'));
model.S(pos,rxn)=0;



% extra h20
rxn=find(contains(model.rxns,'TRS_SULFATE'));
pos = find(ismember(model.mets,'Bacteroid_MNXM2[c]'));
model.S(pos,rxn)=0;
% extra h20
rxn=find(contains(model.rxns,'TRS_SULFATE'));
pos = find(ismember(model.mets,'Leave_MNXM2[c]'));
model.S(pos,rxn)=-0.25;
rxn=find(contains(model.rxns,'TRS_SULFATE'));
pos = find(ismember(model.mets,'Root_MNXM2[c]'));
model.S(pos,rxn)=-0.25;

%% Fix protons

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Change proton count when obvious
for n = 1:length(massChargeBalancing_2)
    test_1 = str2num(strrep(massChargeBalancing_2{n,2}, ' H', ''));
    if test_1 == massChargeBalancing_2{n,3}
        reactionName = split(massChargeBalancing_2{n,1}, '_');
        if strmatch('Bacteroid', reactionName{1,1})
            metaboliteToAdd = 'Bacteroid_MNXM1[c]';
        else
            metaboliteToAdd = strcat(reactionName{1,1}, '_MNXM1[', reactionName{end,1}, ']');
        end
        reactionID = findRxnIDs(model, massChargeBalancing_2{n,1});
        metaboliteID = findMetIDs(model, metaboliteToAdd);
        model.S(metaboliteID, reactionID) = model.S(metaboliteID, reactionID) - massChargeBalancing_2{n,3};
    end
end

%% Fix other things

model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-GDP-FORMING-RXN_c')) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-GDP-FORMING-RXN_c')) - 1;
model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-UDP-FORMING-RXN_c')) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-UDP-FORMING-RXN_c')) - 1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-desaturase_p')) - 1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-hydrolase_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-hydrolase_p')) + 1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_MNXR145556_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_MNXR145556_p')) - 2;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_MNXR145556_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_MNXR145556_p')) + 2;


model.metFormulas{findMetIDs(model, 'Leave_MNXM925[p]')} = 'H25C14O8PN3S';
model.metFormulas{findMetIDs(model, 'Leave_MNXM565[p]')} = 'C14H18N5O11P';
model.metFormulas{findMetIDs(model, 'Leave_MNXM794[p]')} = 'C18H30O2';
model.metFormulas{findMetIDs(model, 'Leave_STARCH[p]')} = 'C6H10O5';
model.metFormulas{findMetIDs(model, 'Leave_MNXM29072[p]')} = 'C30H53N3O9PS';
model.metFormulas{findMetIDs(model, 'Leave_Linolenoyl-ACPs[p]')} = 'C32H54N3O9PS';
model.metFormulas{findMetIDs(model, 'Leave_PLASTOQUINOL-1[p]')} = 'C53H82O2';
model.metFormulas{findMetIDs(model, 'Leave_Photon[p]')} = '';
model.metFormulas{findMetIDs(model, 'Leave_XYLAN[c]')} = 'C5H8O4';
model.metFormulas{findMetIDs(model, 'Leave_MNXM184[p]')} = 'C3H2O3S';


model.metCharges(findMetIDs(model, 'Leave_MNXM169[p]')) = 1;
model.metCharges(findMetIDs(model, 'Leave_MNXM178[p]')) = 2;
model.metCharges(findMetIDs(model, 'Leave_MNXM230[m]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_Pumped-PROTON[m]')) = 1;
model.metCharges(findMetIDs(model, 'Leave_MNXM746[m]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM74[p]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM74[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM2229[m]')) = 0;
model.metCharges(findMetIDs(model, 'Leave_MNXM53428[m]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM184[p]')) = 0;
model.metCharges(findMetIDs(model, 'Leave_MNXM1223[p]')) = 0;
model.metCharges(findMetIDs(model, 'Leave_MNXM925[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM8467[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM96993[c]')) = 0;
model.metCharges(findMetIDs(model, 'Leave_MNXM89576[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM162942[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90751[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM162827[c]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM90665[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM162827[c]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM163775[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90752[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM167414[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_UDP-GLUCOSE[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM794[p]')) = 0;
model.metCharges(findMetIDs(model, 'Leave_MNXM162355[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM90012[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM169[c]')) = 1;
model.metCharges(findMetIDs(model, 'Leave_MNXM178[c]')) = 2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90886[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM90886[p]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM71[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM10916[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90340[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM162943[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM12847[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM29072[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM90878[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM167415[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90879[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM16648[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90880[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM8390[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM90881[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM164671[c]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM90753[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM162944[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM8977[p]')) = 2;
model.metCharges(findMetIDs(model, 'Leave_MNXM9034[p]')) = 1;
model.metCharges(findMetIDs(model, 'Leave_MNXM90667[c]')) = -3;
model.metCharges(findMetIDs(model, 'Leave_MNXM167416[c]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM6439[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM1223[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM2645[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM25602[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM184[p]')) = -2;
model.metCharges(findMetIDs(model, 'Leave_MNXM12845[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM24502[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM23683[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM28031[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM10019[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM979[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM26616[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM7127[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM23766[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM1473[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM28933[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM12844[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM23842[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM26095[p]')) = -1;
model.metCharges(findMetIDs(model, 'Leave_MNXM5723[p]')) = -1;

%% Fix protons again

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Change proton count when obvious
for n = 1:length(massChargeBalancing_2)
    test_1 = str2num(strrep(massChargeBalancing_2{n,2}, ' H', ''));
    if test_1 == massChargeBalancing_2{n,3}
        reactionName = split(massChargeBalancing_2{n,1}, '_');
        if strmatch('Bacteroid', reactionName{1,1})
            metaboliteToAdd = 'Bacteroid_MNXM1[c]';
        else
            metaboliteToAdd = strcat(reactionName{1,1}, '_MNXM1[', reactionName{end,1}, ']');
        end
        reactionID = findRxnIDs(model, massChargeBalancing_2{n,1});
        metaboliteID = findMetIDs(model, metaboliteToAdd);
        model.S(metaboliteID, reactionID) = model.S(metaboliteID, reactionID) - massChargeBalancing_2{n,3};
    end
end

%% Fix other things

model.metCharges(findMetIDs(model, 'Leave_Linoleoyl-ACPs[p]')) = -1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linolenoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linolenoyl-ACP-desaturase_p')) - 1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-desaturase_p')) + 1;
model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-hydrolase_p')) = model.S(findMetIDs(model, 'Leave_MNXM1[p]'), findRxnIDs(model, 'Leave_Linoleoyl-ACP-hydrolase_p')) - 1;

%% Automate fixing of charges of ACPs

for n = 1:length(model.mets)
    if isnan(model.metCharges(n))
        if ~isempty(cell2mat(strfind(model.metNames(n), 'ACPs[')))
            model.metCharges(n) = -1;
        end
    end
end

%% Fix other things

model.metFormulas{findMetIDs(model, 'Leave_Enzyme-N6-dihydrolipoyl-L-lysine[m]')} = 'C15H28N2O2S2';

%% Automate fixing of charges of tRNAs

for n = 1:length(model.mets)
    if isnan(model.metCharges(n))
        if ~isempty(cell2mat(strfind(model.metNames(n), '-tRNAs')))
            if isempty(cell2mat(strfind(model.metNames(n), 'Charged')))
                model.metCharges(n) = -3;
            else
                model.metCharges(n) = -2;
            end
        end
    end
end

%% Fix protons again

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Change proton count when obvious
for n = 1:length(massChargeBalancing_2)
    test_1 = str2num(strrep(massChargeBalancing_2{n,2}, ' H', ''));
    if test_1 == massChargeBalancing_2{n,3}
        reactionName = split(massChargeBalancing_2{n,1}, '_');
        if strmatch('Bacteroid', reactionName{1,1})
            metaboliteToAdd = 'Bacteroid_MNXM1[c]';
        else
            metaboliteToAdd = strcat(reactionName{1,1}, '_MNXM1[', reactionName{end,1}, ']');
        end
        reactionID = findRxnIDs(model, massChargeBalancing_2{n,1});
        metaboliteID = findMetIDs(model, metaboliteToAdd);
        model.S(metaboliteID, reactionID) = model.S(metaboliteID, reactionID) - massChargeBalancing_2{n,3};
    end
end

%% Fix root-shoot transfer reactions

for n = 1:length(model.rxns)
    if strmatch('TSR_WATER', model.rxns{n})
    elseif strmatch('TRS_WATER', model.rxns{n})
    elseif strmatch('TRS_SULFATE', model.rxns{n})
    elseif strmatch('TSR_', model.rxns{n})
        model.S(findMetIDs(model, 'Leave_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), n) - 0.25;
        model.S(findMetIDs(model, 'Root_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Root_MNXM2[c]'), n) - 0.25;
    elseif strmatch('TRS_', model.rxns{n})
        model.S(findMetIDs(model, 'Leave_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), n) - 0.25;
        model.S(findMetIDs(model, 'Root_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Root_MNXM2[c]'), n) - 0.25;
    end
end

%% Fix other things

model.metCharges(findMetIDs(model, 'Leave_BiomassShoot[c]')) = 0;
model.metCharges(findMetIDs(model, 'BiomassPlant[c]')) = 0;
model.metCharges(findMetIDs(model, 'BiomassNodule[c]')) = 0;
model.metCharges(findMetIDs(model, 'BiomassTotal[c]')) = 0;
model.metCharges(findMetIDs(model, 'Root_MNXM230[m]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM169[p]')) = 1;
model.metCharges(findMetIDs(model, 'Root_MNXM178[p]')) = 2;
model.metCharges(findMetIDs(model, 'Root_MNXM746[m]')) = -2;
model.metCharges(findMetIDs(model, 'Root_MNXM74[p]')) = -2;
model.metCharges(findMetIDs(model, 'Root_MNXM74[c]')) = -2;
model.metCharges(findMetIDs(model, 'Root_MNXM2229[m]')) = 0;
model.metCharges(findMetIDs(model, 'Root_MNXM53428[m]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM184[p]')) = 0;
model.metCharges(findMetIDs(model, 'Root_MNXM1223[p]')) = 0;
model.metFormulas{findMetIDs(model, 'Root_MNXM925[p]')} = 'H25C14O8PN3S';
model.metCharges(findMetIDs(model, 'Root_MNXM925[p]')) = -1;
model.metFormulas{findMetIDs(model, 'Root_MNXM565[p]')} = 'C14H18N5O11P';
model.metCharges(findMetIDs(model, 'Root_MNXM96993[c]')) = 0;
model.metCharges(findMetIDs(model, 'Root_MNXM162827[c]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM167414[c]')) = -3;
model.metFormulas{findMetIDs(model, 'Root_Linolenoyl-ACPs[p]')} = 'C32H54N3O9PS';
model.metFormulas{findMetIDs(model, 'Root_MNXM4897[c]')} = 'C6H10O5';
model.metCharges(findMetIDs(model, 'Root_UDP-GLUCOSE[c]')) = -2;
model.metCharges(findMetIDs(model, 'Root_MNXM164671[c]')) = -1;
model.metCharges(findMetIDs(model, 'Root_Linoleoyl-ACPs[p]')) = -1;
model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linolenoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linolenoyl-ACP-desaturase_p')) - 1;
model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-desaturase_p')) + 1;
model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-hydrolase_p')) = model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-hydrolase_p')) - 1;
model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-desaturase_p')) = model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-desaturase_p')) - 1;
model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-hydrolase_p')) = model.S(findMetIDs(model, 'Root_MNXM1[p]'), findRxnIDs(model, 'Root_Linoleoyl-ACP-hydrolase_p')) + 1;
model.metFormulas{findMetIDs(model, 'Root_MNXM184[p]')} = 'C3H2O3S';
model.metCharges(findMetIDs(model, 'Root_MNXM1223[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM184[p]')) = -2;
model.metCharges(findMetIDs(model, 'Root_MNXM169[c]')) = 1;
model.metCharges(findMetIDs(model, 'Root_MNXM178[c]')) = 2;
model.metCharges(findMetIDs(model, 'Root_MNXM10019[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM979[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM26616[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM7127[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM23766[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM1473[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM28933[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM12844[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM23842[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM26095[p]')) = -1;
model.metCharges(findMetIDs(model, 'Root_MNXM5723[p]')) = -1;
model.metFormulas{findMetIDs(model, 'Root_Enzyme-N6-dihydrolipoyl-L-lysine[m]')} = 'C15H28N2O2S2';
model.metCharges(findMetIDs(model, 'Root_BiomassRoot[c]')) = 0;

%% Fix protons again

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Change proton count when obvious
for n = 1:length(massChargeBalancing_2)
    test_1 = str2num(strrep(massChargeBalancing_2{n,2}, ' H', ''));
    if test_1 == massChargeBalancing_2{n,3}
        reactionName = split(massChargeBalancing_2{n,1}, '_');
        if strmatch('Bacteroid', reactionName{1,1})
            metaboliteToAdd = 'Bacteroid_MNXM1[c]';
        else
            metaboliteToAdd = strcat(reactionName{1,1}, '_MNXM1[', reactionName{end,1}, ']');
        end
        reactionID = findRxnIDs(model, massChargeBalancing_2{n,1});
        metaboliteID = findMetIDs(model, metaboliteToAdd);
        model.S(metaboliteID, reactionID) = model.S(metaboliteID, reactionID) - massChargeBalancing_2{n,3};
    end
end

%% Fix other things

model.metCharges(findMetIDs(model, 'Nodule_MNXM230[m]')) = -1;
model.metCharges(findMetIDs(model, 'Nodule_MNXM74[p]')) = -2;
model.metCharges(findMetIDs(model, 'Nodule_MNXM74[c]')) = -2;
model.metCharges(findMetIDs(model, 'Nodule_MNXM746[m]')) = -2;
model.metFormulas{findMetIDs(model, 'Nodule_MNXM565[p]')} = 'C14H18N5O11P';
model.metFormulas{findMetIDs(model, 'Nodule_STARCH[p]')} = 'C6H10O5';
model.metCharges(findMetIDs(model, 'Nodule_MNXM96993[c]')) = 0;
model.metCharges(findMetIDs(model, 'Nodule_MNXM722779[c]')) = -1;

%% Get metabolites from original Bradyrhizobium model

for n = 1:length(model.mets)
    if isnan(model.metCharges(n))
        if ~isempty(cell2mat(strfind(model.mets(n), 'Bacteroid')))
            if isempty(cell2mat(strfind(model.mets(n), 'MNXM')))
                met = strrep(model.mets{n}, 'Bacteroid_', '');
                pos = findMetIDs(brady, met);
                model.metCharges(n) = brady.metCharges(pos);
                model.metFormulas{n} = brady.metFormulas{pos};
            end
        end
    end
end

%% Fix other things

model.metFormulas{findMetIDs(model, 'Bacteroid_pnphis__L[c]')} = 'C7H9N4O5PR';
model.mets{findMetIDs(model, 'Bacteroid_pnphis__L[c]')} = 'Bacteroid_MNXM90093[c]';
model.mets{findMetIDs(model, 'Bacteroid_his__L__p[c]')} = 'Bacteroid_MNXM3704[c]';
model.metFormulas{findMetIDs(model, 'Bacteroid_MNXM3704[c]')} = 'C7H8N4O2R';
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM90093[c]')) = 0;
model = addReaction(model, 'Bacteroid_rxnSymCoF', 'reactionFormula', ...
    ['Bacteroid_MNXM3[c] + Bacteroid_cpdFe4S4[c] + Bacteroid_MNXM151[c] + Bacteroid_MNXM1[c] + Bacteroid_MNXM5782[c] + Bacteroid_MNXM1026[c] + Bacteroid_MNXM249[c] + Bacteroid_MNXM419[c] 	->	Bacteroid_symCoF[c] + Bacteroid_MNXM7[c] + Bacteroid_MNXM1[c] + Bacteroid_MNXM9[c]'],...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
model = removeRxns(model, 'Bacteroid_MNXR96123');
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM74[c]')) = -2;
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM169[e]')) = 1;
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM178[e]')) = 2;
model = addReaction(model, 'Bacteroid_FDXNH', 'reactionFormula', ...
    ['2 Bacteroid_MNXM169[e] + 2 Bacteroid_MNXM1[e] 	->	2 Bacteroid_MNXM178[e] + Bacteroid_MNXM195[e]'],...
    'reversible', false, 'lowerBound', 0, 'upperBound', 1000000);
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM169[c]')) = 1;
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM178[c]')) = 2;
model.metCharges(findMetIDs(model, 'Bacteroid_osucc[c]')) = -3;
model.mets{findMetIDs(model, 'Bacteroid_osucc[c]')} = 'Bacteroid_MNXM1093495[c]';
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM96993[c]')) = 0;
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM722779[e]')) = -1;
model.metCharges(findMetIDs(model, 'Bacteroid_MNXM5782[c]')) = -1;
model.metFormulas{findMetIDs(model, 'Bacteroid_MNXM5782[c]')} = 'C7H9O7';

%% Fix protons again

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Change proton count when obvious
for n = 1:length(massChargeBalancing_2)
    test_1 = str2num(strrep(massChargeBalancing_2{n,2}, ' H', ''));
    if test_1 == massChargeBalancing_2{n,3}
        reactionName = split(massChargeBalancing_2{n,1}, '_');
        if strmatch('Bacteroid', reactionName{1,1})
            metaboliteToAdd = 'Bacteroid_MNXM1[c]';
        else
            metaboliteToAdd = strcat(reactionName{1,1}, '_MNXM1[', reactionName{end,1}, ']');
        end
        reactionID = findRxnIDs(model, massChargeBalancing_2{n,1});
        metaboliteID = findMetIDs(model, metaboliteToAdd);
        model.S(metaboliteID, reactionID) = model.S(metaboliteID, reactionID) - massChargeBalancing_2{n,3};
    end
end

%% Fix EXCT reactions, TRN, TNR reactions

% Fix EXCT reactions
for n = 1:length(model.rxns)
    if strmatch('EXCT_', model.rxns{n})
        if model.S(findMetIDs(model, 'Nodule_MNXM3[c]'), n) == -0.25
            model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), n) - 0.25;
        end
    end
end

% Fix TRN and TNR reactions
for n = 1:length(model.rxns)
    if strmatch('TRN_WATER', model.rxns{n})
    elseif strmatch('TRN_WATER', model.rxns{n})
    elseif strmatch('TRN_', model.rxns{n})
        model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), n) - 0.25;
        model.S(findMetIDs(model, 'Bacteroid_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Bacteroid_MNXM2[c]'), n) - 0.25;
    elseif strmatch('TNR_', model.rxns{n})
        model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), n) - 0.25;
        model.S(findMetIDs(model, 'Bacteroid_MNXM2[c]'), n) = model.S(findMetIDs(model, 'Bacteroid_MNXM2[c]'), n) - 0.25;
    end
end

%% Fix other

model.metFormulas{findMetIDs(model, 'Leave_MNXM4897[c]')} = 'C6H10O5';
model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-GDP-FORMING-RXN_c')) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-GDP-FORMING-RXN_c')) + 1;
model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-UDP-FORMING-RXN_c')) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Leave_CELLULOSE-SYNTHASE-UDP-FORMING-RXN_c')) + 1;
model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Bacteroid_PHB_DEPOLY_c')) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), findRxnIDs(model, 'Bacteroid_PHB_DEPOLY_c')) - 1;
model.S(findMetIDs(model, 'Leave_MNXM1[c]'), findRxnIDs(model, 'Bacteroid_PHB_DEPOLY_c')) = model.S(findMetIDs(model, 'Leave_MNXM1[c]'), findRxnIDs(model, 'Bacteroid_PHB_DEPOLY_c')) + 1;



%% Fix transport reactions

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end

% Fix transport reactions
for n = 76:length(massChargeBalancing_2)
    reactionName = split(massChargeBalancing_2{n,1}, '_');
    pos = findRxnIDs(model, massChargeBalancing_2{n,1});
    if ~isempty(cell2mat(strfind(reactionName(end), 'c')))
        if strmatch('Leave', reactionName{1})
            model.S(findMetIDs(model, 'Leave_MNXM2[c]'), pos) = model.S(findMetIDs(model, 'Leave_MNXM2[c]'), pos) - 0.25;
        elseif strmatch('Root', reactionName{1})
            model.S(findMetIDs(model, 'Root_MNXM2[c]'), pos) = model.S(findMetIDs(model, 'Root_MNXM2[c]'), pos) - 0.25;
        elseif strmatch('Nodule', reactionName{1})
            model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), pos) = model.S(findMetIDs(model, 'Nodule_MNXM2[c]'), pos) - 0.25;
        end
    end
end

% Fix other
model.metFormulas{findMetIDs(model, 'Leave_MNXM794[c]')} = 'C18H30O2';
model.metCharges(findMetIDs(model, 'Leave_MNXM794[c]')) = 0;

% Get list of unbalanced reactions
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
massChargeBalancing = table2cell(table(model.rxns, imBalancedMass, imBalancedCharge));
massChargeBalancing_2 = {};
for n = 1:length(massChargeBalancing)
    if isempty(massChargeBalancing{n,2}) && massChargeBalancing{n,3} == 0
    else
        massChargeBalancing_2 = vertcat(massChargeBalancing_2, massChargeBalancing(n,:));
    end
end
