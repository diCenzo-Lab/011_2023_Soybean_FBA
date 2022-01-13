%% Load models

% Load USDA110 model
m = readCbModel('iYY1101.xml');

% Load Rm1021 model
sino = readCbModel('iGD1348.xml');
for n = 1:length(sino.mets)
    sino.mets{n} = strrep(sino.mets{n}, '[e0]', '[e]');
    sino.mets{n} = strrep(sino.mets{n}, '[c0]', '[c]');
end

%% Modify iYY1101

% Modify nitrogenase reaction to produce a specific NH3
m = addMetabolite(m, 'fixed[c]', 'metName', 'fixed ammonia', 'metFormula', 'H4N', 'charge', 1);
m = changeRxnMets(m, {'nh3[c]'}, {'fixed[c]'}, 'NIT');

% Add a nitrogen conversion reaction
m = addMetabolite(m, 'symCoF[c]', 'metName', 'symbiosis cofactors', 'metFormula', '', 'charge', 0);
m = addReaction(m, 'convFixed', 'reactionFormula', '1.0 fixed[c] + 0.001 symCoF[c] -> 1.0 nh3[c]');

% Add a symbiosis cofactor production reaction based on the version in
% ViNE. Unlike in ViNE, does not contain Mg or adenosylcobalamin.
m = addReaction(m, 'rxnSymCoF', 'reactionFormula', ...
    'h[c] + atp[c] + pheme[c] + gthox[c] + pydxn[c] + mobd[c] + cpdFe4S4[c] + hcit[c] -> symCoF[c]');

% Add missing reactions so that the symbiosis cofactor reaction works
m = addReaction(m, 'MOBDabc', 'reactionFormula', 'atp[c] + h2o[c] + mobd[e] -> adp[c] + h[c] + pi[c] + mobd[c]', ...
    'reversible', false, 'lowerBound', 0);
m = addReaction(m, 'EX_MOBD', 'reactionFormula', '1.0 mobd[e] <=> ', 'lowerBound', -1000);
m = addReaction(m, 'rxnNifU', 'reactionFormula', 'fe2[c] + sufsesh[c] -> cpdFe4S4[c] + sufse[c]', ...
    'reversible', false, 'lowerBound', 0);
m = addReaction(m, 'rxnGFTPcpd00919', 'reactionFormula', 'atp[c] + h2o[c] + hcit[e] -> adp[c] + h[c] + pi[c] + hcit[c]', ...
    'reversible', false, 'lowerBound', 0);
m = addReaction(m, 'EX_HCIT', 'reactionFormula', '1.0 hcit[e] <=> ', 'lowerBound', -1000);

% Update PHB synthesis
m = changeRxnBounds(m, 'HACD1', -1000, 'l');
m = addMetabolite(m, 'phb[c]', 'metName', 'Poly-hydroxybutyrate');
m = addGenes(m, {'blr3732'; 'bll4360'});
m = addReaction(m, 'PHB_syn', 'reactionFormula', '3hbcoa[c] -> coa[c] + phb[c]', ...
    'reversible', false, 'lowerBound', 0, 'geneRule', 'blr3732 or bll4360');
m = addReaction(m, 'EX_PHB', 'reactionFormula', '1.0 phb[c] <=> ', 'reversible', false, 'lowerBound', 0);

%% Save and clear

USDA110 = m;
save('USDA110_model.mat', 'USDA110');
clear

