model=readCbModel('nodulatedPlant.mat')
changeCobraSolver ('ibm_cplex');
[hans]=optimizeCbModel(model)
model = changeRxnBounds(model,'Nodule_R_CO2_tx', -1000000, 'l');
[smit]=optimizeCbModel(model)
model = changeRxnBounds(model,'TRN_CO2_tx', 0, 'b');
[smit]=optimizeCbModel(model)

surfNet(model,'Nodule_R_CO2_tx');
smit.v(3938)
model = changeRxnBounds(model,'Nodule_R_MNXM1026_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_MNXM1026_tx', -1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_MNXM1026_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_MNXM111_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_MNXM111_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_H2O_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_H2O_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_Mg_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_Mg_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_O2_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_O2_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_PROTON_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_PROTON_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_Pi_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_Pi_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_SO4_tx', 1000000, 'u');
model = changeRxnBounds(model,'Nodule_R_SO4_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_N2_tx', -1000000, 'l');
model = changeRxnBounds(model,'Nodule_R_N2_tx', 1000000, 'u');
[mit]=optimizeCbModel(model)
frodo = cell2table(horzcat(model.rxns, num2cell(mit.v)))
writetable(frodo, 'frodo.txt', 'Delimiter', '\t')
 % ureide synthesis reactions
   rxnFormula5='Nodule_CPD-15318[c] <=>  Nodule_RIBOSE-5P[c]' ;
 model=addReaction(model,'Nodule_R_RIOBOSE-5P_c',rxnFormula5, [], 0, -1000000, 1000000, 0);
 
   
  
 
 