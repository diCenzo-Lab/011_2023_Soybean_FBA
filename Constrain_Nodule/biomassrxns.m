
changeCobraSolver ('ibm_cplex');
model=nodulatedPlant;
%model=readCbModel('WorkingSoy.mat');
%% add reactions missing in order for galactose pathway to be complete 
rxnFormula800='Leave_UDP-GLUCOSE[c] <=> Leave_CPD-14553[c]';
model=addReaction(model,'Leave_R_MNXR152255',rxnFormula800, [], 0, -1000000, 1000000, 0); 
rxnFormula200='Leave_PPI[c] + Leave_CPD-14553[c] <=> Leave_UTP[c] + 3 Leave_PROTON[c] + Leave_GALACTOSE-1P[c]';
model=addReaction(model,'Leave_R_MNXR144933',rxnFormula200, [], 0, -1000000, 1000000, 0);  
%% Hand written biomass reactions
%updated leave biomass reaction 
%model=removeRxns(model,{'Root_R_BiomassRoot','Leave_R_BiomassShoot'};
%rxnFormula7='62.6329016 Leave_LEU[c] + 174.768323 Leave_L-ALPHA-ALANINE[c] + 72.569145 Leave_5-OXOPROLINE[c] + 3.08007786 Leave_ATP[c] + 353.554109 Leave_BETA-D-FRUCTOSE[c] + 154.416579 Leave_FUM[c] + 18.3569438 Leave_MET[c] + 69.4414228 Leave_4-AMINO-BUTYRATE[c] + 95.8229334 Leave_CYS[c] + 141.655201 Leave_RIBOSE-5P[c] + 20.1672457 Leave_PUTRESCINE[c] + 57.1713752 Leave_ALPHA-D-GALACTOSE[c] + 1.32402344 Leave_LYS[c] + 36.7705543 Leave_TYR[c] + 1260.91011 Leave_CELLULOSE[c] + 3.69332839 Leave_CHLOROPHYLL-A[p] + 1.98350167 Leave_CHLOROPHYLL-B[p] + 51.5463272 Leave_CIT[c] + 299.40177 Leave_CONIFERYL-ALCOHOL[c] + 28.28150 Leave_COUMARYL-ALCOHOL[c] + 285.799202 Leave_CPD-12513[c] + 35.3059189 Leave_UDP-GLUCOSE[c] + 55.7912626 Leave_CPD-13118[c] + 3.08007786 Leave_CTP[c] + 38.509125 Leave_Charged-ALA-tRNAs[c] + 1.83792352 Leave_Charged-ASN-tRNAs[c] + 35.4641168 Leave_Charged-ASP-tRNAs[c] + 3.70785787 Leave_Charged-CYS-tRNAs[c] + 418.152167 Leave_Charged-GLN-tRNAs[c] + 35.8172645 Leave_Charged-GLT-tRNAs[c] + 37.6753713 Leave_Charged-GLY-tRNAs[c] + 17.8249483 Leave_Charged-ILE-tRNAs[c] + 31.0453613 Leave_Charged-LEU-tRNAs[c] + 230.640704 Leave_Charged-LYS-tRNAs[c] + 17.3139025 Leave_Charged-MET-tRNAs[c] + 22.5154572 Leave_Charged-PHE-tRNAs[c] + 24.5940026 Leave_Charged-PRO-tRNAs[c] + 29.7471671 Leave_Charged-SER-tRNAs[c] + 38.7447422 Leave_Charged-THR-tRNAs[c] + 0.30011943 Leave_Charged-TRP-tRNAs[c] + 27.3168168 Leave_Charged-TYR-tRNAs[c] + 27.1339165 Leave_Charged-VAL-tRNAs[c] + 7.67867026 Leave_DATP[c] + 4.70628177 Leave_DCTP[c] + 4.70628177 Leave_DGTP[c] + 101.17504 Leave_GDP-L-GALACTOSE[c] + 127.839729 Leave_GDP-MANNOSE[c] + 20.3411237 Leave_GLC[c] + 129.417411 Leave_GLT[c] +12.6383967 Leave_GLY[c] + 3.08007786 Leave_GTP[c] + 11.3697239 Leave_ILE[c] + 94.7413012 Leave_L-ASPARTATE[c] + 40.24253493 Leave_LINOLEIC_ACID[c] + 95.9866249 Leave_MAL[c] + 41.3271333 Leave_MYO-INOSITOL[c] + 62.92592514 Leave_OLEATE-CPD[c] + 25.85735708 Leave_PALMITATE[c] + 10.8377686 Leave_PHE[c] + 25.9583867 Leave_PRO[c] + 55.6760858 Leave_SER[c] + 38.82740 Leave_SINAPYL-ALCOHOL[c] + 69.3827709 Leave_STARCH[p] + 44.93943796 Leave_STEARIC_ACID[c] + 3.91278817 Leave_SUCROSE[c] + 6.21652439 Leave_SUC[c] + 42.1031949 Leave_THR[c] + 7.67867026 Leave_TTP[c] + 748.282201 Leave_UDP-D-XYLOSE[c] + 30.4398718 Leave_UDP-L-RHAMNOSE[c] + 3.08007786 Leave_UTP[c] + 14.935839 Leave_VAL[c] + 2460.064405 Leave_WATER[c] -> 38.509125 Leave_ALA-tRNAs[c] + 1.83792352 Leave_ASN-tRNAs[c] + 35.4641168 Leave_ASP-tRNAs[c] + 1.15345 Leave_BiomassShoot[c] + 3.70785787 Leave_CYS-tRNAs[c] + 284.8060315 Leave_GDP[c] + 418.152167 Leave_GLN-tRNAs[c] + 35.8172645 Leave_GLT-tRNAs[c] + 37.6753713 Leave_GLY-tRNAs[c] + 17.8249483 Leave_ILE-tRNAs[c] + 31.0453613 Leave_LEU-tRNAs[c] + 230.640704 Leave_LYS-tRNAs[c] + 17.3139025 Leave_MET-tRNAs[c] + 22.5154572 Leave_PHE-tRNAs[c] + 37.0902 Leave_PPI[c] + 24.5940026 Leave_PRO-tRNAs[c] + 283.1518898 Leave_PROTON[c] + 29.7471671 Leave_SER-tRNAs[c] + 38.7447422 Leave_THR-tRNAs[c] + 0.30011943 Leave_TRP-tRNAs[c] + 27.3168168 Leave_TYR-tRNAs[c] + 1099.827194 Leave_UDP[c] + 27.1339165 Leave_VAL-tRNAs[c]';
%model=addReaction(model,'Leave_R_BiomassShoot',rxnFormula7, [], 0, 0, 1000000, 0);  
% updated root biomass reaction
%rxnFormula999='20.3365501 Root_LEU[c] + 142.758821 Root_L-ALPHA-ALANINE[c] + 38.002186 Root_5-OXOPROLINE[c] + 1.89973151 Root_ATP[c] + 681.769778 Root_BETA-D-FRUCTOSE[c] + 75.8243336 Root_FUM[c] + 3.45465531 Root_MET[c] + 51.9525977 Root_4-AMINO-BUTYRATE[c] + 41.7733828 Root_CYS[c] + 152.149446 Root_RIBOSE-5P[c] + 73.1291161 Root_PUTRESCINE[c] + 102.342059 Root_ALPHA-D-GALACTOSE[c] + 6.71812237 Root_LYS[c] + 21.3099275 Root_TYR[c] + 1492.15663 Root_CELLULOSE[c] + 25.7560117 Root_CIT[c] + 354.311017 Root_CONIFERYL-ALCOHOL[c] + 33.4682248 Root_COUMARYL-ALCOHOL[c] + 23.6161303 Root_UDP-GLUCOSE[c] + 1.89973151 Root_CTP[c] +48.69614 Root_Charged-ALA-tRNAs[c] + 0.243643 Root_Charged-ASN-tRNAs[c] + 36.51174 Root_Charged-ASP-tRNAs[c] + 1.883067 Root_Charged-CYS-tRNAs[c] + 25.49445 Root_Charged-GLN-tRNAs[c] + 28.38151 Root_Charged-GLT-tRNAs[c] +  47.90468 Root_Charged-GLY-tRNAs[c] + 26.94779 Root_Charged-ILE-tRNAs[c] + 36.25194 Root_Charged-LEU-tRNAs[c] + 120.0183 Root_Charged-LYS-tRNAs[c] + 10.92005 Root_Charged-MET-tRNAs[c] + 23.93824 Root_Charged-PHE-tRNAs[c] + 36.74606 Root_Charged-PRO-tRNAs[c] + 35.01283 Root_Charged-SER-tRNAs[c] +  29.83545 Root_Charged-THR-tRNAs[c] + 0.188491 Root_Charged-TRP-tRNAs[c] + 16.03218 Root_Charged-TYR-tRNAs[c] + 37.04172 Root_Charged-VAL-tRNAs[c] + 10.2585357 Root_DATP[c] + 6.28748963 Root_DCTP[c] + 6.28748963 Root_DGTP[c] + 224.36361 Root_GDP-L-GALACTOSE[c] + 60.3337231 Root_GDP-MANNOSE[c] + 36.5521021 Root_GLC[c] + 82.6807767 Root_GLT[c] +  6.15095261 Root_GLY[c] + 1.89973151 Root_GTP[c] + 8.06783348 Root_ILE[c] + 45.6342485 Root_L-ASPARTATE[c] + 18.200026 Root_LINOLEIC_ACID[c] + 120.05663 Root_MAL[c] +   29.7472624 Root_MYO-INOSITOL[c] + 6.21253603 Root_OLEATE-CPD[c] + 24.8147746 Root_PALMITATE[c] + 3.97599951 Root_PHE[c] + 9.98990771 Root_PRO[c] + 50.4590946 Root_SER[c] + 45.9482118 Root_SINAPYL-ALCOHOL[c] + 30.3844679 Root_STEARIC_ACID[c] + 4.43972981 Root_SUCROSE[c] + 5.95321865 Root_SUC[c] + 40.5179491 Root_THR[c] +   10.2585357 Root_TTP[c] + 1066.35667 Root_UDP-D-XYLOSE[c] + 250.408398 Root_UDP-L-RHAMNOSE[c] + 1.89973151 Root_UTP[c] + 12.4966529 Root_VAL[c] + 2227.817768 Root_WATER[c] ->   48.69614 Root_ALA-tRNAs[c] + 0.243643 Root_ASN-tRNAs[c] + 36.51174 Root_ASP-tRNAs[c] + 1.0852144 Root_BiomassRoot[c] + 1.883067 Root_CYS-tRNAs[c] + 284.6973328 Root_GDP[c] + 25.49445 Root_GLN-tRNAs[c] +  28.38151 Root_GLT-tRNAs[c] + 47.90468 Root_GLY-tRNAs[c] + 26.94779 Root_ILE-tRNAs[c] +  36.25194 Root_LEU-tRNAs[c] + 120.0183 Root_LYS-tRNAs[c] + 10.92005 Root_MET-tRNAs[c] + 23.93824 Root_PHE-tRNAs[c] +   40.69097675 Root_PPI[c] + 36.74606 Root_PRO-tRNAs[c] + 1009.88464 Root_PROTON[c] + 35.01283 Root_SER-tRNAs[c] + 29.83545 Root_THR-tRNAs[c] + 0.188491 Root_TRP-tRNAs[c] + 16.03218 Root_TYR-tRNAs[c] + 1340.381201 Root_UDP[c] +  37.04172 Root_VAL-tRNAs[c]';
%model=addReaction(model,'Root_R_BiomassRoot',rxnFormula999, [], 0, 0, 1000000, 0);  

%% importing & adding new biomass reactions 
comp=readcell('stoich.xlsx');
mets=comp(:,1);
leafmets=strcat('Leave_',mets);
rootmets=strcat('Root_',mets);
leafstoich=comp(:,2);
rootstoich=comp(:,3);
model = addReaction(model, 'Leave_R_BiomassShoot', 'metaboliteList', leafmets, ...
    'stoichCoeffList', cell2mat(leafstoich), 'reversible', false, 'lowerBound', 0, ...
    'upperBound', 1000000);
model = addReaction(model, 'Root_R_BiomassRoot', 'metaboliteList', rootmets, ...
    'stoichCoeffList', cell2mat(rootstoich), 'reversible', false, 'lowerBound', 0, ...
    'upperBound', 1000000);
%%
%[metList,stoich] = findMetsFromRxns(model,'Leave_R_BiomassShoot')
%[metList] = findMetsFromRxns(model,'Root_R_BiomassRoot')
%     model = addReaction(model,'Leave_R_BiomassShoot', 'metaboliteList', metList, 'stoichCoeffList', stoich, ...
        %    'reversible', 0, 'lowerBound', 0, 'upperBound', 1000, 'objectiveCoef', 0);



