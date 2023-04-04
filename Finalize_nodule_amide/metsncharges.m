%% UpdatingSoyFBA
changeCobraSolver ('ibm_cplex'); 

soy=readCbModel('soybeangeorge.mat');
brady=readCbModel('USDA110_model.mat');
model=readCbModel('NEWcombinedModel.mat')
% % Import METANETX compound database
metanetxForm = table2cell(readtable('chem_prop.txt', 'Delimiter', '\t','ReadVariableNames', false));
metanetxChem = table2cell(readtable('chem_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));
constrainedNodule=model;
%% adding metFormulas and metCharges to constrainedNodule 


uniquemets=strrep(constrainedNodule.mets,'Leave_','');
uniquemets=strrep(uniquemets,'Root_','');
uniquemets=strrep(uniquemets,'Bacteroid_','');
uniquemets=strrep(uniquemets,'Nodule_','');
uniquemets = erase(uniquemets,'[c]');   
uniquemets = erase(uniquemets,'[p]');
uniquemets = erase(uniquemets,'[x]');
uniquemets = erase(uniquemets,'[m]');
uniquemets = erase(uniquemets,'[e]');
uniquemets=unique(uniquemets);
metFormulas={};
metCharges=[];
for n=1:length(uniquemets)
    n
    met=uniquemets{n};
pos=strmatch(met,soy.mets,'exact');  
pos1=strmatch(met,brady.mets,'exact');
% soybean met labelling
if sum(pos)~=0
    metFormulas{n,1}=soy.metFormulas{pos};
metCharges(n,1)=soy.metCharges(pos);
% brady met labelling 
elseif sum(pos1)~=0
metFormulas{n,1}=brady.metFormulas{pos1};
metCharges(n,1)= brady.metCharges(pos1);
% other labelling ~ MNXM labelling 
else 
 
pos2=strmatch(met,metanetxForm(:,1),'exact');
    pos3=strmatch(met,metanetxChem(:,2),'exact');

if sum(pos2)~=0
    metFormulas{n,1}=metanetxForm{pos2(1),4};
metCharges(n,1)=str2double(metanetxForm{pos2(1),5});
elseif sum(pos3)~=0
    rop=[];
    for i=1:length(pos3)
        pos4=strmatch(metanetxChem{pos3(i),1},metanetxForm(:,3),'exact');
        rop=[rop,pos4];
          if sum(rop)~=0
        break
          
    end
    end
    if sum(rop)~=0
    metFormulas{n,1}= metanetxForm{rop(1),4};
metCharges(n,1)=str2double(metanetxForm(rop(1),5));
else
metFormulas{n,1}='N/A';
metCharges(n,1)= NaN;
    end
else 
metFormulas{n,1}='N/A';
metCharges(n,1)= NaN;
end
end 
end

% adding the formula and charge by hand 
%1
pos=strmatch('MNXM182',uniquemets,'exact');
metFormulas{pos,1}='C6H12O6';
metCharges(pos,1)=0;
%2
pos=strmatch('MNXM231',uniquemets,'exact');
 metFormulas{pos,1}='C6H13NO2';
metCharges(pos,1)=0;
%3
pos=strmatch('MNXM242',uniquemets,'exact');
metFormulas{pos,1}='C5H10O5';
metCharges(pos,1)=0;
%4
pos=strmatch('MNXM2',uniquemets,'exact');
 metFormulas{pos,1}='H2O';
metCharges(pos,1)=0;
%5
pos=strmatch('MNXM323',uniquemets,'exact');
metFormulas{pos,1}='HO3S2';
metCharges(pos,1)=-1;
%6 
pos=strmatch('MNXM370',uniquemets,'exact');
metFormulas{pos,1}='C6H14NO8P';
metCharges(pos,1)=0;
%7 mnxm390
pos=strmatch('MNXM390',uniquemets,'exact');
metFormulas{pos,1}='C6H12O6';
metCharges(pos,1)=0;
%8 mnxm42
pos=strmatch('MNXM42',uniquemets,'exact');
metFormulas{pos,1}='C4H7NO4';
metCharges(pos,1)=0;
%9 MNXM722779
pos=strmatch('MNXM722779',uniquemets,'exact');
 metFormulas{pos,1}='C7H9O7';
metCharges(pos,1)=0;
%33 phb
pos=strmatch('phb',uniquemets,'exact');
 metFormulas{pos,1}='C4H6O2';
metCharges(pos,1)=0;
% 44 3a2opp
% 45 C6H12O6
pos=strmatch('MNXM175',uniquemets,'exact');
 metFormulas{pos,1}='C6H11O9P';
metCharges(pos,1)=0;
% 46
pos=strmatch('glc_A',uniquemets,'exact');
 metFormulas{pos,1}='C6H12O6';
metCharges(pos,1)=0;
% 47 MIGHT BE +- H
pos=strmatch('MNXM5782',uniquemets,'exact');
 metFormulas{pos,1}='C7H8O7';
metCharges(pos,1)=0;
% 48 MNXM23766
pos=strmatch('MNXM23766',uniquemets,'exact');
 metFormulas{pos,1}='C8H13OSR';
metCharges(pos,1)=-2;
% MNXM463
pos=strmatch('MNXM463',uniquemets,'exact');
 metFormulas{pos,1}='C7H14N2O8P';
metCharges(pos,1)=-1;
%  MNXM388
pos=strmatch('MNXM388',uniquemets,'exact');
 metFormulas{pos,1}='C8H13N3O7P';
metCharges(pos,1)=-1;
% MNXM568
pos=strmatch('MNXM568',uniquemets,'exact');
 metFormulas{pos,1}='C8H15N3O8P';
metCharges(pos,1)=-1;
% MNXM33
pos=strmatch('MNXM33',uniquemets,'exact');
 metFormulas{pos,1}='C27H31N9O15P2';
% MNXM565
pos=strmatch('MNXM565',uniquemets,'exact');
 metFormulas{pos,1}='C8H15N3O8P';
metCharges(pos,1)=0;
% MNXM318
pos=strmatch('MNXM318',uniquemets,'exact');
 metFormulas{pos,1}='C20H24N7O6';
%Pumped-PROTON
pos=strmatch('Pumped-PROTON',uniquemets,'exact');
 metFormulas{pos,1}='H';
metCharges(pos,1)=1;
% UDP-GLUCOSE
pos=strmatch('UDP-GLUCOSE',uniquemets,'exact');
 metFormulas{pos,1}='C15H22N2O17P2';
metCharges(pos,1)=0;
% MNXM16
pos=strmatch('MNXM16',uniquemets,'exact');
 metFormulas{pos,1}='C15H23N6O5S';
metCharges(pos,1)=1;
% MNXM8962
pos=strmatch('MNXM8962',uniquemets,'exact');
 metFormulas{pos,1}='C32H57N3O9PS';
metCharges(pos,1)=-1;
% MNXM1289
pos=strmatch('MNXM1289',uniquemets,'exact');
 metFormulas{pos,1}='C6H3O7';
metCharges(pos,1)=-3;
% MNXM7733
pos=strmatch('MNXM7733',uniquemets,'exact');
 metFormulas{pos,1}='C14H27O2S';
metCharges(pos,1)=-1;
%MNXM96993  C6H9NO2S2
pos=strmatch('MNXM96993',uniquemets,'exact');
 metFormulas{pos,1}='C10H14N4O4S2';
metCharges(pos,1)=1;
% MNXM148
pos=strmatch('MNXM148',uniquemets,'exact');
 metFormulas{pos,1}='C10H12N4O4S2';
metCharges(pos,1)=0;
% MNXM91417
pos=strmatch('MNXM91417',uniquemets,'exact');
 metFormulas{pos,1}='C27H37N7O19P3S';
metCharges(pos,1)=-5;
% MNXM507
pos=strmatch('MNXM507',uniquemets,'exact');
 metFormulas{pos,1}='C9H12N3O9P';
metCharges(pos,1)=-2;
% MNXM40588
pos=strmatch('MNXM40588',uniquemets,'exact');
 metFormulas{pos,1}='C6H12O6';
metCharges(pos,1)=0;
% MNXM4637
pos=strmatch('MNXM4637',uniquemets,'exact');
 metFormulas{pos,1}='C6H12O6';
metCharges(pos,1)=0;
% MNXM11871
pos=strmatch('MNXM11871',uniquemets,'exact');
 metFormulas{pos,1}='C13H14N2O3';
metCharges(pos,1)=0;
% MNXM90088
pos=strmatch('MNXM90088',uniquemets,'exact');
 metFormulas{pos,1}='C16H23N5O16P2';
metCharges(pos,1)=-2;
%MNXM5723
pos=strmatch('MNXM5723',uniquemets,'exact');
 metFormulas{pos,1}='C12H23OS';
metCharges(pos,1)=0;
% STARCH
pos=strmatch('STARCH',uniquemets,'exact');
 metFormulas{pos,1}='C12H20O11';
metCharges(pos,1)=0;
% Root_MNXM164034
pos=strmatch('MNXM164034',uniquemets,'exact');
 metFormulas{pos,1}='C4H7O6P';
metCharges(pos,1)=-2;
% MNXM15900
pos=strmatch('MNXM15900',uniquemets,'exact');
 metFormulas{pos,1}='C5H9O8P';
% MNXM746 *unsure_on_charge*
pos=strmatch('MNXM746',uniquemets,'exact');
 metFormulas{pos,1}='C34H32FeN4O4S2';
metCharges(pos,1)=0;
% MNXM16649
pos=strmatch('MNXM16649',uniquemets,'exact');
 metFormulas{pos,1}='C33H44N12O22P3S';
metCharges(pos,1)=-2;
% MNXM164670
pos=strmatch('MNXM164670',uniquemets,'exact');
 metFormulas{pos,1}='C33H41N12O24P3';
metCharges(pos,1)=-3;
% MNXM12055
pos=strmatch('MNXM12055',uniquemets,'exact');
 metFormulas{pos,1}='C6H11O9P';
metCharges(pos,1)=-2;
% MNXM163814
pos=strmatch('MNXM163814',uniquemets,'exact');
 metFormulas{pos,1}='C16H23N5O16P2';
metCharges(pos,1)=-2;
% MNXM7127
pos=strmatch('MNXM7127',uniquemets,'exact');
 metFormulas{pos,1}='C10H19O2S';
metCharges(pos,1)=0;
% MNXM91793
pos=strmatch('MNXM91793',uniquemets,'exact');
 metFormulas{pos,1}='C4H7O2S';
metCharges(pos,1)=-1;
% MNXM10019
pos=strmatch('MNXM10019',uniquemets,'exact');
 metFormulas{pos,1}='C8H15O2S';
metCharges(pos,1)=-2;
% MNXM583
pos=strmatch('MNXM583',uniquemets,'exact');
 metFormulas{pos,1}='C40H56O4';
metCharges(pos,1)=0;
% MNXM743
pos=strmatch('MNXM743',uniquemets,'exact');
 metFormulas{pos,1}='C6H11O9P';
metCharges(pos,1)=-2;
% MNXM74 C3H9O6P
pos=strmatch('MNXM74',uniquemets,'exact');
 metFormulas{pos,1}='C3H5O6P';
% MNXM3224b
pos=strmatch('MNXM3224b',uniquemets,'exact');
 metFormulas{pos,1}='C5H10O5';
metCharges(pos,1)=0;
% Mg-PROTOPORPHYRIN
pos=strmatch('Mg-PROTOPORPHYRIN',uniquemets,'exact');
 metFormulas{pos,1}='C34H30N4O4Mg';
metCharges(pos,1)=-2;
% Mg-PROTOPORPHYRIN-MONOMETHYL-ESTER
pos=strmatch('Mg-PROTOPORPHYRIN-MONOMETHYL-ESTER',uniquemets,'exact');
 metFormulas{pos,1}='C35H35MgN4O4';
metCharges(pos,1)=1;
% Mg+2
pos=strmatch('Mg+2',uniquemets,'exact');
 metFormulas{pos,1}='Mg';
metCharges(pos,1)=2;
% MET-tRNAs
pos=strmatch('MET-tRNAs',uniquemets,'exact');
 metFormulas{pos,1}='C28H34N11O21P3';
metCharges(pos,1)=-3;
% XYLAN
pos=strmatch('XYLAN',uniquemets,'exact');
 metFormulas{pos,1}='C10H15O8';
metCharges(pos,1)=0;
% PLASTOQUINOL-1
pos=strmatch('PLASTOQUINOL-1',uniquemets,'exact');
 metFormulas{pos,1}='C13H18O2';
metCharges(pos,1)=0;
% ARACHIDOYL-COA
pos=strmatch('ARACHIDOYL-COA',uniquemets,'exact');
 metFormulas{pos,1}='C41H70N7O17P3S';
metCharges(pos,1)=-4;
% osucc
pos=strmatch('osucc',uniquemets,'exact');
 metFormulas{pos,1}='C6H3O7';
metCharges(pos,1)=-4;
% his__L__p
pos=strmatch('his__L__p',uniquemets,'exact');
 metFormulas{pos,1}='C6H9N3O2';
metCharges(pos,1)=0;
% IRON(II)
pos=strmatch('IRON(II)',uniquemets,'exact');
 metFormulas{pos,1}='Fe';
metCharges(pos,1)=2;
% IRON(II)
pos=strmatch('IRON(III)',uniquemets,'exact');
 metFormulas{pos,1}='Fe';
metCharges(pos,1)=3;
% Linoleoyl-ACPs
pos=strmatch('Linoleoyl-ACPs',uniquemets,'exact');
 metFormulas{pos,1}='C18H31OS';
metCharges(pos,1)=0;
% Linolenoyl-ACPs
pos=strmatch('Linolenoyl-ACPs',uniquemets,'exact');
 metFormulas{pos,1}='C34H57N3O9PS';
metCharges(pos,1)=0;
% Enzyme-N6-lipoyl-L-lysine
pos=strmatch('Enzyme-N6-lipoyl-L-lysine',uniquemets,'exact');
 metFormulas{pos,1}='C15H26N2O2S2';
metCharges(pos,1)=1;
% Enzyme-N6-dihydrolipoyl-L-lysine
pos=strmatch('Enzyme-N6-dihydrolipoyl-L-lysine',uniquemets,'exact');
 metFormulas{pos,1}='C14H26N2O2S2';
metCharges(pos,1)=1;
% Enzyme-N6-S-acetyldihydrolipoyl-L-lysine
pos=strmatch('Enzyme-N6-S-acetyldihydrolipoyl-L-lysine',uniquemets,'exact');
 metFormulas{pos,1}='C17H30N2O3S2';
metCharges(pos,1)=1;
% MNXM453
pos=strmatch('MNXM453',uniquemets,'exact');
 metFormulas{pos,1}='C8H13N2O9P';
metCharges(pos,1)=-2;
% MNXM722712b
pos=strmatch('MNXM722712b',uniquemets,'exact');
 metFormulas{pos,1}='C5H9O8P';
metCharges(pos,1)=-2;
% MNXM147310
pos=strmatch('MNXM147310',uniquemets,'exact');
 metFormulas{pos,1}='HS2R';
metCharges(pos,1)=0;
% MNXM147309
pos=strmatch('MNXM147309',uniquemets,'exact');
 metFormulas{pos,1}='HSR';
metCharges(pos,1)=0;
% photon
pos=strmatch('Photon',uniquemets,'exact');
metCharges(pos,1)=0; 
% MNXM89621
pos=strmatch('MNXM89621',uniquemets,'exact');
 metFormulas{pos,1}='C6H10O12P2';
% % 3a2opp
% pos=strmatch('3a2opp',uniquemets,'exact');
%  metFormulas{pos,1}='C4H8NO7P';
% metCharges(pos,1)=-2;
for n=1:length(constrainedNodule.mets)
    n
    met = constrainedNodule.mets{n};
met = erase(met,'Leave_');
met = erase(met,'Root_');
met = erase(met,'Bacteroid_');
met = erase(met,'Nodule_');
met = erase(met,'[c]');   
met = erase(met,'[p]');
met = erase(met,'[x]');
met = erase(met,'[m]');
met = erase(met,'[e]');
pos=strmatch(met,uniquemets,'exact');
constrainedNodule.metFormulas{n,1}=metFormulas{pos};
constrainedNodule.metCharges(n,1)= metCharges(pos);
end
% identifying mets left over 

metsleft={};
for n=1:length(constrainedNodule.mets)
    if isempty(constrainedNodule.metFormulas{n}) || contains(constrainedNodule.metFormulas(n), ...
            'N/A')
        metsleft=[metsleft,constrainedNodule.mets(n)]
    else 
    end
end 
metsleft=transpose(metsleft)

% identifying lefover charges
charges=constrainedNodule.mets(find(isnan(constrainedNodule.metCharges)))
metsleft=transpose(metsleft)


constrainedNodule.metFormulas=strrep(constrainedNodule.metFormulas,'*','');
model=constrainedNodule;
save('combined_model_with_formulas','model')

% 
% 
% %%
% % adding the formula and charge by hand 
% %1
% pos=strmatch('Bacteroid_MNXM182[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H12O6';
% constrainedNodule.metCharges(pos,1)=0;
% %2
% pos=strmatch('Bacteroid_MNXM231[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H13NO2';
% constrainedNodule.metCharges(pos,1)=0;
% %3
% pos=strmatch('Bacteroid_MNXM242[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C5H10O5';
% constrainedNodule.metCharges(pos,1)=0;
% %4
% pos=strmatch('Bacteroid_MNXM2[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='H2O';
% constrainedNodule.metCharges(pos,1)=0;
% %5
% pos=strmatch('Bacteroid_MNXM323[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='HO3S2';
% constrainedNodule.metCharges(pos,1)=-1;
% %6 
% pos=strmatch('Bacteroid_MNXM370[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H14NO8P';
% constrainedNodule.metCharges(pos,1)=0;
% %7 mnxm390
% pos=strmatch('Bacteroid_MNXM390[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H12O6';
% constrainedNodule.metCharges(pos,1)=0;
% %8 mnxm42
% pos=strmatch('Bacteroid_MNXM42[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C4H7NO4';
% constrainedNodule.metCharges(pos,1)=0;
% %9 MNXM722779
% pos=strmatch('Bacteroid_MNXM722779[e]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C7H9O7';
% constrainedNodule.metCharges(pos,1)=0;
% %10 cpdFe4S4
% %11 hcit
% % pos=strmatch('Bacteroid_hcit[c]',constrainedNodule.mets,'exact');
% %   constrainedNodule.metFormulas{pos,1}='C7H7O7';
% % constrainedNodule.metCharges(pos,1)=-3;
% %12 mobd
% % pos=strmatch('Bacteroid_mobd[c]',constrainedNodule.mets,'exact');
% %   constrainedNodule.metFormulas{pos,1}='H2MoO4';
% % constrainedNodule.metCharges(pos,1)=0;
% %13 symCoF
% %14 photon
% %15 MNXM182
% pos=strmatch('Nodule_MNXM182[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H12O6';
% constrainedNodule.metCharges(pos,1)=0;
% %16 _MNXM231
% pos=strmatch('Nodule_MNXM231[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H13NO2';
% constrainedNodule.metCharges(pos,1)=0;
% %17 MNXM242
% pos=strmatch('Nodule_MNXM242[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C5H10O5';
% constrainedNodule.metCharges(pos,1)=0;
% %18-21 MNXM2
% pos=strmatch('Nodule_MNXM2[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='H2O';
% constrainedNodule.metCharges(pos,1)=0;
% pos=strmatch('Nodule_MNXM2[m]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='H2O';
% constrainedNodule.metCharges(pos,1)=0;
% pos=strmatch('Nodule_MNXM2[p]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='H2O';
% constrainedNodule.metCharges(pos,1)=0;
% pos=strmatch('Nodule_MNXM2[x]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='H2O';
% constrainedNodule.metCharges(pos,1)=0;
% % 22 MNXM323
% pos=strmatch('Nodule_MNXM323[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='HO3S2';
% constrainedNodule.metCharges(pos,1)=-1;
% %23 MNXM370
% pos=strmatch('Nodule_MNXM370[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H14NO8P';
% constrainedNodule.metCharges(pos,1)=0;
% %24 MNXM390
% pos=strmatch('Nodule_MNXM390[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C6H12O6';
% constrainedNodule.metCharges(n,1)=0;
% %25-27 MNXM42
% pos=strmatch('Nodule_MNXM42[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C4H7NO4';
% constrainedNodule.metCharges(pos,1)=0;
% pos=strmatch('Nodule_MNXM42[m]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C4H7NO4';
% constrainedNodule.metCharges(pos,1)=0;
% pos=strmatch('Nodule_MNXM42[p]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C4H7NO4';
% constrainedNodule.metCharges(pos,1)=0;
% %28 MNXM722779
% pos=strmatch('Nodule_MNXM722779[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C7H9O7';
% constrainedNodule.metCharges(pos,1)=0;
% %29 MNXM746
% %33 phb C12H20O6
% pos=strmatch('Bacteroid_phb[c]',constrainedNodule.mets,'exact');
%   constrainedNodule.metFormulas{pos,1}='C4H6O2';
% constrainedNodule.metCharges(pos,1)=0;
% % %% generating multiple database IDs
% % 
% % metanetx=model.mets(find(contains(model.mets,'MNXM')));
% % 
% % Other=setdiff(model.mets,metanetx);
% % load('mappedCompounds.mat')
% % model.metIDs={};leftovers=[];
% % for n=1:length(model.mets)
% %  n   
% % boo=contains(model.mets{n},'MNXM')
% %     if boo~=0
% % met = model.mets{n};
% % met = erase(met,'Leave_');
% % met = erase(met,'Root_');
% % met = erase(met,'Bacteroid_');
% % met = erase(met,'Nodule_');
% % met = erase(met,'[c]');
% % met = erase(met,'[p]');
% % met = erase(met,'[x]');
% % met = erase(met,'[m]');
% % met = erase(met,'[e]');
% % 
% % 
% % pos=strmatch(met,metacyc_trunca2(:,1),'exact');  
% % pos1=strmatch(met,metanetxChem(:,2),'exact');
% % 
% % if sum(pos) ~=0
% % keggMet = metacyc_trunca2{pos(1),2};
% % fields = splitString(keggMet, '\;');
% % legg=fields{1};
% % leggfields = splitString(legg, '\_');
% % model.metIDs{n,1}=leggfields{1};
% % elseif sum(pos1) ~=0
% % poos=find(contains(metanetxChem(pos1,1),'kegg'))
% %     if sum(pos1) > 1 && sum(poos)~=0
% %    % here 
% %    keggMet1 = metanetxChem{pos1(poos),1};
% %    keggMet1 = erase(keggMet1,'kegg:')
% % model.metIDs{n,1}=keggMet1;
% %     else
% % keggMet1 = metanetxChem{pos1(1),1};
% % model.metIDs{n,1}=keggMet1;
% %     end
% % else
% %     model.metIDs{n,1}='N/A';
% %     leftovers=[leftovers,model.mets(n)];
% % 
% % end
% %     else 
% % model.metIDs{n,1}='N/A';
% %     leftovers=[leftovers,model.mets(n)];
% % 
% % end
% % end 
% % 
% % 
% % 
% % % Extract just the METACYC compound names
% % metaMets = {};
% % x = 0;
% % for n = 1:length(metanetxChem)
% %     if strmatch('metacyc', metanetxChem{n,1})
% %         x = x + 1;
% %         metaMets(x,:) = metanetxChem(n,:);
% %     end
% % end
% % % 
% % % % Import METANETX compound database
% % % metanetxChem = table2cell(readtable('chem_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));
% % % 
% % % 
% % % 
% % % %% Change metabolite names
% % % 
% % % % Extract just the BIGG compound names
% % % seedMets = {};
% % % x = 0;
% % % for n = 1:length(metanetxChem)
% % %     if strmatch('bigg', metanetxChem{n,1})
% % %         x = x + 1;
% % %         seedMets(x,:) = metanetxChem(n,:);
% % %     end
% % % end
% % % 
% % % 
% % % % Extract just the METACYC compound names
% % % mnxmMets = {};
% % % x = 0;
% % % for n = 1:length(metanetxChem)
% % %     if strmatch('MNXM', metanetxChem{n,1})
% % %         x = x + 1;
% % %         mnxmMets(x,:) = metanetxChem(n,:);
% % %     end
% % % end
% % % 
% % % % Extract model compound names
% % % cpdNames = cell(length(model.mets), 3);
% % % for n = 1:length(model.mets)
% % %     if strmatch('Leave_', model.mets{n})
% % %         cpdNames{n,1} = 'Leave_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'Leave_', '');
% % %     elseif strmatch('Root_', model.mets{n})
% % %         cpdNames{n,1} = 'Root_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'Root_', '');
% % %     elseif strmatch('Nodule_', model.mets{n})
% % %         cpdNames{n,1} = 'Nodule_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'Nodule_', '');
% % %     elseif strmatch('NoduleI_', model.mets{n})
% % %         cpdNames{n,1} = 'NoduleI_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'NoduleI_', '');
% % %     elseif strmatch('NoduleIId_', model.mets{n})
% % %         cpdNames{n,1} = 'NoduleIId_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIId_', '');
% % %     elseif strmatch('NoduleIIp_', model.mets{n})
% % %         cpdNames{n,1} = 'NoduleIIp_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIIp_', '');
% % %     elseif strmatch('NoduleIZ_', model.mets{n})
% % %         cpdNames{n,1} = 'NoduleIZ_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIZ_', '');
% % %     elseif strmatch('NoduleIII_', model.mets{n})
% % %         cpdNames{n,1} = 'NoduleIII_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIII_', '');
% % %     elseif strmatch('BacteroidIId_', model.mets{n})
% % %         cpdNames{n,1} = 'BacteroidIId_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIId_', '');
% % %     elseif strmatch('BacteroidIIp_', model.mets{n})
% % %         cpdNames{n,1} = 'BacteroidIIp_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIIp_', '');
% % %     elseif strmatch('BacteroidIZ_', model.mets{n})
% % %         cpdNames{n,1} = 'BacteroidIZ_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIZ_', '');
% % %     elseif strmatch('Bacteroid_', model.mets{n})
% % %         cpdNames{n,1} = 'Bacteroid_';
% % %         cpdNames{n,2} = strrep(model.mets{n}, 'Bacteroid_', '');
% % %     else
% % %         cpdNames{n,2} = model.mets{n};
% % %     end
% % %     splitString = strsplit(cpdNames{n,2}, '[');
% % %     if length(splitString) == 2
% % %         cpdNames{n,2} = splitString{1};
% % %         cpdNames{n,3} = splitString{2};
% % %     end
% % % end
% % % 
% % % % Change compound names, where possible
% % % for n = 1:length(model.mets)
% % %     if strmatch('Bacteroid', cpdNames{n,1})
% % %         cpd = ['bigg:' cpdNames{n,2}];
% % %         pos = strmatch(cpd, seedMets(:,1), 'exact');
% % %         if ~isempty(pos)
% % %             cpdNames{n,2} = seedMets{pos,2};
% % %         end
% % %     else
% % %         cpd = ['metacyc:' cpdNames{n,2}];
% % %         pos = strmatch(cpd, metaMets(:,1), 'exact');
% % %         if ~isempty(pos)
% % %             cpdNames{n,2} = metaMets{pos,2};
% % %         end
% % %     end
% % % end
% % % for n = 1:length(model.mets)
% % %     cpd = strcat(cpdNames{n,1}, cpdNames{n,2});
% % %     if ~isempty(cpdNames{n,3})
% % %         cpdTemp = strcat(cpd, '_', reactionNames{n,3});
% % %         if strmatch(cpdTemp, model.mets, 'exact')
% % %             model.mets{n,1} = strcat(cpd, 'b', '[', cpdNames{n,3});
% % %         else
% % %             model.mets{n,1} = strcat(cpd, '[', cpdNames{n,3});
% % %         end        
% % %     else
% % %         if strmatch(cpdTemp, model.mets, 'exact')
% % %             model.mets{n,1} = strcat(cpd, 'b');
% % %         else
% % %             model.mets{n,1} = cpd;
% % %         end        
% % %     end
% % % end
