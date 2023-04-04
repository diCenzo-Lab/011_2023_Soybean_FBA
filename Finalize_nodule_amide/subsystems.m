model=readCbModel('Combined_model_with_databases.mat');
rxns=model.rxns;
rxns = erase(rxns,'Root_');
rxns = erase(rxns,'Leave_');
rxns = erase(rxns,'Nodule_');
rxns = erase(rxns,'Bacteroid_');
for n=1:length(rxns)
    if contains(rxns(n),'MNXR')~=0
rxns(n) = erase(rxns(n),'_c');
rxns(n) = erase(rxns(n),'_p');
rxns(n) = erase(rxns(n),'_m');
rxns(n) = erase(rxns(n),'_x');
    else
    end
end

metanetxChem = table2cell(readtable('reac_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));
depreciated= table2cell(readtable('reac_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));
KEGGS={};
for n=1:length(rxns)

pos1=strmatch(rxns(n),metanetxChem(:,2),'exact');
pos1=strmatch(rxns(n),depreciated(:,2),'exact');
if sum(pos1) ~=0
    rxns1=metanetxChem(pos1,1);
   kegg=find(contains(rxns1,'kegg'));
   if kegg~=0
    keggrxn=rxns1(kegg(1));
KEGGS(n)=keggrxn;
   else 
       KEGGS{n}='N/A';

   end
else
KEGGS{n}='N/A';
end
end
KEGGS=transpose(KEGGS)
KEGGS = erase(KEGGS,'keggR:');
%%

for n=1:length(model.rxns)
 if  contains(KEGGS(n),'R')~=0
chr = ['https://rest.kegg.jp/get/',KEGGS{n}];
url = convertCharsToStrings(chr)

code = webread(url);
str = extractHTMLText(code);
newStr = extractAfter(str,"PATHWAY ");

newStr2 = extractBefore(newStr,"BRITE KEGG Orthology (KO)");
newStr2 = extractBefore(newStr,"ORTHOLOGY");
model.subSystems{n,1}=newStr2
 else
     model.subSystems{n,1}='';

 end
end
model.rxns=strrep(model.rxns,'Leave_','Shoot_')
model.rxnNames=strrep(model.rxnNames,'Leave_','Shoot_')
model.mets=strrep(model.mets,'Leave_','Shoot_')
model.mets=strrep(model.mets,'Leave_','Shoot_')
model.genes=strrep(model.genes,'Leave_','Shoot_')
model.metNames=strrep(model.metNames,'Leave_','Shoot_')
    model.grRules=strrep(model.grRules,'Leave_','Shoot_')


save('combined_model_with_databases_subs.mat','model')
% 
% pat = digitsPattern
% maps = extract(newStr2,pat)
% maps = strcat('rn',maps);
% sub={};
% for i=1:(length(maps))
% Str = extractAfter(newStr2,maps{i});
% if n<length(maps)
% newsy = extractBefore(Str,maps{i+1});
% else 
%     newsy=Str;
% end
% maps1=strcat(maps{i},newsy)
% model.subSystems{n,i}=maps1
% end




