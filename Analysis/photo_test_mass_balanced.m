function photo_test
% code to determine whether altering rate of Co2 assimilation alters the
% cost of N fixation
clear 
close ALL
changeCobraSolver ('glpk')

%model=readCbModel('GeorgeUreideMay22.mat');
%model=readCbModel('NEWcombinedModel.mat')
model=readCbModel('mass_charge_balanced_model.mat')


model = changeRxnBounds(model,'Root_CO2_tx', 0, 'u'); 
model = changeRxnBounds(model,'Leave_RXN490-3650_p',2000000,'u');
% model_ureide=model;
% model_ureide = changeRxnBounds(model_ureide, 'TNR_ASN', 0, 'b');
% model_ureide = changeRxnBounds(model_ureide, 'TNR_GLN', 0, 'b');
% model= changeRxnBounds(model, 'TNR_S-ALLANTOIN', 0, 'b');
% model = changeRxnBounds(model, 'TNR_ALLANTOATE', 0, 'b');
% sol_ureide = optimizeCbModel(model_ureide);
nay=[];noo=[];moo=[];nfix=[];
%for n=10000:10000:140000;
n=10000

    model = changeRxnBounds(model,'Leave_CO2_tx', n, 'u');     
    RIPE=optimizeCbModel(model);
           noo=[noo,n];

rgr=RIPE.f*24/1000;
CN=0;rgr=0;N=0;
nay=[nay,CN];moo=[moo,rgr];nfix=[nfix,N];
    for n=20000:5000:200000;
aa=[];
    model = changeRxnBounds(model,'Leave_CO2_tx', n, 'u');     
    RIPE=optimizeCbModel(model);

%     model_ureide = changeRxnBounds(model_ureide,'Leave_CO2_tx', n, 'u'); 
%     RIPE1=optimizeCbModel(model_ureide);

    while RIPE.f==0
        a=n+100;
        
    model = changeRxnBounds(model,'Leave_CO2_tx',a, 'u');     
    RIPE=optimizeCbModel(model);
      %  CN=0;rgr=0;
%     elseif RIPE1.f==0
%         CN1=0;
    aa=[aa,a];
end
    % carbon cost
    if isempty(aa)==0
         noo=[noo,aa(end)]  
    elseif  isempty(aa)~=0 
        noo=[noo,n];
    end 

rgr=RIPE.f*24/1000;
Nit=find(contains(model.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N=RIPE.v(Nit); % finds the flux of nitrogenase reaction
co2=find(contains(model.rxns,'Nodule_CO2_tx')); % finds the reaction number for CO2 flux out of nodule
C=RIPE.v(co2); % finds the flux of CO2 out of nodule
CN=12.*(abs(C)./abs(N.*2))./14;

% Nit1=find(contains(model_ureide.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
% N1=RIPE1.v(Nit1); % finds the flux of nitrogenase reaction
% co21=find(contains(model_ureide.rxns,'Nodule_CO2_tx')); % finds the reaction number for CO2 flux out of nodule
% C1=RIPE1.v(co21); % finds the flux of CO2 out of nodule
% CN1=12.*(abs(C1)./abs(N1.*2))./14;
nfix=[nfix,N];    
nay=[nay,CN];moo=[moo,rgr]
end 
% figure(1)         
%             plot(noo,nay,'LineWidth',4)            
%             ylabel('Carbon cost of N fixation gC/gN','FontSize',40)     
%             xlabel('CO_2 assimilation rate \mu mol/g/hr ','FontSize',40)
% hold on, drawnow
%             plot(noo,moo,'LineWidth',4)            
%             legend('amide','ureide','Location','Best');
% 
%               set(gca,'LineWidth',2,'FontSize',40)
%                       set(gcf, 'PaperUnits', 'inches'); 
%  x_width=24 ;y_width=15;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%  print('Costvphoto','-depsc','-loose');

 figure(1)         
yyaxis left 
            plot(noo,nay,'LineWidth',4)            
            ylabel('Carbon cost of N fixation gC/gN','FontSize',40)     
            xlabel('CO_2 assimilation rate \mu mol/g/hr ','FontSize',40)

   yyaxis right 
            plot(noo,moo,'LineWidth',4)            
            ylabel('RGR g/g/d','FontSize',40)

              set(gca,'LineWidth',2,'FontSize',40)
                      set(gcf, 'PaperUnits', 'inches'); 
 x_width=24 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('Costvphoto','-depsc','-loose');

 figure(2)         

            plot(noo,nfix*2,'LineWidth',4)            
            ylabel('N fixation rate \mu mol/g/hr','FontSize',40)
            xlabel('CO_2 assimilation rate \mu mol/g/hr ','FontSize',40)

              set(gca,'LineWidth',2,'FontSize',40)
                      set(gcf, 'PaperUnits', 'inches'); 
 x_width=18 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('CostvphotoNfix','-depsc','-loose');
end