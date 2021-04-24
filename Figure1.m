clear all

% Figure 1. Feasibility and impact of rapidly responsive vaccination strategies


[num, txt, raw] = xlsread('districts.xlsx','Delhi','O1:Q134');
% datecell =  txt(2:end,1); % Form 01/12/2020 to 12/04/2021
% Cases = num(:,1);         % Form 01/12/2020 to 12/04/2021
% TPR = 100*num(:,2);       % Form 01/12/2020 to 12/04/2021

datecell =  txt(18:118,1); % Form 17/12/2020 to 27/03/2021
Cases = num(17:117,1);     % Form 17/12/2020 to 27/03/2021
TPR = 100*num(17:117,2);   % Form 17/12/2020 to 27/03/2021

figure; fs = 14; lw = 2; %1.5;
subplot(1,2,1)
h1= plot(datetime(datecell),Cases,'linewidth',lw)
ylabel('Confirmed daily cases','fontsize',fs)
yyaxis right
h2 =plot(datetime(datecell),TPR,'linewidth', lw)
ylabel('TPR (%)','fontsize',fs)
tstart = datetime(datecell(1));
tend = datetime(datecell(101));
xlim([tstart tend]);
hold on;
h3 =plot(datetime(datecell),repmat(0.5,1,size(datetime(datecell),1)),'-.k','linewidth',lw)
legend([h1,h2],'Daily cases','TPR','fontsize',10);


load simulate_sev.mat

% Use function 'linespecer':It basically provides colour selections that are suitable for publication.
cols = linspecer(9); 

%figure; fs = 14; lw = 2.0;transparency1=0.1;transparency2=0.1; 
xl0= 350;
yl0 = 13;
subplot(1,2,2)
plot([TVS,TVS],[0,yl0],'-.k','linewidth',1.0); hold on;
plot([TVS+30,TVS+30],[0,yl0],'-.k','linewidth',1.0); hold on;
% Find percentiles for incidence
mort0_prc = prctile(mort0,[2.5,50,97.5],2);                
y = permute(mort0_prc,[2,1,3]);
h1=plot(y(2,:),'Color',cols(1,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:), y(1,:),cols(1,:), 'None', 1, 0.15); hold on;

mort_prc = prctile(mort,[2.5,50,97.5],2);                
y = permute(mort_prc,[2,1,3]);
h2=plot(y(2,:),'Color',cols(2,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:), y(1,:),cols(2,:), 'None', 1, 0.15); hold on;
%ylabel({'Symptomatic cases per thousand population'},'fontsize',fs); 

mort_prc1 = prctile(mort1,[2.5,50,97.5],2);                
y = permute(mort_prc1,[2,1,3]);
h3=plot(y(2,:),'Color',cols(3,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:), y(1,:),cols(2,:), 'None', 1, 0.15); hold on;
%ylabel({'Symptomatic cases per thousand population'},'fontsize',fs); 
ylabel({'Deaths per 100,000 population'},'fontsize',fs); 

xlabel('Days (second wave)')
xlim([50 250]);
ylim([0 yl0]);
legend([h1,h2,h3],'No intervention','Responsive vaccination','Responsive vaccination + NPI','fontsize',10);







% % Plot with full dataset for Delhi
% 
% [num, txt, raw] = xlsread('districts.xlsx','Delhi','J1:M353');
% datecell =  txt(8:end,1); % Form 02/05/2020 to 12/04/2021
% Cases = num(6:end,1);         % Form 26/04/2020 to 12/04/2021
% TPR = 100*num(6:end,3);       % Form 26/04/2020 to 12/04/2021
% 
% figure; fs = 14; lw = 1.5;
% h1= plot(datetime(datecell),Cases,'linewidth',lw)
% ylabel('Confirmed daily cases','fontsize',fs)
% yyaxis right
% h2 =plot(datetime(datecell),TPR,'linewidth', lw)
% ylabel('TPR (%)','fontsize',fs)
% tstart = datetime(datecell(1));
% tend = datetime(datecell(end));
% xlim([tstart tend]);
% hold on;
% h3 =plot(datetime(datecell),repmat(0.5,1,size(datetime(datecell),1)),'-.k','linewidth',lw)
% legend([h1,h2],'Daily cases','TPR','fontsize',fs);
% 
