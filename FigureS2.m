% Use function 'linespecer':It basically provides colour selections that are suitable for publication.
cols = linspecer(9); 


% Plot of Cases averted vs TPR
figure; fs = 14; fs1 = 16; lw = 1.5;transparency1=0.1;transparency2=0.1; 

subplot(1,2,1)  % Figure.S1 (A)
load simulate2_sev.mat
h1=plot(th_TPR*100,x1(2,:),'Color',cols(3,:),'linewidth', lw); hold on;
jbfill(th_TPR*100, x1(3,:), x1(1,:),cols(3,:), 'None', 1, 0.3); hold on;
xlabel({'Interventions are triggered at TPR(%)'},'fontsize',fs1); 
ylabel({'Percent deaths averted'},'fontsize',fs1);
xlim([0.5 3]); % 15.66

load simulate2p5_sev.mat
h2=plot(th_TPR*100,x1(2,:),'Color',cols(4,:),'linewidth', lw); hold on;
jbfill(th_TPR*100, x1(3,:), x1(1,:),cols(4,:), 'None', 1, 0.3); hold on;
xlim([0.5 3]); % 15.66
ylim([0 70]);
legend([h1, h2],'R0 = 2.0','R0 = 2.5','fontsize',fs)
title('Symptomatic disease preventing vaccine','fontsize',fs1)

subplot(1,2,2) % Figure.S1 (B)
load simulate2_sus.mat
h1=plot(th_TPR*100,x1(2,:),'Color',cols(3,:),'linewidth', lw); hold on;
jbfill(th_TPR*100, x1(3,:), x1(1,:),cols(3,:), 'None', 1, 0.3); hold on;
xlabel({'Interventions are triggered at TPR(%)'},'fontsize',fs1); 
ylabel({'Percent deaths averted'},'fontsize',fs1);
xlim([0.5 3]); % 15.66

load simulate2p5_sus.mat
h2=plot(th_TPR*100,x1(2,:),'Color',cols(4,:),'linewidth', lw); hold on;
jbfill(th_TPR*100, x1(3,:), x1(1,:),cols(4,:), 'None', 1, 0.3); hold on;
xlim([0.5 3]); % 15.66
ylim([0 70]);
legend([h1, h2],'R0 = 2.0','R0 = 2.5','fontsize',fs)
title('Infection preventing vaccine','fontsize',fs1)

