% Data for Figure1(B)
% First run "Setup_model.m" by choosing appropriate values for p.c1 and p.c3   

clear all; load Model_setup_India

% --- Make samples for each parameter -------------------------------------

nsam = 250;
xs = repmat(prm.bounds(1,:),nsam,1) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,size(prm.bounds,2));

p.seropos = 0.25*ones(1,3);  % To adjust background sero-prevalence in India before 2nd wave

opts = odeset('NonNegative',[1:i.nstates], 'Refine', 64, 'AbsTol', 1e-10, 'RelTol', 1e-10);

R0  = 2.0; %2.5; % as required
% Assign the value of 'Prevalence of non-COVID symptoms' (Obtained from Find_Pn.m)
Pn = 5298.58; % (Pn = 5298.58 for R0 = 2.0);   (Pn =8119.59 for R0 = 2.5) 

tf  = 600;    % Simulation duration

mk = round(nsam/20);

% --- Specify the order of priority for vaccination -----------------------
covrg = 0.75; ndays = 1*30; vinit = 0;  % 75% coverage among all adults above 18y      

% Rate if vaccinating uniformly among all adults (coverage 75%)
r.vacc = -1/ndays*log(1-covrg)*([[0 1 1];[0 1 1]]); 


for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
                    
        p1 = p; r1 = r;
        r1.vacc = 0*r.vacc;            % Rate of vaccination
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
        [t,soln0] = ode15s(geq, [0:1:tf], r1.init, opts);
        
         % --- Record symptomatic incidence (total)
        inc0(:,ii) = sum(diff(soln0(:,i.aux.inc),1),2); % Daily incidence
        cinc0(:,ii) = sum(soln0(:,i.aux.inc),2);        % Cumulative incidence 
        % --- Record cumulative deaths (total)
        mor0(:,ii) = sum(soln0(:,i.aux.mort),2);        % Cumulative deaths
        % --- Record daily deaths 
        mort0(:,ii) = sum(diff(soln0(:,i.aux.mort),1),2); % Daily deaths
         
       % Find TPR
        TPR(:,ii) = sum(diff(soln0(:,i.aux.inc),1),2)./(sum(diff(soln0(:,i.aux.inc),1),2)+Pn); 
        
end
        TPR_prc = prctile(TPR,[2.5,50,97.5],2); 

 
th_TPR = 0.005; % (For threshold TPR = 0.5%)  

    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
                    
        p1 = p; r1 = r;
        r1.vacc = 0*r.vacc;            % Rate of vaccination
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
      
        TVS = find(TPR_prc(:,2)>th_TPR,1,'first');

        % Run the baseline till vaccination starts
        [t1,soln1] = ode15s(geq, [0:1:TVS], r1.init, opts);

        p2 = p1; r2 = r1;
        r2.vacc = r.vacc;           % Rate of vaccination
        
        % --- Perform the simulation with vaccination
        % Accelerated vaccination to cover 75% of over-18-year-olds within 1 month of exceeding a TPR of 0.5%, 
        M2 = make_model2(p2, r2, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M2, i, s, p2, r2, agg, sel, prm);
        [t2,soln2] = ode15s(geq, [TVS:1:tf], soln1(end,:), opts);
         
        soln3 = [soln1; soln2(2:end,:)]; % Solution with vaccination intervention only
        
        % Reduction of beta due to NPI (such as use of mask etc)
        p3 = p1; r3 = r1;
        r3.vacc = r.vacc;           % Rate of vaccination
        r3.beta = r.beta*(1-0.25);
        
        % --- Perform the simulation with vaccination+NPI
        M3 = make_model2(p3, r3, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M3, i, s, p3, r3, agg, sel, prm);
        [t4,soln4] = ode15s(geq, [TVS:1:TVS+30], soln1(end,:), opts);

        % Back to the normal position after 30 days(Removing NPI)
        p4 = p3; r4 = r3;
        r4.beta = r.beta;
        % --- Perform the simulation with vaccination
        M4 = make_model2(p4, r4, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M4, i, s, p4, r4, agg, sel, prm);
        [t5,soln5] = ode15s(geq, [TVS+30:1:tf], soln4(end,:), opts);

         soln6 = [soln1; soln4(2:end,:); soln5(2:end,:)];
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        % --- Record symptomatic incidence  and mortality with vaccination 
        
        inc(:,ii) = sum(diff(soln3(:,i.aux.inc),1),2);   % Daily incidence
        cinc(:,ii) = sum(soln3(:,i.aux.inc),2);          % Cumulative incidence 
      
        mort(:,ii) = sum(diff(soln3(:,i.aux.mort),1),2); % Daily deaths
        mor(:,ii) = sum(soln3(:,i.aux.mort),2);          % cumulative deaths  
         
        % Cases averted
        CA(:,ii) = sum((inc0(:,ii) - inc(:,ii)),1);
        % Percentage reduction of cumulative incidence 
        Red_cumcases(:,ii) = 100*((cinc0(end,ii)- cinc(end,ii))./cinc0(end,ii));
         % Percentage reduction of cumulative deaths 
        Red_Deaths(:,ii) = 100*((mor0(end,ii)- mor(end,ii))./mor0(end,ii));
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         % --- Record symptomatic incidence and mortality with vaccination+NPI 
         
        inc1(:,ii) = sum(diff(soln6(:,i.aux.inc),1),2);   % Daily incidence
        cinc1(:,ii) = sum(soln6(:,i.aux.inc),2);          % Cumulative incidence 
        % --- Record deaths with vaccination (total)
        mort1(:,ii) = sum(diff(soln6(:,i.aux.mort),1),2); % Daily deaths
        mor1(:,ii) = sum(soln6(:,i.aux.mort),2);          % cumulative deaths  
        
        % Cases averted
        CA1(:,ii) = sum((inc0(:,ii) - inc1(:,ii)),1);
        % Percentage reduction of cumulative cases  
        Red_cumcases1(:,ii) = 100*((cinc0(end,ii)- cinc1(end,ii))./cinc0(end,ii));
         % Percentage reduction of cumulative deaths 
        Red_Deaths1(:,ii) = 100*((mor0(end,ii)- mor1(end,ii))./mor0(end,ii));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                 
    end
    fprintf('\n');

% Save data as below:    
% save simulate_sus;   % For infection preventing vaccine
save simulate_sev;     % For disease prevanting vaccine

% Impact of responsive vaccination
RedDeaths_prc = prctile(Red_Deaths,[2.5,50,97.5],2);    % Percentage reduction of cumulative deaths
RedCumcase_prc = prctile(Red_cumcases,[2.5,50,97.5],2); % Percentage reduction of cumulative incidence 
CA_prc = prctile(CA,[2.5,50,97.5],2);                   % Total Cases averted

% Impact of responsive vaccination + NPI
RedDeaths_prc1 = prctile(Red_Deaths1,[2.5,50,97.5],2);   % Percentage reduction of cumulative deaths
RedCumcase_prc1 = prctile(Red_cumcases1,[2.5,50,97.5],2);% Percentage reduction of cumulative incidence  
CA_prc1 = prctile(CA1,[2.5,50,97.5],2);                  % Total Cases averted


% Use function 'linespecer':It basically provides colour selections that are suitable for publication.
cols = linspecer(9); 

figure; fs = 14; lw = 2.0;transparency1=0.1;transparency2=0.1; 
xl0= 350;
yl0 = 13;
%subplot(1,2,1)
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
xlim([0 xl0]);
ylim([0 yl0]);
legend([h1,h2,h3],'No intervention','Responsive vaccination','Responsive vaccination + NPI','fontsize',fs);



