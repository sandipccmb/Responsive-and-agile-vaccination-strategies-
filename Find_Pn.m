% First run "Setup_model.m" by choosing appropriate values for p.c1 and p.c3 run this code and save the data as described in line 73-75.  

clear all; load Model_setup_India

% --- Make samples for each parameter -------------------------------------

nsam = 250;
xs = repmat(prm.bounds(1,:),nsam,1) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,size(prm.bounds,2));

p.seropos = 0.0*ones(1,3);  % Sero-positivity in India before 1st epidemic

opts = odeset('NonNegative',[1:i.nstates], 'Refine', 64, 'AbsTol', 1e-10, 'RelTol', 1e-10);

R0  = 2.5; %2.0; 
tf  = 600; % Simulation duration

mk = round(nsam/20);

% First find the baseline scenario
    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
                 
        p1 = p; r1 = r;
        p1.vacc1  = [0, 0, 0];           % Proportion vaccinated in general population
        p1.vacc2  = [0, 0, 0];           % Proportion vaccinated among the risk group
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
        [t,soln0] = ode15s(geq, [0:1:tf], r1.init, opts);
        
        inc(:,ii) = sum(diff(soln0(:,i.aux.inc),1),2); % Daily incidence
        % --- Record cumulative deaths (total)
        mor0(:,ii) = sum(soln0(:,i.aux.mort),2); 
         
                                  
    end
    fprintf('\n');
    
% Use function 'linespecer':It basically provides colour selections that are suitable for publication.
cols = linspecer(9); 

figure; fs = 14; lw = 1.5;transparency1=0.1;transparency2=0.1; 
xl0= 350;

% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                
y = permute(inc_prc,[2,1,3]);
h1=plot(y(2,:),'Color',cols(2,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:), y(1,:),cols(2,:), 'None', 1, 0.3); hold on;
ylabel({'Daily symptomatic cases per 100,000 population'}); %
xlim([0 xl0]);

%peak value with indices
[val, ind] = max(y(2,1:350));
peak_pre = y(:,ind)';
TPR = 0.1;             % TPR = 10% during the peak epidemic of the 1st wave in India

% Pn = Prevalence of non-COVID symptoms (at peak);
% Pc = Actual Prevalence with COVID symptoms (model estimated peak_pre)
% TPR (at peak) = Pc/(Pc+Pn)  
Pn = peak_pre*(1-TPR)./TPR

 
% Estimated value of Pn (per 100,000 population)
%Pn = [1751.87, 5298.58, 8512.74];  % R0 = 2 per thousand population and TPR = 0.1  
%Pn = [2615.42, 8119.59, 13420.39]; % R0 = 2.5 per thousand population and TPR = 0.1  

 
