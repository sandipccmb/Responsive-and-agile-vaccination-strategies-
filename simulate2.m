% (Data generation for Figure S1)
% First run "Setup_model.m" by choosing appropriate values for p.c1 and p.c3   

clear all; load Model_setup_India

% --- Make samples for each parameter -------------------------------------
nsam = 250;
xs = repmat(prm.bounds(1,:),nsam,1) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,size(prm.bounds,2));

p.seropos = 0.25*ones(1,3);      % To adjust background sero-prevalence in India

opts = odeset('NonNegative',[1:i.nstates], 'Refine', 64, 'AbsTol', 1e-10, 'RelTol', 1e-10);

R0  = 2.0; %2.5
% Assign the value of 'Prevalence of non-COVID symptoms' (Obtained from Find_Pn.m)
Pn = 5298.58; % (Pn = 5298.58 for R0 = 2.0);   (Pn =8119.59 for R0 = 2.5) 

tf  = 600;    % Simulation duration

mk = round(nsam/20);

% --- Specify the order of priority for vaccination -----------------------
covrg = 0.75; ndays = 1*30; vinit = 0;

% Rate if vaccinating uniformly among all adults
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
        inc0(:,ii,1) = sum(diff(soln0(:,i.aux.inc),1),2); % Daily incidence
        cinc0(:,ii,1) = sum(soln0(:,i.aux.inc),2);        % Cumulative incidence 
        % --- Record cumulative deaths (total)
        mor0(:,ii,1) = sum(soln0(:,i.aux.mort),2); 
         
        TPR(:,ii) = sum(diff(soln0(:,i.aux.inc),1),2)./(sum(diff(soln0(:,i.aux.inc),1),2)+Pn); 
        
end
        TPR_prc = prctile(TPR,[2.5,50,97.5],2); 

% Threshold TPR (here varied from 0.5% to 3% in 5 discreate steps)
        th_TPR = linspace(0.005,0.03,5);

for ip = 1:length(th_TPR)
    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
                    
        p1 = p; r1 = r;
        r1.vacc = 0*r.vacc;            % Rate of vaccination
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
       
        TVS = find(TPR_prc(:,2)>th_TPR(ip),1,'first');

        % Run the baseline till vaccination starts
        [t1,soln1] = ode15s(geq, [0:1:TVS], r1.init, opts);

        
        p2 = p1; r2 = r1;
        r2.vacc = r.vacc;           % Rate of vaccination
        
        % --- Perform the simulation with vaccination
        M2 = make_model2(p2, r2, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M2, i, s, p2, r2, agg, sel, prm);
        [t2,soln2] = ode15s(geq, [TVS:1:tf], soln1(end,:), opts);
         
        soln3 = [soln1; soln2(2:end,:)];
             
        % Reduction of beta due to NPI (such as use of mask)
        p3 = p1; r3 = r1;
        r3.vacc = r.vacc;              % Rate of vaccination
        r3.beta = r.beta*(1-0.25);
        
        % --- Perform the simulation with vaccination
        M3 = make_model2(p3, r3, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M3, i, s, p3, r3, agg, sel, prm);
        [t4,soln4] = ode15s(geq, [TVS:1:TVS+30], soln1(end,:), opts);

        % Back to the normal position
        p4 = p3; r4 = r3;
        r4.beta = r.beta;
        % --- Perform the simulation with vaccination
        M4 = make_model2(p4, r4, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M4, i, s, p4, r4, agg, sel, prm);
        [t5,soln5] = ode15s(geq, [TVS+30:1:tf], soln4(end,:), opts);

         soln6 = [soln1; soln4(2:end,:); soln5(2:end,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        % --- Record symptomatic incidence  and mortality with vaccination 

        inc(:,ii,ip) = sum(diff(soln3(:,i.aux.inc),1),2); % Daily incidence
        cinc(:,ii,ip) = sum(soln3(:,i.aux.inc),2);        % Cumulative incidence 
        mor(:,ii,ip) = sum(soln3(:,i.aux.mort),2);        % Cumulative deaths
        
        % Cases averted
        CA(:,ii,ip) = sum((inc0(:,ii,1) - inc(:,ii,ip)),1);
        % Percentage reduction of cumulative cases  
        Red_cumcases(:,ii,ip) = 100*((cinc0(end,ii,1)- cinc(end,ii,ip))./cinc0(end,ii,1));
        % Percentage reduction of cumulative deaths 
        Red_Deaths(:,ii,ip) = 100*((mor0(end,ii,1)- mor(end,ii,ip))./mor0(end,ii,1));
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         % --- Record symptomatic incidence and mortality with vaccination + NPI       
 
        inc1(:,ii,ip) = sum(diff(soln6(:,i.aux.inc),1),2); % Daily incidence
        cinc1(:,ii,ip) = sum(soln6(:,i.aux.inc),2);        % Cumulative incidence 
        mor1(:,ii,ip) = sum(soln6(:,i.aux.mort),2);        % Cumulative deaths
        
        % Cases averted
        CA1(:,ii,ip) = sum((inc0(:,ii,1) - inc1(:,ii,ip)),1);
        % Percentage reduction of cumulative cases  
        Red_cumcases1(:,ii,ip) = 100*((cinc0(end,ii,1)- cinc1(end,ii,ip))./cinc0(end,ii,1));
        % Percentage reduction of cumulative deaths 
        Red_Deaths1(:,ii,ip) = 100*((mor0(end,ii,1)- mor1(end,ii,ip))./mor0(end,ii,1));
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    fprintf('\n');
end
  
% Impact of responsive vaccination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Percentage reduction of cumulative deaths
RedDeaths_prc = prctile(Red_Deaths,[2.5,50,97.5],2); 
x = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
x(:,ip) = RedDeaths_prc(:,:,ip);
end

 % Cases averted
CA_prc = prctile(CA,[2.5,50,97.5],2); 
y = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
y(:,ip) = CA_prc(:,:,ip);
end

% Percentage reduction of cumulative incidence 
RedCumcase_prc = prctile(Red_cumcases,[2.5,50,97.5],2); 
z = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
z(:,ip) = RedCumcase_prc(:,:,ip);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impact of responsive vaccination + NPI
% Percentage reduction of cumulative deaths
RedDeaths_prc1 = prctile(Red_Deaths1,[2.5,50,97.5],2); 
x1 = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
x1(:,ip) = RedDeaths_prc1(:,:,ip);
end

% Percentage reduction of cumulative incidence 
RedCumcase_prc1 = prctile(Red_cumcases1,[2.5,50,97.5],2); 
z1 = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
z1(:,ip) = RedCumcase_prc1(:,:,ip);
end

 % Cases averted
CA_prc1 = prctile(CA1,[2.5,50,97.5],2); 
y1 = zeros(3,length(th_TPR));
for ip = 1:length(th_TPR)
y1(:,ip) = CA_prc1(:,:,ip);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save data as bellow for Figure S1 and Table S1: 

% Severity reducing vaccine  (p.c1 = 0; p.c3 = 0.6)
% save simulate2p5_sev;     % R0 = 2.5;
% save simulate2_sev;       % R0 = 2.0

% Infection reducing vaccine  (p.c1 = 0.6; p.c3 = 0)
% save simulate2p5_sus;     % R0 = 2.5;
 save simulate2_sus;       % R0 = 2.0;




