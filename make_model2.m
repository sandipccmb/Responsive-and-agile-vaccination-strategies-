% v2: Version to reflect updated model, i.e. with urbrur, age groups, and risk groups

function M = make_model2(p, r, i, s, gps, prm)


% -------------------------------------------------------------------------
% --- Get the linear rates ------------------------------------------------

m = zeros(i.nstates);
for iv = 1:length(gps.vac)
    vac = gps.vac{iv};
    for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};

               getaddr = @(st) i.(st).(vac).(risk).(age);
                
                U  = getaddr('U');                                         % Uninfected
                E  = getaddr('E');                                         % Exposed
                A  = getaddr('A');                                         % Asymptomatic
                P  = getaddr('P');                                         % Presymptomatic
                MS = getaddr('MS');                                        % Mild symptomatic
                R  = getaddr('R');                                         % Recovered
                
                % --- Development of symptoms (Mild or Severe)
                source  = P;
                destin  = MS;
                rate    = r.eta;
                m(destin, source) = m(destin, source) + rate;
                               
                % --- Recovering from disease
                sources = [A         MS];
                destin  = R;
                rate    = r.gamma;
                m(destin, sources) = m(destin, sources) + rate;              

        end
    end
end

for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};
            
                % --- Incubation (un-vaccinated)
                source  = i.E.v0.(risk).(age);
                destins = [i.A.v0.(risk).(age),      i.P.v0.(risk).(age)];
                rates   = [(1-p.sympto),             p.sympto]*r.incub;
                m(destins, source) = m(destins, source) + rates';
                
                % --- Incubation (vaccinated - intermediate)
                source  = i.E.vi.(risk).(age);
                destins = [i.A.vi.(risk).(age),      i.P.vi.(risk).(age)];
                rates   = [(1-p.sympto),             p.sympto]*r.incub;
                m(destins, source) = m(destins, source) + rates';
              
        end
end

for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};
            
                % --- Incubation
                source  = i.E.vf.(risk).(age);
                destins = [i.A.vf.(risk).(age),   i.P.vf.(risk).(age)];
                rates   = [1-p.sympto*(1-p.c3),   p.sympto*(1-p.c3)]*r.incub;
   
                m(destins, source) = m(destins, source) + rates';
        end
end


for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};
            
                % --- Vaccination (Rate of vaccination)
                source  = i.U.v0.(risk).(age);
                destin  = i.U.vi.(risk).(age);
                rate    = r.vacc(ir,ia);     
                m(destin, source) = m(destin, source) + rate;
                
                % --- Vaccination (Per-capita rate of vaccine induced immunity )
                source  = i.U.vi.(risk).(age);
                destin  = i.U.vf.(risk).(age);
                rate    = r.satu;            
                m(destin, source) = m(destin, source) + rate;
        end
end


M.lin = sparse(m - diag(sum(m,1)));


% -------------------------------------------------------------------------
% --- Get the nonlinear rates ---------------------------------------------

for ia = 1:length(gps.age)
    age = gps.age{ia};
      
        m = zeros(i.nstates);
        for iv = 1:length(gps.vac)
            vac = gps.vac{iv};
            for ir = 1:length(gps.risk)
                risk = gps.risk{ir};
                
                getaddr = @(st) i.(st).(vac).(risk).(age);
                U  = getaddr('U');                                         % Uninfected
                E  = getaddr('E');                                         % Exposed
                R  = getaddr('R');                                         % Recovered
                
                m(E, U) = 1;
                m(E, R) = 1-p.imm;
            end
        end
        
        m(:,s.vf) = m(:,s.vf)*(1-p.c1);
        M.nlin.(age) = sparse(m - diag(sum(m,1)));
%    end
end

% -------------------------------------------------------------------------
% --- Getting force-of-infection ------------------------------------------
% p.c1: Relative susceptibility
% p.c2: Relative infectiousness
% p.c: Relative infectiousness of asymptomatic vs symptomatic infection

% First, construct basic building block of repeating contact matrix, that can be weighted as necessary
templ = zeros(3,i.nstates);
templ(:,s.infectious) = repmat(prm.contact,1,length(s.infectious)/length(gps.age));

% Adjust for asymptomatics being less infectious
templ(:,s.A) = templ(:,s.A)*p.c;

% % Adjust for vaccinated people being less infectious
% templ(:,s.vf) = templ(:,s.vf)*(1-p.c2);

% Correct for population numbers
tmp2  = prm.N(:)';
den   = repmat(tmp2, size(templ,1), length(s.infectious)/length(gps.age)); 
templ(:,s.infectious) = templ(:,s.infectious)./den;

m = r.beta*templ;
M.lam = sparse(m);

% --- Get the mortality rates
m = zeros(1,i.nstates);
ii = 1;
for ir = 1:length(gps.risk)
    risk = gps.risk{ir};
    for ia = 1:length(gps.age)
        age = gps.age{ia};
        inds = intersect(s.MS,intersect(s.(risk),s.(age)));
        m(inds) = r.mu(ii); ii = ii+1;
    end
end
M.mortvec = m';
