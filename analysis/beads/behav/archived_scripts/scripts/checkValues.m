function [new_AQs, totalf, totalfc, totalcp, totalpar] = checkValues(allfrontal, allfc, allcp, allpar, all_AQdiffs)

% this function exctracts both AQ difference values and all the EEG
% averaged arrays and checks wether they are of the same length (i.e., if
% the model was drawing as much as the participants). In case the model's
% AQs or the EEG data is less, we add zeros

% OUTPUT: new arrays with all data being the of same length 

% ---------------------

% extract all EEG data from structures and AQ differences and check whether
% they are of the same length

cond_counter        = 2;
counter             = 1;

f               = nan(1,1);
lf              = nan(1,1);
rf              = nan(1,1);

fc              = nan(1,1);
lfc             = nan(1,1);
rfc             = nan(1,1);

cp              = nan(1,1);
lcp             = nan(1,1);
rcp             = nan(1,1);

p               = nan(1,1);
lp              = nan(1,1);
rp              = nan(1,1);
new_AQs         = nan(1,1);

% extract cells
% frontal 
tf               = allfrontal.f;
tlf              = allfrontal.lf;
trf              = allfrontal.rf;

% frontocentral
tfc              = allfc.fc;
tlfc             = allfc.lfc;
trfc             = allfc.rfc;

% centroparietal
tcp              = allcp.cp;
tlcp             = allcp.lcp;
trcp             = allcp.rcp;

% parieatl
tp              = allpar.p;
tlp             = allpar.lp;
trp             = allpar.rp;

for cond = 1:cond_counter
    
    % all frontal
    cond_f      = tf{1,cond};
    cond_lf     = tlf{1,cond};
    cond_rf     = trf{1,cond};
    
    % all frontocentral
    cond_fc     = tfc{1,cond};
    cond_lfc    = tlfc{1,cond};
    cond_rfc    = trfc{1,cond};
    
    % all cp
    cond_cp     = tcp{1,cond};
    cond_lcp    = tlcp{1,cond};
    cond_rcp    = trcp{1,cond};
    
    % all p 
    cond_p      = tp{1,cond};
    cond_lp     = tlp{1,cond};
    cond_rp     = trp{1,cond};
    
    cond_AQ     = all_AQdiffs{1,cond};
    
    for ii = 1:length(cond_f) % for trials 
        
        % frontal
        this_f = cond_f{1,ii}(:,1);
        this_lf = cond_lf{1,ii}(:,1);
        this_rf = cond_rf{1,ii}(:,1);
        
        % fc
        this_fc = cond_fc{1,ii}(:,1);
        this_lfc = cond_lfc{1,ii}(:,1);
        this_rfc = cond_rfc{1,ii}(:,1);
        
        % cp
        this_cp = cond_cp{1,ii}(:,1);
        this_lcp = cond_lcp{1,ii}(:,1);
        this_rcp = cond_rcp{1,ii}(:,1);
        
        % par
        this_p = cond_p{1,ii}(:,1);
        this_lp = cond_lp{1,ii}(:,1);
        this_rp = cond_rp{1,ii}(:,1);
        
        this_AQ = cond_AQ{1,ii}
        
        for iii = 1:length(this_f) % for draws 
            
            % frontal
            f(counter,1) = this_f(iii);
            lf(counter,1) = this_lf(iii);
            rf(counter,1) = this_rf(iii);
            
            % fc
            fc(counter,1) = this_fc(iii);
            lfc(counter,1) = this_lfc(iii);
            rfc(counter,1) = this_rfc(iii);
            
            % cp
            cp(counter,1) = this_cp(iii);
            lcp(counter,1) = this_lcp(iii);
            rcp(counter,1) = this_rcp(iii);
            
            % par
            p(counter,1) = this_p(iii);
            lp(counter,1) = this_lp(iii);
            rp(counter,1) = this_rp(iii);
            
            new_AQs(counter,1) = this_AQ(iii);
            
            counter = counter + 1
        end

    end % end of ii (sequence loop)
    
end % end of condition

% store all frontal
totalf.f    = f;
totalf.lf   = lf;
totalf.rf   = rf;

% store all fc
totalfc.fc  = fc;
totalfc.lfc = lfc;
totalfc.rfc = rfc;

% store all cp
totalcp.cp  = cp;
totalcp.lcp = lcp;
totalcp.rcp = rcp;

% store all par
totalpar.p   = p;
totalpar.lp = lp;
totalpar.rp = rp;



end