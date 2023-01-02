function [new_AQs, totalf, totalfc, totalcp, totalpar] = checkValuesV2(allfrontal, allfc, allcp, allpar, all_AQdiffs, sub)

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

% create new arrays for all EEG data
% frontal
f               = nan(1,1);
lf              = nan(1,1);
rf              = nan(1,1);

% fc
fc              = nan(1,1);
lfc             = nan(1,1);
rfc             = nan(1,1);

% cp
cp              = nan(1,1);
lcp             = nan(1,1);
rcp             = nan(1,1);

% par
p               = nan(1,1);
lp              = nan(1,1);
rp              = nan(1,1);

% AQ difference values
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
    
    % extract condition AQ difference values
    cond_AQ     = all_AQdiffs{1,cond};
    
    for ii = 1:length(cond_f) % for trials 
        
        % frontal
        this_f      = cond_f{1,ii}(:,1);
        this_lf     = cond_lf{1,ii}(:,1);
        this_rf     = cond_rf{1,ii}(:,1);
        
        % fc
        this_fc     = cond_fc{1,ii}(:,1);
        this_lfc    = cond_lfc{1,ii}(:,1);
        this_rfc    = cond_rfc{1,ii}(:,1);
        
        % cp
        this_cp     = cond_cp{1,ii}(:,1);
        this_lcp    = cond_lcp{1,ii}(:,1);
        this_rcp    = cond_rcp{1,ii}(:,1);
        
        % par
        this_p      = cond_p{1,ii}(:,1);
        this_lp     = cond_lp{1,ii}(:,1);
        this_rp     = cond_rp{1,ii}(:,1);
        
        % now exactract AQs
        % if there are 
        this_AQ     = cond_AQ{1,ii};
        
        % if 
        if size(this_AQ,2) ~= 1
            this_AQ = this_AQ'; % transpose
        end
        
        if sub == 36 && cond == 2 && ii == 14
            
            this_AQ(1:5,:) = [];
        end
        
        % remove nans from AQs (if there are any)
        this_AQ(isnan(this_AQ(:,1)),:) = [];
        
        % check new length of AQ array and if inconsistent with length of
        % EEG arrays, remove last values (from EEG arrays)
        tmp_len = length(this_f);
        
        if tmp_len ~= length(this_AQ)
            
            % compute difference between them
            tmp_diff = tmp_len - length(this_AQ)
            
            % frontal
            new_f   = this_f(tmp_diff+1:end);
            new_lf  = this_lf(tmp_diff+1:end);
            new_rf  = this_rf(tmp_diff+1:end);
            
            % fronto-central
            new_fc      = this_fc(tmp_diff+1:end);
            new_lfc     = this_lfc(tmp_diff+1:end);
            new_rfc     = this_rfc(tmp_diff+1:end);
            
            % centro-parietal
            new_cp      = this_cp(tmp_diff+1:end);
            new_lcp     = this_lcp(tmp_diff+1:end);
            new_rcp     = this_rcp(tmp_diff+1:end);
            
            % parietal
            new_p       = this_p(tmp_diff+1:end);
            new_lp      = this_lp(tmp_diff+1:end);
            new_rp      = this_rp(tmp_diff+1:end);
            
        else
            % frontal
            new_f       = this_f;
            new_lf      = this_lf;
            new_rf      = this_rf;
            
            % fc
            new_fc      = this_fc;
            new_lfc     = this_lfc;
            new_rfc     = this_rfc;
            
            % cp
            new_cp      = this_cp;
            new_lcp     = this_lcp;
            new_rcp     = this_rcp;
            
            % par
            new_p       = this_p;
            new_lp      = this_lp;
            new_rp      = this_rp;
            
        end
        
        % loop over draws 
        for iii = 1:length(new_f) % for draws 
            
            % frontal
            f(counter,1)    = new_f(iii);
            lf(counter,1)   = new_lf(iii);
            rf(counter,1)   = new_rf(iii);
            
            % fc
            fc(counter,1)   = new_fc(iii);
            lfc(counter,1)  = new_lfc(iii);
            rfc(counter,1)  = new_rfc(iii);
            
            % cp
            cp(counter,1)   = new_cp(iii);
            lcp(counter,1)  = new_lcp(iii);
            rcp(counter,1)  = new_rcp(iii);
            
            % par
            p(counter,1)    = new_p(iii);
            lp(counter,1)   = new_lp(iii);
            rp(counter,1)   = new_rp(iii);
            
            new_AQs(counter,1) = this_AQ(iii);
            
            counter = counter + 1

        end % end of draws loop (iii)
        
    end % end of trials loop (ii)
    
end % end of conditions loop

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