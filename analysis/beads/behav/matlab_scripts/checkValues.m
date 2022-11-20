function [new_AQs, f, fc, cp, par] = checkValues(allfrontal, allfc, allcp, allpar, all_AQdiffs)

% this function exctracts both AQ difference values and all the EEG
% averaged arrays and checks wether they are of the same length (i.e., if
% the model was drawing as much as the participants). In case the model's
% AQs or the EEG data is less, we add zeros

% OUTPUT: new arrays with all data being the of same length 

% ---------------------

% extract all EEG data from structures
% frontal
ft          = allfrontal.f;
lft         = allfrontal.lf;
rft         = allfrontal.rf;

f_one       = ft{1,1}(:,1);
f_two       = ft{1,2}(:,1);
fcon        = cat(1,f_one, f_two); % frontal concantinated

fl_one       = lft{1,1}(:,1);
fl_two       = lft{1,2}(:,1);
flcon        = cat(1,fl_one, fl_two); % left frontal concantinated

fr_one       = rft{1,1}(:,1);
fr_two       = rft{1,2}(:,1);
frcon        = cat(1,fr_one, fr_two); % right frontal concantinated

% frnontocentral
fct         = allfc.fc;
lfct        = allfc.lfc;
rfct        = allfc.rfc;

fc_one       = fct{1,1}(:,1);
fc_two       = fct{1,2}(:,1);
fccon        = cat(1,fc_one, fc_two); % frontocentral concantinated

fcl_one       = lfct{1,1}(:,1);
fcl_two       = lfct{1,2}(:,1);
fclcon        = cat(1,fcl_one, fcl_two); % left frontocentral concantinated

fcr_one       = rfct{1,1}(:,1);
fcr_two       = rfct{1,2}(:,1);
fcrcon        = cat(1,fcr_one, fcr_two); % right frontal central concantinated

% centroparietal
cpt         = allcp.cp;
lcpt        = allcp.lcp;
rcpt        = allcp.rcp;

cp_one       = cpt{1,1}(:,1);
cp_two       = cpt{1,2}(:,1);
cpcon        = cat(1,cp_one, cp_two); % centroparietal concantinated

cpl_one       = lcpt{1,1}(:,1);
cpl_two       = lcpt{1,2}(:,1);
cplcon        = cat(1,cpl_one, cpl_two); % left centroparietal concantinated

cpr_one       = rcpt{1,1}(:,1);
cpr_two       = rcpt{1,2}(:,1);
cprcon        = cat(1,cpr_one, cpr_two); % right centroparietal concantinated

% parietal
pt          = allpar.p;
lpt         = allpar.lp;
rpt         = allpar.rp;

p_one       = pt{1,1}(:,1);
p_two       = pt{1,2}(:,1);
pcon        = cat(1,p_one, p_two); % parietal concantinated

pl_one       = lpt{1,1}(:,1);
pl_two       = lpt{1,2}(:,1);
plcon        = cat(1,pl_one, pl_two); % left parietal concantinated

pr_one       = rpt{1,1}(:,1);
pr_two       = rpt{1,2}(:,1);
prcon        = cat(1,pr_one, pr_two); % parietal concantinated

% extract condition AQ differences 
aqs_one     = all_AQdiffs{1,1}';
aqs_two     = all_AQdiffs{1,2}';
AQs         = cat(1, aqs_one, aqs_two); % concantinate

if length(fcon) < length(AQs)
    
    diff_values = length(AQs) - length(fcon); % compute difference between lengths
    oldend      = length(fcon);
    
    % fill EEG arrays with zeros
    for i = 1:diff_values 
        
        % fill frontal arrays
        fcon(oldend+i)      = zeros;
        flcon(oldend+i)     = zeros;
        frcon(oldend+i)     = zeros;
        
        % fill frontocentral arrays
        fccon(oldend+i)     = zeros;
        fclcon(oldend+i)    = zeros;
        fcrcon(oldend+i)    = zeros;
        
        % fill centroparietal
        cpcon(oldend+i)     = zeros;
        cplcon(oldend+i)    = zeros;
        cprcon(oldend+i)    = zeros;
        
        % fill parietal
        pcon(oldend+i)      = zeros;
        plcon(oldend+i)     = zeros;
        prcon(oldend+i)     = zeros;

    end
    
elseif length(fcon) > length(AQs)
    
    diff_values = length(fcon) - length(AQs);
    oldend      = length(AQs);
    
    % fill AQS differences array with zeros
    for i = 1:diff_values 
        
        AQs(oldend+i)     = zeros;
        
    end
end % end of if statement

% re arrange all arrays four outputs
new_AQs     = AQs;
f.f         = fcon;
f.fl        = flcon;
f.fr        = frcon;

fc.fc       = fccon;
fc.lfc      = fclcon;
fc.rfc      = fcrcon;

cp.cp       = cpcon;
cp.lcp      = cplcon;
cp.rcp      = cprcon;

par.p       = pcon;
par.lp      = plcon;
par.rp      = prcon;

end