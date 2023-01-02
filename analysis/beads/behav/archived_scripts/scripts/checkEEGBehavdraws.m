function[allfrontal, allfc, allcp, allpar] = checkEEGBehavdraws(frontal, frontcent, centpar, par, sub_drawinfo, sub)

% splits eeg data into sequences and conditions 
c       = 1; % counter
e       = 1; % easy trials
d       = 1; % difficult trials


trls    = unique(sub_drawinfo(:,2));

if sub == 36
    
    sub_drawinfo(112:116,:) = []
end

for t = 1:length(trls)
    
    % what trial is it? how many draws to that trial?
    this_trl    = find(sub_drawinfo(:,2)== t);
    
%     if sub == 36 % there was a problem with sub 36 (in trial 27 first 5 draws where not recorded)
%         
%         if t == 27
%             this_trl = this_trl(1:5,1)
%         end
%     end
%    
    
    % get this_trials all eeg data 
    % frontal
    thiseeg_f   = frontal.f((this_trl),:);
    thiseeg_lf  = frontal.lf((this_trl),:);
    thiseeg_rf  = frontal.rf((this_trl),:);
    
    % frontocentral
    thiseeg_fc  = frontcent.fc((this_trl),:);
    thiseeg_lfc = frontcent.lfc((this_trl),:);
    thiseeg_rfc = frontcent.rfc((this_trl),:);
    
    % centroparietal
    thiseeg_cp  = centpar.cp((this_trl),:);
    thiseeg_lcp = centpar.lcp((this_trl),:);
    thiseeg_rcp = centpar.rcp((this_trl),:);
    
    % parietal 
    thiseeg_p   = par.p((this_trl),:);
    thiseeg_lp  = par.lp((this_trl),:);
    thiseeg_rp  = par.rp((this_trl),:);
    
    % store this-sequence events and eeg in correct condition cell
    if thiseeg_f(1,2) == 1 || thiseeg_f(1,2) == 2
        
        % update all frontal
        all_eegf{1,1}{1,e}  = thiseeg_f;
        all_eeglf{1,1}{1,e} = thiseeg_lf;
        all_eegrf{1,1}{1,e} = thiseeg_rf;
        
        % update all frontocentral
        all_eegfc{1,1}{1,e}     = thiseeg_fc;
        all_eeglfc{1,1}{1,e}    = thiseeg_lfc;
        all_eegrfc{1,1}{1,e}    = thiseeg_rfc;
        
        % update all cp
        all_eegcp{1,1}{1,e}     = thiseeg_cp;
        all_eeglcp{1,1}{1,e}    = thiseeg_lcp;
        all_eegrcp{1,1}{1,e}    = thiseeg_rcp;
        
        % update all p
        all_eegp{1,1}{1,e}     = thiseeg_p;
        all_eeglp{1,1}{1,e}    = thiseeg_lp;
        all_eegrp{1,1}{1,e}    = thiseeg_rp;
        
        
        e                       = e + 1;
    else
        all_eegf{1,2}{1,d}  = thiseeg_f;
        all_eeglf{1,2}{1,d} = thiseeg_lf;
        all_eegrf{1,2}{1,d} = thiseeg_rf;
        
        % update all frontocentral
        all_eegfc{1,2}{1,d}     = thiseeg_fc;
        all_eeglfc{1,2}{1,d}    = thiseeg_lfc;
        all_eegrfc{1,2}{1,d}    = thiseeg_rfc;
        
        % update all cp
        all_eegcp{1,2}{1,d}     = thiseeg_cp;
        all_eeglcp{1,2}{1,d}    = thiseeg_lcp;
        all_eegrcp{1,2}{1,d}    = thiseeg_rcp;
        
        % update all p
        all_eegp{1,2}{1,d}     = thiseeg_p;
        all_eeglp{1,2}{1,d}    = thiseeg_lp;
        all_eegrp{1,2}{1,d}    = thiseeg_rp;
   
        d                       = d + 1;
    end
    
end % end of trials loop

% store cells in structures 
allfrontal.f = all_eegf;
allfrontal.lf = all_eeglf;
allfrontal.rf = all_eegrf;

allfc.fc = all_eegfc;
allfc.lfc = all_eeglfc;
allfc.rfc = all_eegrfc;

allcp.cp = all_eegcp;
allcp.lcp = all_eeglcp;
allcp.rcp = all_eegrcp;

allpar.p = all_eegp;
allpar.lp = all_eeglp;
allpar.rp = all_eegrp;

end