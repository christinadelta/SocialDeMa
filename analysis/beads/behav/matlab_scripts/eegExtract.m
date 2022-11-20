function [allfrontal, allfc, allcp, allpar] = eegExtract(sub_eeg)

% this function is part of the BEADS formal analysis
% Inputs: - a double (subject's number of draws)
%         - structure with cropped EEG data

% second input (sub_eeg structure) contains data (channels x samples x
% epochs), conditions and events
% we need all these to seperate data for different channels and to re-order
% them (condition-wise), and average by epochs/draws

% lastly, the function splits the data into the two conditions (0.6, 0.8) 

% Output: 4 structures based on scalp sites (frontal, fc, cp, parietal) 
% all 4 structs contain: all (site) averaged data, left and right (sites)
% averaged data. These will be used to run regressions (with AQ difference as predictor)

% ---------------------------------------------

data        = sub_eeg.cropped_struct.data; 
events      = sub_eeg.cropped_struct.events; 

% separete MEEG data based on channels 
fdata       = data([1,2,3,4,17,18,19,20,21],:,:);       % frontal channles
lfdata      = data([1,2,3,4,17],:,:);                   % left frontal channels
rfdata      = data(18:21,:,:);                          % right frontal channels

fcdata      = data([5,6,7,25,22,23,24],:,:);            % frontocentral channels
lfcdata     = data([5,6,7,25],:,:);                     % left frontocentral channels
rfcdata     = data([22,23,24],:,:);                     % right frontocentral channels

cpdata      = data([8,9,10,16,26,27,28],:,:);           % centroparietal channels
lcpdata     = data([8,9,10,16],:,:);                    % left centroparietal channels
rcpdata     = data([26,27,28],:,:);                     % right centroparietal channels

pdata       = data([11,12,13,14,15,29,30,31,32],:,:);   % parietal channels
lpdata      = data(11:15,:,:);                          % left parietal channels
rpdata      = data(29:32,:,:);                          % right parietal channels

%% average the 3d matricies by epochs (3rd dimension)

% average data by epochs (this should result in an array with lenght
% equal subdraw). One data point for each epoch/draw
% % % frontal
av_f(:,1)       = mean(reshape(fdata, [], size(fdata,3)));
av_f(:,2)       = events;
av_lf(:,1)      = mean(reshape(lfdata, [], size(lfdata,3)));
av_lf(:,2)      = events;
av_rf(:,1)      = mean(reshape(rfdata, [], size(rfdata,3)));
av_rf(:,2)      = events;

% % % frontocentral
av_fc(:,1)      = mean(reshape(fcdata, [], size(fcdata,3)));
av_fc(:,2)      = events;
av_lfc(:,1)     = mean(reshape(lfcdata, [], size(lfcdata,3)));
av_lfc(:,2)     = events;
av_rfc(:,1)     = mean(reshape(rfcdata, [], size(rfcdata,3)));
av_rfc(:,2)     = events;

% % % centro-parietal
av_cp(:,1)      = mean(reshape(cpdata, [], size(cpdata,3)));
av_cp(:,2)      = events;
av_lcp(:,1)     = mean(reshape(lcpdata, [], size(lcpdata,3)));
av_lcp(:,2)     = events;
av_rcp(:,1)     = mean(reshape(rcpdata, [], size(rcpdata,3)));
av_rcp(:,2)     = events;

% % % parietal
av_p(:,1)       = mean(reshape(pdata, [], size(pdata,3)));
av_p(:,2)       = events;
av_lp(:,1)      = mean(reshape(lpdata, [], size(lpdata,3)));
av_lp(:,2)      = events;
av_rp(:,1)      = mean(reshape(rpdata, [], size(rpdata,3)));
av_rp(:,2)      = events;

%% split data into conditions 

% frontal sites 
tmp_f           = find(av_f(:,2) < 3); % find frontal events 1,2
cond_f{1,1}     = av_f((tmp_f),:); clear tmp_f
tmp_f           = find(av_f(:,2) > 2); % find frontal events 3,4
cond_f{1,2}     = av_f((tmp_f),:); clear tmp_f

tmp             = find(av_lf(:,2) < 3); % find left frontal events 1,2
cond_lf{1,1}    = av_lf((tmp),:); clear tmp
tmp             = find(av_lf(:,2) > 2); % find left frontal events 3,4
cond_lf{1,2}    = av_lf((tmp),:); clear tmp

tmp             = find(av_rf(:,2) < 3); % find right frontal events 1,2
cond_rf{1,1}    = av_rf((tmp),:); clear tmp
tmp             = find(av_rf(:,2) > 2); % find right frontal events 3,4
cond_rf{1,2}    = av_rf((tmp),:); clear tmp

% frontocentral sites 
tmp             = find(av_fc(:,2) < 3); % find frontocentral events 1,2
cond_fc{1,1}    = av_fc((tmp),:); clear tmp
tmp             = find(av_fc(:,2) > 2); % find frontocentral events 3,4
cond_fc{1,2}    = av_fc((tmp),:); clear tmp

tmp             = find(av_lfc(:,2) < 3); % find left frontocentral events 1,2
cond_lfc{1,1}   = av_lfc((tmp),:); clear tmp
tmp             = find(av_lfc(:,2) > 2); % find left frontocentral events 3,4
cond_lfc{1,2}   = av_lfc((tmp),:); clear tmp

tmp             = find(av_rfc(:,2) < 3); % find right frontocentral events 1,2
cond_rfc{1,1}   = av_rfc((tmp),:); clear tmp
tmp             = find(av_rfc(:,2) > 2); % find right frontocentral events 3,4
cond_rfc{1,2}   = av_rfc((tmp),:); clear tmp

% centroparietal sites 
tmp             = find(av_cp(:,2) < 3); % find centroparietal events 1,2
cond_cp{1,1}    = av_cp((tmp),:); clear tmp
tmp             = find(av_cp(:,2) > 2); % find centroparietal events 3,4
cond_cp{1,2}    = av_cp((tmp),:); clear tmp

tmp             = find(av_lcp(:,2) < 3); % find left centroparietal events 1,2
cond_lcp{1,1}   = av_lcp((tmp),:); clear tmp
tmp             = find(av_lcp(:,2) > 2); % find left centroparietal events 3,4
cond_lcp{1,2}   = av_lcp((tmp),:); clear tmp    
    
tmp             = find(av_rcp(:,2) < 3); % find right centroparietal events 1,2
cond_rcp{1,1}   = av_rcp((tmp),:); clear tmp
tmp             = find(av_rcp(:,2) > 2); % find right centroparietal events 3,4
cond_rcp{1,2}   = av_rcp((tmp),:); clear tmp        

% parietal sites 
tmp             = find(av_p(:,2) < 3); % find parietal events 1,2
cond_p{1,1}     = av_p((tmp),:); clear tmp_f
tmp             = find(av_p(:,2) > 2); % find parietal events 3,4
cond_p{1,2}     = av_p((tmp),:); clear tmp_f

tmp             = find(av_lp(:,2) < 3); % find left parietal events 1,2
cond_lp{1,1}    = av_lp((tmp),:); clear tmp
tmp             = find(av_lp(:,2) > 2); % find left parietal events 3,4
cond_lp{1,2}    = av_lp((tmp),:); clear tmp

tmp             = find(av_rp(:,2) < 3); % find right parietal events 1,2
cond_rp{1,1}    = av_rp((tmp),:); clear tmp
tmp             = find(av_rp(:,2) > 2); % find right parietal events 3,4
cond_rp{1,2}    = av_rp((tmp),:); clear tmp

% add them all to sites-specific structures
allfrontal.f    = cond_f;
allfrontal.lf   = cond_lf;
allfrontal.rf   = cond_rf;

allfc.fc        = cond_fc;
allfc.lfc       = cond_lfc;
allfc.rfc       = cond_rfc;

allcp.cp        = cond_cp;
allcp.lcp       = cond_lcp;
allcp.rcp       = cond_rcp;

allpar.p        = cond_p;
allpar.lp       = cond_lp;
allpar.rp       = cond_rp;


end