function [allfrontal, allfc, allcp, allpar] = eegExtract(sub_eeg, sub_drawinfo, sub)

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
fdata       = data([1,2,3,4,18,19,20,21,22],:,:);       % frontal channles
lfdata      = data([1,2,3,4,18],:,:);                   % left frontal channels
rfdata      = data(19:22,:,:);                          % right frontal channels

fcdata      = data([5,6,7,23,24,25,26],:,:);            % frontocentral channels
lfcdata     = data([5,6,7,26],:,:);                     % left frontocentral channels
rfcdata     = data([23,24,25],:,:);                     % right frontocentral channels

cpdata      = data([8,9,10,17,27,28,29],:,:);           % centroparietal channels
lcpdata     = data([8,9,10,17],:,:);                    % left centroparietal channels
rcpdata     = data([27,28,29],:,:);                     % right centroparietal channels

pdata       = data([11,12,13,14,15,16,30,31,32,33,34],:,:);   % parietal channels
lpdata      = data(11:16,:,:);                          % left parietal channels
rpdata      = data(30:34,:,:);                          % right parietal channels

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

%% store in structures

frontal.f       = av_f;
frontal.lf      = av_lf;
frontal.rf      = av_rf;

frontcent.fc     = av_fc;
frontcent.lfc   = av_lfc;
frontcent.rfc   = av_rfc;

centpar.cp      = av_cp;
centpar.lcp     = av_lcp;
centpar.rcp     = av_rcp;

par.p           = av_p;
par.lp          = av_lp;
par.rp          = av_rp;

%% split EEG data in sequences and conditions 

% based on draw info split EEG data into sequences and conditions and check
% whether length for draws and for events is the same 
[allfrontal, allfc, allcp, allpar] = checkEEGBehavdraws(frontal, frontcent, centpar, par, sub_drawinfo, sub);

end