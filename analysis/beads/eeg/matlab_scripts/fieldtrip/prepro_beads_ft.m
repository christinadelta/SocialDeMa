%%% Preprocessing TEST of the beads data with Fieldtrip

% file path
sub_path = '/Users/christinadelta/Desktop/os_data/beads';
sub_data = fullfile(sub_path, 'sub_08_beads_block_01.bdf');

% check the biosemi64 layout
% cfglayout = [];
% cfglayout.layout = 'biosemi64.lay';
% lo = ft_layoutplot(cfglayout);

%% load the data  and define trials 
cfg                     = [];
cfg.dataset             = sub_data;
cfg.trialfun            = 'trialfun_beads';
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = [1 2 3 4];
cfg.trialdef.prestim    = 0.2;
cfg.trialdef.poststim   = 0.8;
cfg                     = ft_definetrial(cfg);

% in this part of analysis draws will be used for epoching. Thus on
% every sequence I focus only on the draws. Every draw in a sequence is defined as a trial! 
% Draw epochs are stored in the "trl" cell created using the
% "trialfun_beads" function
draws           = cfg.trials{2};
conds           = length(draws); % easy[80:20] and difficult[60:40] conditions 

%% preprocess 

% baseline correction 
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];

% filtering 
% cfg.lpfilter        = 'yes';
% cfg.lpfreq          = 1;

% re-referencing 
cfg.reref           = 'yes';
cfg.refchannel      = {'EXG1' 'EXG2'};

% for preprocessing
% cfg.trialdef.eventtype = 'STATUS';
% cfg.trialdef.eventvalue = [1 2 3 4];
% cfg.trialdef.prestim = 0.2;
% cfg.trialdef.poststim = 0.8;


% extract trls from the cfg structure. each trl contains the draws [trials]
% to be epoched for each condition
for i = 1:conds
    
    this_cond{i}    = draws{i};
    
    % extract the current condition cell 
    tmp_cell        = this_cond{1};
    length_cell     = length(tmp_cell);
    
    % extract trls from that condition cell
    for j = 1:length_cell
        
        this_trl    = tmp_cell{j};
        cfg.trl     = this_trl;
        
        % preprocess this_trl
        data        = ft_preprocessing(cfg);
        
        clear this_trl cfg.trl
    end
  
end

cfg                 = [];
ft_databrowser(cfg, data)


