subj = 'EC151';
rootdir = '/Users/mattleonard/Documents/Research/data';

%% FLAGS & PARAMS

preprocess_flag = 1;        % whether to preprocess data
    VNS_parms = [25 0.4 250];  % stim parameters used (frequency, amplitude, pw)
    VNS_duty_cycle = [2,28,0];  % timing of VNS duty cycle (ramp,duration,ramp)
    EKG_ch = [43];         % indices of EKG channels
    load_range = [];            % if empty, load all files in h5, otherwise load specified files
    notch_60Hz_flag = 1;        % whether to notch 60Hz + harmonics
    notch_VNS_flag = 1;         % whether to notch VNS artifact
    CARflag = 0;                % whether to CAR the data
    fsDs = 400;                 % desired downsampled fs
    ERP_times = [-10 10];       % window for ERPs
    h5_file_info = [];      % which h5 file(s) to use
    save_preproc_flag = 1;      % whether to save output
    outfile_name = [subj '_VNS_ERPs_' ...
        num2str(VNS_parms(1)) '_' ...
        num2str(VNS_parms(2)) '_' ...
        num2str(VNS_parms(3))]; % name of outfile

zscore_flag = 0;            % whether to z-score data (will only do once)
find_bad_data_flag = 0;     % whether to reject bad trials/channels
    thresh = 10;                % STD threshold for bad trials
    nBadTrialThresh = 20;       % nTrials threshold for bad channels
    remove_badTrial_flag = 1;   % whether to remove bad trials

plot_ERP_flag = 0;          % whether to plot ERP data
    ds_flag = 0;                % whether to downsample data further
    plot_avg_ERP_flag = 1;      % whether to plot average ERPs for all chans
    plot_raster_ERP_flag = 1;   % whether to plot rasters for all chans
    plot_avg_ERP_region_flag = 0; % whether to plot ERPs for anatomical region
    plot_single_chan_flag = 0;  % whether to plot ERP and raster for single elec

test_stats_flag = 0;        % whether to run statistics
    test_pre_post_dist_flag = 1;    % whether to compare pre- vs. post- VNS distributions
    test_ERP_timecourse_flag = 0;   % whether to test ERP timecourses
    plot_brain_flag = 1;            % whether to plot results on brain
    testStat = 'ranksum';           % which test to use ('ranksum', 'signrank', 'ttest', 'ttest2')
    twins = [-10 -5 ; 0 5];         % which time windows to use

trls = [];                  % which trials to use
elec = 6;                   % which elec to plot

recType = 'TDT';       % 'clinical' or 'TDT'
ch = [1 2 3 4 5 6 7 ...
    8 9 11 12 13 14 ...
    33 34 35 36 37 38 ...
    39 40 43 44 45 47 48 ...
    65 66 67 68 69 70 ...
    71 72 73 74 75 76 ...
    77 78 81 82 83 84];                 % indices of ECoG channels to use

%% SET UP DIRECTORIES

dat_dir = [rootdir filesep subj filesep 'VNS/h5/' num2str(VNS_parms(1)) '_' num2str(VNS_parms(2)) '_' num2str(VNS_parms(3))];
brain_dir = [rootdir filesep '..' filesep 'pia/data_store2/imaging/subjects'];
out_dir = [rootdir filesep subj filesep 'VNS'];

ch = [ch EKG_ch];

%% RUN VNS_process_raw_v2

if preprocess_flag
    VNS_dat = VNS_process_raw_v2(subj, dat_dir, brain_dir, EKG_ch, ch,...
        'load_range',load_range,...
        'fsDs',fsDs,...
        'ERP_times',ERP_times',...
        'recType',recType,...
        'stim_freq',VNS_parms(1),...
        'VNS_duty_cycle', VNS_duty_cycle,...
        'notch_60Hz_flag',notch_60Hz_flag,...
        'notch_VNS_flag',notch_VNS_flag,...
        'CARflag',CARflag,...
        'save_raw_erps',true,...
        'save_preproc_flag',save_preproc_flag,...
        'out_dir',out_dir,...
        'outfile_name',outfile_name,...
        'load_env',false,'notch_data',false);
elseif ~exist('VNS_dat','var')
    fprintf('Loading VNS_dat....\n');
    load([rootdir filesep subj filesep 'VNS' filesep outfile_name '.mat']);
end

%% REMOVE UNDESIRED BLOCKS

if ~isempty(h5_file_info) & ~exist('remove_bad_blocks_flag','var') & length(unique(VNS_dat.file_info)) > 1
    fprintf('Removing undesired blocks\n');
    VNS_dat.raw_erps(:,:,find(~strcmpi(VNS_dat.file_info,h5_file_info))) = [];
    VNS_dat.file_info(find(~strcmpi(VNS_dat.file_info,h5_file_info))) = [];
    remove_bad_blocks_flag = 1;
end

%% Z-SCORE DATA

if zscore_flag
    if isfield(VNS_dat,'zscore_flag')
        fprintf('Data have already been z-scored.\n');
    else
        VNS_dat = VNS_zscore_data(VNS_dat);
    end
end

%% REJECT BAD DATA

if find_bad_data_flag
    if isfield(VNS_dat,'badTrials')
        fprintf('Data have already been rejected.\n');
    else
        VNS_dat = VNS_find_bad_data(VNS_dat,thresh,nBadTrialThresh,h5_file_info,remove_badTrial_flag);
    end
end

%% PLOT ERP DATA

if plot_ERP_flag
    VNS_plot_ERPs(VNS_dat,trls,elec,...
        'ds_flag',ds_flag,...
        'plot_avg_ERP_flag',plot_avg_ERP_flag,...
        'plot_raster_ERP_flag',plot_raster_ERP_flag,...
        'plot_avg_ERP_region_flag',plot_avg_ERP_region_flag,...
        'plot_single_chan_flag',plot_single_chan_flag);
end

%%

if test_stats_flag
    [pval,statVal] = VNS_ERP_stats(VNS_dat,subj,brain_dir,...
        'recType',recType,...
        'test_pre_post_dist_flag',test_pre_post_dist_flag,...
        'test_ERP_timecourse_flag',test_ERP_timecourse_flag,...
        'plot_brain_flag',plot_brain_flag,...
        'testStat',testStat,...
        'twins',twins);
end

%%

