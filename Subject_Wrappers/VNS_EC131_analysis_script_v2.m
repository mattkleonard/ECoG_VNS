% Issues: 
% 2. Time windows for VNS clipped at beginning end of block
% - Solution: Increase window size to +-30
% changed make_vns_erps: raw_erp_win = [-10 10]; --> raw_erp_win = [-30
% 30]; Also find_vns_stim_on.m  stim_duration_thresh = 30; % minimum time in seconds of stim

% 3. Possible edge artifacts by doing hilb on individual blocks (need to cut out data?)
% - Potential Solution: Add script that cuts out VNS times with a buffer on
% each side, concatenating stims that are separated into different files
% using the raw data. 
% - Current Solution: NaNed out the first and last 20 samples from each
% portion.


% Need to output everything at each stage and see if it makes sense

% VNS process data script:
root_dir = '/Users/mattleonard/Documents/Research/data/EC131/VNS/EC131_QuietTime'; % '/Users/mattleonard/Documents/Research/pia/data_store1/human/prcsd_data/EC151/EC151_B46';
%dat_dir = [root_dir '/EC131_759989b0'];
dat_dir = [root_dir filesep 'h5'];
brain_dir = '/Users/mattleonard/Documents/Research/pia/data_store2/imaging/subjects';
patient_code = 'EC131';
out_dir = '/Users/mattleonard/Documents/Research/data/EC131';

% parameters
recType = 'clinical';
load([brain_dir '/' patient_code '/elecs/' recType '_elecs_all.mat']);

ch = 1:136;
vns_ch = [137 138];
ch = [ch vns_ch];

load_range = [];

CARflag = 0;

fsIn_fsOut = [1024 1024];
VNS_duty_cycle = [2,26,2];

%% 
% %% View bad channels:
% Pre-processed data - No hilbert transform
VNS_raw = VNS_process_raw(patient_code, dat_dir, brain_dir, vns_ch, ch, [],'load_full_dir',load_range,'load_RAW_erps',true,'load_env',false,'notch_data',false, 'VNS_duty_cycle', VNS_duty_cycle,'fsIn_fsOut',fsIn_fsOut,'recType',recType,'CARflag',CARflag);
%[is_bad_chan] = plot_spectra_debug(VNS_raw, 1024, 0.6); % atm, all
%channels are considered to be good
% visually_inspect_channels(VNS_raw); also another (probably better) method
% of detecting bad channels


% started at about 9:05am, ended at about
clear VNS_raw

% %% Find bad times:
VNS_dat = VNS_process_raw(patient_code, dat_dir, brain_dir,137 ,1:136, [],'load_env',true,'load_full_dir',[1]);
visually_inspect_channels(VNS_dat);
VNS_dat.good_channels = ~is_bad_chan;
% 
VNS_dat = artifact_detector_3(VNS_dat, 10);

save([root_dir '/VNS_dat_env.mat'],'VNS_dat' ,'-v7.3')

%% Can start from here if file exists
load('/Users/willschuerman/Documents/Research/Data/EC131/VNS_dat_env.mat');
out_table_mean = generate_block_analysis(VNS_dat,30);  % check the arguments on the paircoeff, make sure I am doing it right
out_table_mean = generate_block_analysis(VNS_dat,30,[],[],true); % plot only significant erps...

[kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10,[-30 60]); % plot all channels

