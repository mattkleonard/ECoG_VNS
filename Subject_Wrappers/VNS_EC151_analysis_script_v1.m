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
root_dir = '/Users/willschuerman/Documents/Research/Data/EC151/EC151_B32';
%dat_dir = [root_dir '/EC131_759989b0'];
dat_dir = [root_dir filesep 'h5'];
brain_dir = '/Users/willschuerman/Documents/Research/Data/EC151/BrainPlot';
patient_code = 'EC151';
out_dir = root_dir;
% 
% %% View bad channels:
% Pre-processed data - No hilbert transform

% TO DO: Look at the notching scripts, and all the steps individually 
% Look at raw data, VNS=off.
% figures suggest that notching isn't working right. 


% what does stimulation artifact look like on all channels, esp one used
% for reference?
% Use white matter electrodes?
% If notching doesn't work, can use white matter as CAR, subtract 10 random
% white matter electrodes from all other ones. Potentially better than
% notching. 

% 
VNS_raw = VNS_process_raw(patient_code, dat_dir, brain_dir, 3,[1,2,4:85], [],'load_RAW_erps',true,'load_env',false,'notch_data',true, 'VNS_duty_cycle', [2,28,0]);
[is_bad_chan] = plot_spectra_debug(VNS_raw, 3052, 0.6); % FIND EXACT SAMPLING RATE
visually_inspect_channels(VNS_raw); %also another (probably better) method
% of detecting bad channels
is_bad_chan = ones(84,1);
is_bad_chan([6:14, 16, 19:36, 42:72, 74:81, 84]) = 0;
clear VNS_raw

VNS_dat = VNS_process_raw(patient_code, dat_dir, brain_dir,1 ,2:85, [],'load_env',true,'VNS_duty_cycle', [2,28,0]);
visually_inspect_channels(VNS_dat);
VNS_dat.good_channels = ~is_bad_chan;
% % %% Find bad times:
VNS_dat = artifact_detector_3(VNS_dat, 8);
find_boolean_on(VNS_dat.is_bad_time_theta)

save([root_dir '/VNS_dat_B46_env.mat'],'VNS_dat' ,'-v7.3')

%% Can start from here if file exists
%load('/Users/willschuerman/Documents/Research/Data/EC131/VNS_dat_env.mat');
out_table_mean = generate_block_analysis(VNS_dat,10);  % check the arguments on the paircoeff, make sure I am doing it right
out_table_mean = generate_block_analysis(VNS_dat,30,[],[],true); % plot only significant erps...

[kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10,[-30 60]); % plot all channels
[ks_spectrograms, time_axis, f_axis] = ks_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10,[-30 60]); % plot all channels
