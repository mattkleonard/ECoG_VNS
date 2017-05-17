function VNS_dat = VNS_process_raw(patient_code, root_dir, brain_dir, ekg_ch,electrodes, varargin)

% This is the primary function for loading raw hdf5 files and converting
% them into useful ECoG measurements such as hilbert transformed envelope,
% phase, ERPs.
% The output is a structured array (VNS_dat) that contains these measurements
% as well as information about the timing of VNS stimulation.
% This data structure is commonly used in subsequent data analysis functions.
%
% This function is equivalent to VNS_prelim_script_9, however the name has been
% changed to simplify further processing.
%
% Determines when VNS Stimulation is on/off based on EKG
% Runs basic analyses comparing the ECoG signal with/without stimulation
%
% The data from a block is loaded, processed into the various bands and
% downsampled block by block.
% This saves memory space and speeds up the time required to
% generate these analyses
% (roughly speaking the computational time of the band processing is
% n*log(n) - if we have n blocks of length m, the comparison is between
% n(m*log(m)), vs (nm)*log(nm), or nmlog(m) + nm(log(n)).
%
% Additionally, the block data are now preallocated in memory such that
% matrcies no longer grow through concatenation.
%
% as a default the electrodes are sorted such that all electrodes are
% listed by the alphabetical order of their anatomy. This was justifed
% because a large portion of the electrodes were non-grid electrodes.

% Edited on  4/20/2017 by Will Schuerman
%
%
%% Inputs:
%
% patient_code - the code for the patient, e.g., EC131
%
% root_dir - the directory where raw data is stored
% (EX: '/Users/changlab/Documents/data/VNS/EC131/EC131_05104e39')
%
% ekg_ch - CH containing EKG data can vary per subject (EX: 137)
%
% electrodes - list of relevant electrodes (E:x 1:136 or [44, 45])
%
%% Variable Inputs
% 1 relevant_time_window % 2x2 matrix of start times in hours minutes
% [hours_start, minutes_start; hours_stop, minutes_stop]
% will only use files whose times start within the specified windows
% This should be left empty if the full directory is desired.
%
%% Flag Inputs:
% 'notch_data' - default true (Notch VNS Harmonics)
%
% 'stim_freq' - default 25 (the frequency of VNS stimulation pulses)
%
% 'VNS_duty_cycle' - default [30,2] (the timing in seconds of the power of
% a single pulse of VNS stimulation)
%
% 'load_full_dir' - default load all .h5 files in root_dir,
% input range of files to load (ex 1:4:30) or (1:25);
%
% 'load_phase' - default false (load phase data)
%
% 'load_env' - default true (load envelope data)
%
% 'load_vns_band' - default false (load vns harmonics band)
%
% 'load_RAW' - default false  (loads all raw data for channels (will be huge) - shorten)
%
% 'load_RAW_erps' - default false (loads raw data as erps (also big)
%
% 'load_peak_freqs' - default false (loads the peak frequency for each band
% for the off and on intervals)
%
%% Output:
% VNS_dat - the standard data structure for subsequent VNS analysis. By
% default contains timing information on Stimulations

apply_twin_roi = false;
if length(varargin)>0
    if ~isempty(varargin{1});
        apply_twin_roi = true;
        twin_roi = varargin{1};
        load_full_dir = true;
    end
end

%% Load Flags:
% notch data: notch data to remove the VNS harmonics
if ~isempty(find(strcmpi(varargin,'notch_data')))
    notch_flag = varargin{find(strcmpi(varargin,'notch_data'))+1};
else
    notch_flag = true;
end

% stimulation frequency:
if ~isempty(find(strcmpi(varargin,'stim_freq')))
    stim_freq = varargin{find(strcmpi(varargin,'stim_freq'))+1};
else
    stim_freq = 25;
end

% VNS Duty cycle:
if ~isempty(find(strcmpi(varargin, 'VNS_duty_cycle')))
    VNS_duty_cycle = varargin{find(strcmpi(varargin, 'VNS_duty_cycle'))+1};
else
    VNS_duty_cycle = [30,2];
end


% Loads all blocks of Clinical data in Root_dir (otherwise specify a range)
if ~isempty(find(strcmpi(varargin,'load_full_dir')))
    load_full_dir = false;
    load_range = varargin{find(strcmpi(varargin,'load_full_dir'))+1};
else
    load_full_dir = true;
end

% Load ECoG Analytic Amplidutes of Theta - High Gamma bands:
if ~isempty(find(strcmpi(varargin,'load_env')))
    load_envelope = varargin{find(strcmpi(varargin,'load_env'))+1};
else
    load_envelope = true;       % Load the bands for envelope
end

% Load ECoG Phase Angle of Theta - High Gamma bands:
if ~isempty(find(strcmpi(varargin,'load_phase')))
    load_phase = varargin{find(strcmpi(varargin,'load_phase'))+1};
else
    load_phase = false;          % Load the bands for phase data
end

load_vns_phase = false;

% Load the intensity of ecog bandpassed for VNS harmonics
if ~isempty(find(strcmpi(varargin,'load_vns_band')))
    load_vns_band_env = varargin{find(strcmpi(varargin,'load_vns_band'))+1};
else
    load_vns_band_env = false;
end

% load raw data:
% Load the combined Raw data for all blocks (HUGE)
if ~isempty(find(strcmpi(varargin,'load_RAW')))
    load_raw = varargin{find(strcmpi(varargin,'load_raw'))+1};
else
    load_raw = false;
end
% generate erps of the RAW data about onsets
if ~isempty(find(strcmpi(varargin,'load_RAW_erps')))
    load_raw_erps = varargin{find(strcmpi(varargin,'load_raw_erps'))+1};
else
    load_raw_erps = false;
end

% load peak frequency from the raw:
if ~isempty(find(strcmpi(varargin,'load_peak_freqs')));
    load_peak_freqs = varargin{find(strcmpi(varargin,'load_peak_freqs'))+1};
else
    load_peak_freqs = false;
end

blank_singular_anatomy = false; % ignore the names of anatomical regions that are found only once
load_stim_freq_raw = false; % load the stimulation frequency harmonics





%% Load anatomy file:
anatomy_file = [brain_dir filesep patient_code '/elecs/clinical_elecs_all.mat'];
if exist(anatomy_file) == 2
    load([brain_dir filesep patient_code '/elecs/clinical_elecs_all.mat']);
    anatomy_elecs = anatomy(electrodes,4);
else
    anatomy_elecs = strsplit(num2str(1:length(electrodes)));
end

%% Blank Anatomy with single type
anatomy_list = unique(anatomy_elecs);
if blank_singular_anatomy
    for i = 1:length(anatomy_list)
        is_anatomy = strcmpi(anatomy_elecs,anatomy_list{i});
        if sum(is_anatomy) == 1
            anatomy_elecs(is_anatomy) = {'single electrode in region'};
        end
    end
end
[anatomy_elecs, electrode_order] = sort(anatomy_elecs);
for i = 1:length(anatomy_elecs)
    anatomy_elecs{i} = strrep(anatomy_elecs{i}, '_',' ');
end
electrodes = electrodes(electrode_order); % group electrodes by anatomy


%% Load list of h5 filenames
filenames = dir(root_dir);
filenames(1:2) = []; % on windows, only two files
%filenames(1:4) = []; % on mac, this needs an extra 2 files????

for k = 1:length(filenames)
    test_files{k} = filenames(k).name;
end
for k = 1:length(test_files)
    tmp = test_files{k};
    is_h5(k) = strcmpi(tmp((end-2):end),'.h5');
end
test_files(~is_h5) = [];

if ~load_full_dir
    %file_range = 1:ceil(length(test_files)/2); % tentative reduced list...
    %file_range = [1 1:37 37];
    %file_range = 31:length(test_files);
    %file_range = 1:4:length(test_files);
    file_range = load_range;
    test_files = test_files(file_range);
end

%% Find Timing Data;
%twin = [13 15; 14 35]; % 2ma_500us specified in military time. Row 1 hours, mins (start), Row2 hours, mins_stop
%twin = [9 15; 10 35]; % 1ma_250us
%twin = [14 45; 16 0]; %2.25ma_250us
%twin = [8 0; 9 14]; % baseline
% twin = [11 0;12 30]; % 10hz
if apply_twin_roi
    t_start = twin_roi(1,1)*60+twin_roi(1,2); %Minutes since midnight
    t_stop = twin_roi(2,1)*60+twin_roi(2,2); %minutes since midnight

    is_in_twin = false(size(test_files));
    for k = 1:length(test_files)
        timestamps = hdf5read([root_dir filesep test_files{k}], 'timestamp vector');
        tsvec_1 = datevec((timestamps(1) - (7*3600))/86400 + datenum(1970,1,1));
        tsvec_1 =tsvec_1(4)*60+tsvec_1(5); % minutes since midnight
        if (tsvec_1 >= t_start) & (tsvec_1 <= t_stop)
            is_in_twin(k) = true;
        end
    end
    test_files(~is_in_twin) = [];
end

timestamps = hdf5read([root_dir filesep test_files{1}], 'timestamp vector');
duration = timestamps(end) - timestamps(1);
%tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));
%duration=(tsvec(end,5)*60+tsvec(end,6))-(tsvec(1,5)*60+tsvec(1,6));%total time
%fp = round(length(timestamps)/duration); %samples/sec (must be integer)
fileInfo = h5info([filesep root_dir filesep test_files{1}]);
fileInfo = h5info([filesep root_dir filesep test_files{k}]);
[p,q] = rat(1024/round(fileInfo.Datasets(1).Attributes(2).Value));

%% Process Data Parameters:
%fp = 1024;
fs_in = 1024;
fs_out = 100;  % down sample frequency
theta_freqs = [4,7];
alpha_freqs = [8,15];
beta_freqs = [18, 30];
low_gamma_freqs = [33 55];
high_gamma_freqs = [70 150];
vns_freqs = [24 26];

%% Initialize Data Arrays (Will save alot of memory)
dat_length_full = zeros(length(test_files),1); % array of length of each test block at full sample rate
dat_length_ds = zeros(length(test_files),1);    % array of lenght of each test block at downsampled sample rate
for k = 1:length(test_files)
    timestamps = hdf5read([root_dir filesep test_files{k}], 'timestamp vector');
    timestamps_ds = resample(timestamps,fs_out,fs_in);
    dat_length_full(k) = length(timestamps);
    dat_length_ds(k) = length(timestamps_ds);
end
% % % onset_inds = [1; (cumsum(ds_length)+1)];
% % % offset_inds = [cumsum(ds_length)];
% % % offset_inds = offset_inds(2:end);
% all_data = zeros(length(electrodes), sum(full_length));






%% Initialize Data:

is_stim_on = false(1, sum(dat_length_ds));
blk_index = zeros(1, sum(dat_length_ds));

if load_envelope
    ecog_theta_env = zeros(length(electrodes), sum(dat_length_ds));
    ecog_alpha_env = zeros(length(electrodes), sum(dat_length_ds));
    ecog_beta_env = zeros(length(electrodes), sum(dat_length_ds));
    ecog_lg_env = zeros(length(electrodes), sum(dat_length_ds));
    ecog_hg_env = zeros(length(electrodes), sum(dat_length_ds));
end
if load_phase
    ecog_theta_ang = zeros(length(electrodes), sum(dat_length_ds));
    ecog_alpha_ang = zeros(length(electrodes), sum(dat_length_ds));
    ecog_beta_ang = zeros(length(electrodes), sum(dat_length_ds));
    ecog_lg_ang = zeros(length(electrodes), sum(dat_length_ds));
    ecog_hg_ang = zeros(length(electrodes), sum(dat_length_ds));
end

if load_vns_band_env
    ecog_vns_env = zeros(length(electrodes), sum(dat_length_ds));
    ecog_vns_raw = zeros(length(electrodes), sum(dat_length_ds));
end

if load_vns_phase
    vns_phase = zeros(1, sum(dat_length_ds));
end

if load_raw
    all_data = zeros(length(electrodes), sum(dat_length_full));
    vns_raw = zeros(1, sum(dat_length_full));
    is_stim_on_full = false(1,sum(dat_length_full));
    cum_data_length = cumsum(dat_length_full);
    file_onset_inds_full = 1+[0; cum_data_length(1:(end-1))]; % starting indecies of each block on downsampled data
    file_offset_ind_full = cum_data_length(1:end); % indicates the index where the block data begins within aggregated data.
end

if load_raw_erps
    raw_erps = [];
end

if load_peak_freqs
    freq_bands = [4,7; 8 15; 18 30; 33 55; 70 150];
    peak_freq_pre_all = [];
    peak_freq_post_all = [];
end







%% Set start/end indices for each block:
cum_data_length = cumsum(dat_length_ds);
file_onset_inds_ds = 1+[0; cum_data_length(1:(end-1))]; % starting indecies of each block on downsampled data
file_offset_ind_ds = cum_data_length(1:end); % indicates the index where the block data begins within aggregated data.

%% Load Data:
for k = 1:length(test_files)
    % Time Data
    disp(k);
    disp(datetime('now'));
    timestamps = hdf5read([root_dir filesep test_files{k}], 'timestamp vector');
    tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));
        %duration=(tsvec(end,5)*60+tsvec(end,6))-(tsvec(1,5)*60+tsvec(1,6));%total time
    duration = (timestamps(end)-timestamps(1));
    fileInfo = h5info([filesep root_dir filesep test_files{k}]);
    [p,q] = rat(1024/round(fileInfo.Datasets(1).Attributes(2).Value)); % does it need an integer value? %samples/sec
    
    % ECoG Data
    data=hdf5read([root_dir filesep test_files{k}],'ECoG Array');
    data = resample(data,p,q);
    if notch_flag
        data(electrodes,:) = applyLineNoiseNotch_VNS_Harmonics(data(electrodes,:),fs_in);
    end
    % Notch60Harmonics
    data(electrodes,:) = apply60HzNotch_filter(data(electrodes,:), fs_in);

    data_z = gdivide(gsubtract(data(electrodes,:),mean(data(electrodes,:),2)),std(data(electrodes,:),[],2));


    %% Generate Analytic Amplitude (Envelope Data)
    if load_envelope
        % Theta
        ecog_block = load_ecog_block(data_z, theta_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_theta_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %alpha
        ecog_block = load_ecog_block(data_z, alpha_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_alpha_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %beta
        ecog_block = load_ecog_block(data_z, beta_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_beta_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %low gamma
        ecog_block = load_ecog_block(data_z, low_gamma_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_lg_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %high gamma
        ecog_block = load_ecog_block(data_z, high_gamma_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_hg_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
    end

    %% Generate Phase Angle Data:
    if load_phase         % Load the bands for phase data
        % Theta
        ecog_block = load_ecog_block(data_z, theta_freqs, false, fs_in, fs_out, dat_length_ds(k));
        ecog_theta_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %alpha
        ecog_block = load_ecog_block(data_z, alpha_freqs, false, fs_in, fs_out, dat_length_ds(k));
        ecog_alpha_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %beta
        ecog_block = load_ecog_block(data_z, beta_freqs, false, fs_in, fs_out, dat_length_ds(k));
        ecog_beta_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %low gamma
        ecog_block = load_ecog_block(data_z, low_gamma_freqs, false, fs_in, fs_out, dat_length_ds(k));
        ecog_lg_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;

        %high gamma
        ecog_block = load_ecog_block(data_z, high_gamma_freqs, false, fs_in, fs_out, dat_length_ds(k));
        ecog_hg_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
    end
    %% Generate Narrow Band VNS stimulation
    if load_vns_band_env
        % get data on all harmonics...
        data_filt = applyLineNoiseBandpass_VNS_Harmonics(data_z,fs_in,25);
        data_raw = data_filt;
        data_intensity = abs(hilbert(data_filt')');
        ecog_block = zeros(size(data_intensity,1), 1-file_onset_inds_ds(k)+file_offset_ind_ds(k));
        ecog_block_raw = zeros(size(data_intensity,1), 1-file_onset_inds_ds(k)+file_offset_ind_ds(k));
        for i = 1:size(data_intensity,1)
            ecog_block(i,:) = resample(data_intensity(i,:),fs_out,fs_in);
        end
        for i = 1:size(data_raw,1)
            ecog_block_raw(i,:) = resample(data_raw(i,:), fs_out, fs_in);
        end
        %ecog_block = load_ecog_block(data_z, vns_freqs, true, fs_in, fs_out, dat_length_ds(k));
        ecog_vns_raw(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block_raw;
        ecog_vns_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
    end
    if load_vns_phase
        vns_blk = data(ekg_ch,:);
        vns_phase_blk = get_stim_phase(vns_blk, fs_in, fs_out, dat_length_ds(k));
        if size(ecog_block,2) ~= size(vns_phase_blk,2)
            a = 1;
        end
%             if (size(ecog_block,2) - size(vns_phase_blk,2)) == 1
%                 vns_phase_ds(:,file_onset_inds_ds(k):(file_offset_ind_ds(k)-1)) = vns_phase_blk;
%             elseif (size(ecog_block,2) - size(vns_phase_blk,2)) == -1
%                 vns_phase_ds(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = vns_phase_blk(1:(end-1));
%             else
%                 warning('downsample_mismatch\n')
%                 size(ecog_block,2) - size(vns_phase_blk,2)
%             end
%         end
        vns_phase_ds(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = vns_phase_blk;
    end

    % Add combined data (bad idea)
    if load_raw
        all_data(:,file_onset_inds_full(k):file_offset_ind_full(k)) = data_z;
        vns_raw(:,file_onset_inds_full(k):file_offset_ind_full(k)) = data(ekg_ch,:)';
    end

    %% Generate Stiulation ON/OFF Data
    ekg_block = data(ekg_ch,:);
    is_stim_on_blk = find_vns_stim_on(ekg_block, fs_in, fs_out, stim_freq, VNS_duty_cycle);
    is_stim_on(file_onset_inds_ds(k):file_offset_ind_ds(k)) = is_stim_on_blk;



    %% Load FileBlockNumber
    blk_index(file_onset_inds_ds(k):file_offset_ind_ds(k)) = k;

    if load_raw
        is_stim_on_blk_full = (ecg_data_vns_smth > on_thresh);
        is_stim_on_full(file_onset_inds_full(k):file_offset_ind_full(k)) = is_stim_on_blk_full;
    end
    %% Load RAW data ERPs:
    if load_raw_erps
        stim_duration_thresh = 2; % minimum time in seconds of stim
        [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_out);               % stim:
        onsets_raw = round(onsets*fs_in/fs_out);

        raw_erp_win = [-10 10];
        [erps_blk, time_axis] = make_vns_erps(data_z, onsets_raw,fs_in, raw_erp_win);
        raw_erps = cat(3,raw_erps, erps_blk);
    end

    if load_peak_freqs
        stim_duration_thresh = 2; % minimum time in seconds of stim
        [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_out);               % stim:
        onsets_raw = round(onsets*fs_in/fs_out);

        raw_erp_win = [-30 30];
        [erps_blk, time_axis] = make_vns_erps(data_z, onsets_raw,fs_in, raw_erp_win);

        for i = 1:size(erps_blk,3);
            dat_pre = squeeze(erps_blk(:,time_axis<0,i));
            dat_post = squeeze(erps_blk(:,time_axis>0,i));
            peak_freqs_pre = zeros(size(erps_blk,1),size(freq_bands,1));
            peak_freqs_post = zeros(size(erps_blk,1),size(freq_bands,1));
            for j = 1:size(erps_blk,1)
                % Peak Freq Pre:
                [pxx,f_axis] = pwelch(dat_pre(j,:),[],[],[],fs_in);
                 for h = 1:size(freq_bands,1)
                    in_range = (f_axis >= freq_bands(h,1)) & (f_axis <= freq_bands(h,2));
                    [~, max_ind] = max(pxx(in_range));
                    f_axis_range = f_axis(in_range);
                    peak_freqs_pre(j,h) = f_axis_range(max_ind);
                 end

                 % Peak Freq Post
                 [pxx,f_axis] = pwelch(dat_post(j,:),[],[],[],fs_in);
                 for h = 1:size(freq_bands,1)
                    in_range = (f_axis >= freq_bands(h,1)) & (f_axis <= freq_bands(h,2));
                    [~, max_ind] = max(pxx(in_range));
                    f_axis_range = f_axis(in_range);
                    peak_freqs_post(j,h) = f_axis_range(max_ind);
                 end
            end
            peak_freq_pre_all = cat(3, peak_freq_pre_all, peak_freqs_pre);
            peak_freq_post_all = cat(3, peak_freq_post_all, peak_freqs_post);
        end
    end


    a = 1;
end



%% After processes:
%fp = round(median(fp));
is_stim_on = logical(is_stim_on);



%% get event onsets:
stim_duration_thresh = 2; % minimum time in seconds of stim
[onsets, offsets] = find_boolean_on(is_stim_on, stim_duration_thresh*fs_out);               % stim:
[onsets_nostim, offsets_nostim] = find_boolean_on(~is_stim_on, stim_duration_thresh*fs_out); % non-stim
if load_raw
    [onsets_full, offsets_full] = find_boolean_on(is_stim_on_full, stim_duration_thresh*fs_in);               % stim:
    [onsets_nostima_full, offsets_nostim_full] = find_boolean_on(~is_stim_on_full, stim_duration_thresh*fs_in); % non-stim
end





%% Set Strucutred Array
VNS_dat.dat_dir = root_dir; % the directory containing the raw data
VNS_dat.sampFreq = fs_out; % sample frequency
VNS_dat.stim_onsets_inds = onsets;
VNS_dat.stim_offset_inds = offsets;
VNS_dat.electrode_order = electrode_order;
VNS_dat.is_stim_on = is_stim_on;
VNS_dat.blk_index = blk_index; % index of list of blocks that gave rise to current data.
VNS_dat.anatomy_elecs = anatomy_elecs;
if load_envelope
    VNS_dat.ecog_theta_env = ecog_theta_env;
    VNS_dat.ecog_alpha_env = ecog_alpha_env;
    VNS_dat.ecog_beta_env = ecog_beta_env;
    VNS_dat.ecog_lg_env = ecog_lg_env;
    VNS_dat.ecog_hg_env = ecog_hg_env;
end
if load_phase
    VNS_dat.ecog_theta_ang = ecog_theta_ang;
    VNS_dat.ecog_alpha_ang = ecog_alpha_ang;
    VNS_dat.ecog_beta_ang = ecog_beta_ang;
    VNS_dat.ecog_lg_ang = ecog_lg_ang;
    VNS_dat.ecog_hg_ang = ecog_hg_ang;
end
if load_vns_band_env
    VNS_dat.ecog_vns_env = ecog_vns_env;
    VNS_dat.ecog_vns_raw = ecog_vns_raw;
end
if load_vns_phase
    VNS_dat.vns_phase_ds = vns_phase_ds
end
if load_raw
    VNS_dat.all_data = all_data;
    VNS_dat.sampFreq_raw = fs_in;
    VNS_dat.vns_raw = vns_raw;
    VNS_dat.is_stim_on_full = is_stim_on_full;
    VNS_dat.onsets_raw = onsets_raw;
end
if load_raw_erps
    VNS_dat.raw_erps = raw_erps;
end


if load_peak_freqs
    VNS_dat.peak_freq_pre_all = peak_freq_pre_all;
    VNS_dat.peak_freq_post_all = peak_freq_post_all;
end



%% Delete Bad Channels:
bad_chan_list = [];
%bad_chan_list = [1:19 41 43 60 71];
%bad_chan_list = [1:19 41 43 60 70 71]; % used in first set of recordings
bad_chan_list = [1:19 41 43 60 70 71 89 120]; % used in 9/4 recordings
good_channels = true(length(anatomy_elecs),1); % boolean reports true for bad channels
good_channels(bad_chan_list) = false;
VNS_dat.good_channels = good_channels;


%% DELETE BAD TRIALS;
%bad_trials_list = [16, 27, 65, 96];
% % % bad_trials_list = [16, 27, 65, 96]; % list of stimulation onsets with bad stuff happening:
% % % if ~load_full_dir
% % %     bad_trials_list(bad_trials_list > length(VNS_dat.stim_onsets_inds)) = [];
% % % end
% % % VNS_dat.stim_onsets_inds(bad_trials_list) = [];
% % % VNS_dat.stim_offset_inds(bad_trials_list) = [];



end
