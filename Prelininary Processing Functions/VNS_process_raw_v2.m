function VNS_dat = VNS_process_raw_v2(subj, rootdir, brain_dir, ekg_ch, electrodes, varargin)

%% Inputs:
%
% subj - the code for the patient, e.g., EC131
%
% rootdir - the directory where raw data is stored
% (Ex: '/Users/changlab/Documents/data/VNS/EC131/EC131_05104e39')
%
% ekg_ch - CH containing EKG data can vary per subject (Ex: [137 138])
%
% electrodes - list of relevant electrodes (Ex: 1:136 or [44, 45])
%
%% Flag Inputs:
%
% 'load_range'      - default [] (load all h5 files in rootdir)
% 'fsDs'            - default 1024 (don't downsample)
% 'recType'         - default 'clinical'
% 'ERP_times'       - default [-10 10] (seconds around VNS onset)
% 'notch_60Hz_flag' - default true
% 'notch_VNS_flag'  - default true
% 'CARflag'         - default false
% 'stim_freq'       - default 25 (frequency of VNS stimulation)
% 'VNS_duty_cycle'  - default [2,26,2] (ramp,duration,ramp)
% 'save_raw_erps'   - default false

%% Output:
% VNS_dat - the standard data structure for subsequent VNS analysis. By
% default contains timing information on Stimulations

%% Load list of h5 filenames

% which h5 files to use
if ~isempty(find(strcmpi(varargin,'load_range')));
    load_range = varargin{find(strcmpi(varargin,'load_range'))+1};
else
    load_range = [];
end

filenames = dir(rootdir);
filenames(1:2) = []; % on windows, only two files

test_files = {filenames.name};
for k = 1:length(test_files)
    tmp = test_files{k};
    is_h5(k) = strcmpi(tmp((end-2):end),'.h5');
end
test_files(~is_h5) = [];

if ~isempty(load_range)
    test_files = test_files(load_range);
end

fileInfo = h5info([filesep rootdir filesep test_files{1}]);
fs = fileInfo.Datasets(1).Attributes(2).Value;

%% Load Flags

% get sampling rate info
if ~isempty(find(strcmpi(varargin,'fsDs')));
    fsDs = varargin{find(strcmpi(varargin,'fsDs'))+1};
else
    fsDs = 1024; % don't downsample by default
end

% get electrode array type
if ~isempty(find(strcmpi(varargin,'recType')));
    recType = varargin{find(strcmpi(varargin,'recType'))+1};
else
    recType = 'clinical';
end

% ERP window to use
if ~isempty(find(strcmpi(varargin,'ERP_times')));
    ERP_times = varargin{find(strcmpi(varargin,'ERP_times'))+1};
else
    ERP_times = [-10 10]; % default
end

% whether to notch 60Hz + harmonics
if ~isempty(find(strcmpi(varargin,'notch_60Hz_flag')));
    notch_60Hz_flag = varargin{find(strcmpi(varargin,'notch_60Hz_flag'))+1};
else
    notch_60Hz_flag = true; % notch 60Hz by default
end

% notch VNS artifact: notch data to remove the VNS harmonics
if ~isempty(find(strcmpi(varargin,'notch_VNS_flag')))
    notch_VNS_flag = varargin{find(strcmpi(varargin,'notch_VNS_flag'))+1};
else
    notch_VNS_flag = true; % notch VNS artifact by default
end

% whether to CAR the data
if ~isempty(find(strcmpi(varargin,'CARflag')));
    CARflag = varargin{find(strcmpi(varargin,'CARflag'))+1};
else
    CARflag = false;
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
    VNS_duty_cycle = [2,26,2];
end

% generate erps of the RAW data about onsets
if ~isempty(find(strcmpi(varargin,'save_raw_erps')))
    save_raw_erps = varargin{find(strcmpi(varargin,'save_raw_erps'))+1};
    raw_erps = [];
else
    save_raw_erps = false;
end

% whether to save output
if ~isempty(find(strcmpi(varargin,'save_preproc_flag')))
    save_preproc_flag = varargin{find(strcmpi(varargin,'save_preproc_flag'))+1};
else
    save_preproc_flag = false;
end

% out directory
if ~isempty(find(strcmpi(varargin,'out_dir')))
    out_dir = varargin{find(strcmpi(varargin,'out_dir'))+1};
end

% outfile name (and ask if you want to overwrite existing file)
if ~isempty(find(strcmpi(varargin,'outfile_name')))
    outfile_name = varargin{find(strcmpi(varargin,'outfile_name'))+1};
end
if exist([out_dir filesep outfile_name '.mat'])
    cont_flag = input('Outfile already exists. Overwrite? (y/n) ','s');
    if strcmpi(cont_flag,'n')
        return;
    elseif strcmpi(cont_flag,'y')
        fprintf('OVERWRITING EXISTING .mat FILE!\n');
    end
end

%% Process Data Parameters:
theta_freqs = [4,7];
alpha_freqs = [8,15];
beta_freqs = [18, 30];
low_gamma_freqs = [33 55];
high_gamma_freqs = [70 150];

%% Load anatomy file:
anatomy_file = [brain_dir filesep subj '/elecs/' recType '_elecs_all.mat'];
if exist(anatomy_file) == 2
    load([brain_dir filesep subj '/elecs/' recType '_elecs_all.mat']);
    anatomy_elecs = anatomy(electrodes,4);
else
    anatomy_elecs = strsplit(num2str(1:length(electrodes)));
end

%% Find Timing Data and get dsFactor

timestamps = hdf5read([rootdir filesep test_files{1}], 'timestamp vector');
duration = timestamps(end) - timestamps(1);
[p,q] = rat(fsDs/fs);

%% Initialize Data Arrays (Will save alot of memory)
dat_length_full = zeros(length(test_files),1); % array of length of each test block at full sample rate
dat_length_ds = zeros(length(test_files),1);    % array of length of each test block at downsampled sample rate
for k = 1:length(test_files)
    timestamps = hdf5read([rootdir filesep test_files{k}], 'timestamp vector');
    timestamps_ds = resample(timestamps,p,q);
    dat_length_full(k) = length(timestamps);
    dat_length_ds(k) = length(timestamps_ds);
end

%% Set start/end indices for each block:
cum_data_length = cumsum(dat_length_ds);
file_onset_inds_ds = 1+[0; cum_data_length(1:(end-1))]; % starting indecies of each block on downsampled data
file_offset_ind_ds = cum_data_length(1:end); % indicates the index where the block data begins within aggregated data.
file_info = [];

%% Load Data:
for k = 1:length(test_files)
    fprintf('\n');
    % Time Data
    fprintf('Processing block [%d] of [%d]\n',k,length(test_files));
    timestamps = hdf5read([rootdir filesep test_files{k}], 'timestamp vector');
    tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));
    duration = (timestamps(end)-timestamps(1));
    fileInfo = h5info([filesep rootdir filesep test_files{k}]);
    
    % ECoG Data
    fprintf('Loading h5 file....');
    tic;
    data_full = h5read([rootdir filesep test_files{k}],'/ECoG Array');
    toc
    
    % create NaN matrix with correct size and then resample data
    tmp = resample(data_full(1,:),p,q);
    data = NaN(length(electrodes),size(tmp,2));
    clear tmp;
    textprogressbar('Resampling data on electrodes: ');
    for i = 1:length(electrodes)
        textprogressbar((i/length(electrodes))*100);
        data(electrodes(i),:) = resample(data_full(electrodes(i),:),p,q);
    end
    fprintf('\n');
    clear data_full;
    
    % Notch 60Hz + Harmonics
    if notch_60Hz_flag
%         fprintf('Applying line noise filters\n');
        textprogressbar('Applying line noise filters: ');
        for i = 1:length(electrodes)
%             fprintf('Electrode [%d] of [%d]\n',i,length(electrodes));
            textprogressbar((i/length(electrodes))*100);
            data(electrodes(i),:) = apply60HzNotch_filter(data(electrodes(i),:), fsDs);
        end
        fprintf('\n');
    end
    
    % Notch VNS artifact + Harmonics
    if notch_VNS_flag
%         fprintf('Applying VNS notch filters\n');
        textprogressbar('Applying VNS notch filters: ');
        for i = 1:length(electrodes)
%             fprintf('Electrode [%d] of [%d]\n',i,length(electrodes));
            textprogressbar((i/length(electrodes))*100);
            if ~ismember(electrodes(i),ekg_ch)
                data(electrodes(i),:) = applyLineNoiseNotch_VNS_Harmonics(data(electrodes(i),:),fsDs,'stim_freq',stim_freq);
            end
        end
        fprintf('\n');
    end
    
    if CARflag
        fprintf('Applying CAR\n');
        data = data - nanmean(data,1);
    end
    
%     data_z = data;
    %     data_z = gdivide(gsubtract(data(electrodes,:),mean(data(electrodes,:),2)),std(data(electrodes,:),[],2));
    
    %% FIND VNS events
    if length(ekg_ch) > 1
        ekg_block = data(electrodes(find(electrodes == ekg_ch(1))),:) - data(electrodes(find(electrodes == ekg_ch(2))),:);
    else
        ekg_block = data(electrodes(find(electrodes == ekg_ch)),:);
    end
    [is_stim_on_blk,onsets{k},offsets{k}] = find_vns_stim_on(ekg_block,fsDs,...
        'stim_freq',stim_freq,...
        'VNS_duty_cycle',VNS_duty_cycle,...
        'ERP_times',ERP_times,...
        'debug_flag',1,...
        'convolve_ideal_stim_flag',0);
    is_stim_on(file_onset_inds_ds(k):file_offset_ind_ds(k)) = is_stim_on_blk;
    
    if length(onsets{k}) > 0
        file_info = [file_info ; cellstr(repmat(test_files{k}(7:14),length(onsets{k}),1))];
    end
        
    %% Load FileBlockNumber
    blk_index(file_onset_inds_ds(k):file_offset_ind_ds(k)) = k;
    
    %% Load RAW data ERPs:
    if save_raw_erps
        [erps_blk, time_axis] = make_vns_erps(data,onsets{k},fsDs,ERP_times);
        raw_erps = cat(3,raw_erps, erps_blk);
    end

end

%% CREATE VNS_dat STRUCTURE
VNS_dat.dat_dir = rootdir; % the directory containing the raw data
VNS_dat.sampFreq = fsDs; % sample frequency
VNS_dat.stim_onsets_inds = onsets;
VNS_dat.stim_offset_inds = offsets;
VNS_dat.file_info = file_info;
VNS_dat.is_stim_on = is_stim_on;
VNS_dat.blk_index = blk_index; % index of list of blocks that gave rise to current data.
VNS_dat.anatomy_elecs = anatomy_elecs;
VNS_dat.EKG_ch = ekg_ch;
VNS_dat.time_axis = time_axis;
if save_raw_erps
    VNS_dat.raw_erps = raw_erps;
end
if save_preproc_flag
    fprintf('Saving outfile....\n');
    save([out_dir filesep outfile_name '.mat'],'VNS_dat','-v7.3');
end




%%
%%%%

% 
% 
% %% Load Flags:
% % Load ECoG Analytic Amplidutes of Theta - High Gamma bands:
% if ~isempty(find(strcmpi(varargin,'load_env')))
%     load_envelope = varargin{find(strcmpi(varargin,'load_env'))+1};
% else
%     load_envelope = true;       % Load the bands for envelope
% end
% 
% % Load ECoG Phase Angle of Theta - High Gamma bands:
% if ~isempty(find(strcmpi(varargin,'load_phase')))
%     load_phase = varargin{find(strcmpi(varargin,'load_phase'))+1};
% else
%     load_phase = false;          % Load the bands for phase data
% end
% 
% load_vns_phase = false;
% 
% % Load the intensity of ecog bandpassed for VNS harmonics
% if ~isempty(find(strcmpi(varargin,'load_vns_band')))
%     load_vns_band_env = varargin{find(strcmpi(varargin,'load_vns_band'))+1};
% else
%     load_vns_band_env = false;
% end
% 
% % load raw data:
% % Load the combined Raw data for all blocks (HUGE)
% if ~isempty(find(strcmpi(varargin,'load_RAW')))
%     load_raw = varargin{find(strcmpi(varargin,'load_raw'))+1};
% else
%     load_raw = false;
% end
%
% % load peak frequency from the raw:
% if ~isempty(find(strcmpi(varargin,'load_peak_freqs')));
%     load_peak_freqs = varargin{find(strcmpi(varargin,'load_peak_freqs'))+1};
% else
%     load_peak_freqs = false;
% end
% 
% blank_singular_anatomy = false; % ignore the names of anatomical regions that are found only once
% load_stim_freq_raw = false; % load the stimulation frequency harmonics
% 
% %% Initialize Data:
% 
% is_stim_on = false(1, sum(dat_length_ds));
% blk_index = zeros(1, sum(dat_length_ds));
% 
% if load_envelope
%     ecog_theta_env = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_alpha_env = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_beta_env = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_lg_env = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_hg_env = zeros(length(electrodes), sum(dat_length_ds));
% end
% if load_phase
%     ecog_theta_ang = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_alpha_ang = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_beta_ang = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_lg_ang = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_hg_ang = zeros(length(electrodes), sum(dat_length_ds));
% end
% 
% if load_vns_band_env
%     ecog_vns_env = zeros(length(electrodes), sum(dat_length_ds));
%     ecog_vns_raw = zeros(length(electrodes), sum(dat_length_ds));
% end
% 
% if load_vns_phase
%     vns_phase = zeros(1, sum(dat_length_ds));
% end
% 
% if load_raw
%     all_data = zeros(length(electrodes), sum(dat_length_full));
%     vns_raw = zeros(1, sum(dat_length_full));
%     is_stim_on_full = false(1,sum(dat_length_full));
%     cum_data_length = cumsum(dat_length_full);
%     file_onset_inds_full = 1+[0; cum_data_length(1:(end-1))]; % starting indecies of each block on downsampled data
%     file_offset_ind_full = cum_data_length(1:end); % indicates the index where the block data begins within aggregated data.
% end
% 
% if load_raw_erps
%     raw_erps = [];
% end
% 
% if load_peak_freqs
%     freq_bands = [4,7; 8 15; 18 30; 33 55; 70 150];
%     peak_freq_pre_all = [];
%     peak_freq_post_all = [];
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %% Load Data:
% for k = 1:length(test_files)
%     % Time Data
%     fprintf('Processing block [%d] of [%d]\n\n',k,length(test_files));
%     disp(datetime('now'));
%     timestamps = hdf5read([rootdir filesep test_files{k}], 'timestamp vector');
%     tsvec = datevec((timestamps - (7*3600))/86400 + datenum(1970,1,1));
%         %duration=(tsvec(end,5)*60+tsvec(end,6))-(tsvec(1,5)*60+tsvec(1,6));%total time
%     duration = (timestamps(end)-timestamps(1));
%     fileInfo = h5info([filesep rootdir filesep test_files{k}]);
%     [p,q] = rat(1024/fileInfo.Datasets(1).Attributes(2).Value); % MKL FIXED 5/17/17
%     
%     % ECoG Data - MKL CHANGED TO RESAMPLE IN LOOP
%     fprintf('Loading h5 file....');
%     tic;
%     data_full = h5read([rootdir filesep test_files{k}],'/ECoG Array');
%     toc
%     fprintf('Done\n');
%     tmp = resample(data_full(1,:),p,q);
%     data = NaN(size(data_full,1),size(tmp,2));
%     clear tmp;
%     for i = 1:length(electrodes)
%         fprintf('Resampling electrode [%d] of [%d]\n',i,length(electrodes));
%         data(electrodes(i),:) = resample(data_full(electrodes(i),:),p,q);
%     end
%     clear data_full;
% %     data=hdf5read([root_dir filesep test_files{k}],'ECoG Array');
% %     data = resample(data,p,q); % MKL: THIS DOESN'T WORK BECAUSE COLUMN-WISE
%     
%     if notch_flag % Notch VNS Harmonics - MKL CHANGED TO FILTER BY ELECTRODE
%         fprintf('Applying VNS notch filters\n');
%         for i = 1:length(electrodes)
%             fprintf('Electrode [%d] of [%d]\n',i,length(electrodes));
%             data(electrodes(i),:) = applyLineNoiseNotch_VNS_Harmonics(data(electrodes(i),:),fs_in);
%         end
% %         data(electrodes,:) = applyLineNoiseNotch_VNS_Harmonics(data(electrodes,:),fs_in);
%     end
%     
%     % Notch60Harmonics - MKL CHANGED TfO FILTER BY ELECTRODE
%     fprintf('Applying line noise filters\n');
%     for i = 1:length(electrodes)
%         fprintf('Electrode [%d] of [%d]\n',i,length(electrodes));
%         data(electrodes(i),:) = apply60HzNotch_filter(data(electrodes(i),:), fs_in);
%     end
% 
%     if CARflag
%         fprintf('Applying CAR\n');
%         data = data - nanmean(data,1);
%     end
%     
% %     data_z = data;
% %     data_z = gdivide(gsubtract(data(electrodes,:),mean(data(electrodes,:),2)),std(data(electrodes,:),[],2));
% 
%     %% FIND VNS events
%     if length(ekg_ch) > 1
%         ekg_block = data(electrodes(find(electrodes == ekg_ch(1))),:) - data(electrodes(find(electrodes == ekg_ch(2))),:);
%     else
%         ekg_block = data(electrodes(find(electrodes == ekg_ch)),:);
%     end
%     is_stim_on_blk = find_vns_stim_on(ekg_block, fs_in, fs_out, stim_freq, VNS_duty_cycle); % MKL NO DOWNSAMPLE
%     is_stim_on(file_onset_inds_ds(k):file_offset_ind_ds(k)) = is_stim_on_blk;
% 
% 
% 
%     %% Load FileBlockNumber
%     blk_index(file_onset_inds_ds(k):file_offset_ind_ds(k)) = k;
% 
%     if load_raw
%         is_stim_on_blk_full = (ecg_data_vns_smth > on_thresh);
%         is_stim_on_full(file_onset_inds_full(k):file_offset_ind_full(k)) = is_stim_on_blk_full;
%     end
%     %% Load RAW data ERPs:
%     if load_raw_erps
%         stim_duration_thresh = 2; % minimum time in seconds of stim
%         [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_in);               % stim:
%         onsets_raw = onsets;
% %         [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_out);               % stim:
% %         onsets_raw = round(onsets*fs_in/fs_out);
% 
%         raw_erp_win = [-10 10];
%         [erps_blk, time_axis] = make_vns_erps(data_z, onsets_raw,fs_in, raw_erp_win);
%         raw_erps = cat(3,raw_erps, erps_blk);
%     end
%     
% 
% 
% 
%     %% Generate Analytic Amplitude (Envelope Data)
%     if load_envelope
%         % Theta
%         fprintf('Performing hilbert transform for theta band\n');
%         ecog_block = load_ecog_block(data_z(electrodes,:), theta_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_theta_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %alpha
%         fprintf('Performing hilbert transform for alpha band\n');
%         ecog_block = load_ecog_block(data_z(electrodes,:), alpha_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_alpha_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %beta
%         fprintf('Performing hilbert transform for beta band\n');
%         ecog_block = load_ecog_block(data_z(electrodes,:), beta_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_beta_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %low gamma
%         fprintf('Performing hilbert transform for low gamma band\n');
%         ecog_block = load_ecog_block(data_z(electrodes,:), low_gamma_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_lg_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %high gamma
%         fprintf('Performing hilbert transform for high gamma band\n');
%         ecog_block = load_ecog_block(data_z(electrodes,:), high_gamma_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_hg_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
%     end
% 
%     %% Generate Phase Angle Data:
%     if load_phase         % Load the bands for phase data
%         % Theta
%         ecog_block = load_ecog_block(data_z, theta_freqs, false, fs_in, fs_out, dat_length_ds(k));
%         ecog_theta_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %alpha
%         ecog_block = load_ecog_block(data_z, alpha_freqs, false, fs_in, fs_out, dat_length_ds(k));
%         ecog_alpha_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %beta
%         ecog_block = load_ecog_block(data_z, beta_freqs, false, fs_in, fs_out, dat_length_ds(k));
%         ecog_beta_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %low gamma
%         ecog_block = load_ecog_block(data_z, low_gamma_freqs, false, fs_in, fs_out, dat_length_ds(k));
%         ecog_lg_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
% 
%         %high gamma
%         ecog_block = load_ecog_block(data_z, high_gamma_freqs, false, fs_in, fs_out, dat_length_ds(k));
%         ecog_hg_ang(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
%     end
%     %% Generate Narrow Band VNS stimulation
%     if load_vns_band_env
%         % get data on all harmonics...
%         data_filt = applyLineNoiseBandpass_VNS_Harmonics(data_z,fs_in,25);
%         data_raw = data_filt;
%         data_intensity = abs(hilbert(data_filt')');
%         ecog_block = zeros(size(data_intensity,1), 1-file_onset_inds_ds(k)+file_offset_ind_ds(k));
%         ecog_block_raw = zeros(size(data_intensity,1), 1-file_onset_inds_ds(k)+file_offset_ind_ds(k));
%         for i = 1:size(data_intensity,1)
%             ecog_block(i,:) = resample(data_intensity(i,:),fs_out,fs_in);
%         end
%         for i = 1:size(data_raw,1)
%             ecog_block_raw(i,:) = resample(data_raw(i,:), fs_out, fs_in);
%         end
%         %ecog_block = load_ecog_block(data_z, vns_freqs, true, fs_in, fs_out, dat_length_ds(k));
%         ecog_vns_raw(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block_raw;
%         ecog_vns_env(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = ecog_block;
%     end
%     if load_vns_phase
%         vns_blk = data(ekg_ch,:);
%         vns_phase_blk = get_stim_phase(vns_blk, fs_in, fs_out, dat_length_ds(k));
%         if size(ecog_block,2) ~= size(vns_phase_blk,2)
%             a = 1;
%         end
% %             if (size(ecog_block,2) - size(vns_phase_blk,2)) == 1
% %                 vns_phase_ds(:,file_onset_inds_ds(k):(file_offset_ind_ds(k)-1)) = vns_phase_blk;
% %             elseif (size(ecog_block,2) - size(vns_phase_blk,2)) == -1
% %                 vns_phase_ds(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = vns_phase_blk(1:(end-1));
% %             else
% %                 warning('downsample_mismatch\n')
% %                 size(ecog_block,2) - size(vns_phase_blk,2)
% %             end
% %         end
%         vns_phase_ds(:,file_onset_inds_ds(k):file_offset_ind_ds(k)) = vns_phase_blk;
%     end
% 
%     % Add combined data (bad idea)
%     if load_raw
%         all_data(:,file_onset_inds_full(k):file_offset_ind_full(k)) = data_z;
%         vns_raw(:,file_onset_inds_full(k):file_offset_ind_full(k)) = data(ekg_ch,:)';
%     end
% 
%     %% Generate Stimulation ON/OFF Data
%     if length(ekg_ch) > 1
%         ekg_block = data(electrodes(find(electrodes == ekg_ch(1))),:) - data(electrodes(find(electrodes == ekg_ch(2))),:);
%     else
%         ekg_block = data(electrodes(find(electrodes == ekg_ch)),:);
%     end
%     is_stim_on_blk = find_vns_stim_on(ekg_block, fs_in, fs_out, stim_freq, VNS_duty_cycle); % MKL NO DOWNSAMPLE
% %     is_stim_on_blk = find_vns_stim_on(ekg_block, fs_in, fs_out, stim_freq, VNS_duty_cycle);
%     is_stim_on(file_onset_inds_ds(k):file_offset_ind_ds(k)) = is_stim_on_blk;
% 
% 
% 
%     %% Load FileBlockNumber
%     blk_index(file_onset_inds_ds(k):file_offset_ind_ds(k)) = k;
% 
%     if load_raw
%         is_stim_on_blk_full = (ecg_data_vns_smth > on_thresh);
%         is_stim_on_full(file_onset_inds_full(k):file_offset_ind_full(k)) = is_stim_on_blk_full;
%     end
%     %% Load RAW data ERPs:
%     if load_raw_erps
%         stim_duration_thresh = 2; % minimum time in seconds of stim
%         [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_in);               % stim:
%         onsets_raw = onsets;
% %         [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_out);               % stim:
% %         onsets_raw = round(onsets*fs_in/fs_out);
% 
%         raw_erp_win = [-10 10];
%         [erps_blk, time_axis] = make_vns_erps(data_z, onsets_raw,fs_in, raw_erp_win);
%         raw_erps = cat(3,raw_erps, erps_blk);
%     end
% 
%     if load_peak_freqs
%         stim_duration_thresh = 2; % minimum time in seconds of stim
%         [onsets, offsets] = find_boolean_on(is_stim_on_blk, stim_duration_thresh*fs_out);               % stim:
%         onsets_raw = round(onsets*fs_in/fs_out);
% 
%         raw_erp_win = [-30 30];
%         [erps_blk, time_axis] = make_vns_erps(data_z, onsets_raw,fs_in, raw_erp_win);
% 
%         for i = 1:size(erps_blk,3);
%             dat_pre = squeeze(erps_blk(:,time_axis<0,i));
%             dat_post = squeeze(erps_blk(:,time_axis>0,i));
%             peak_freqs_pre = zeros(size(erps_blk,1),size(freq_bands,1));
%             peak_freqs_post = zeros(size(erps_blk,1),size(freq_bands,1));
%             for j = 1:size(erps_blk,1)
%                 % Peak Freq Pre:
%                 [pxx,f_axis] = pwelch(dat_pre(j,:),[],[],[],fs_in);
%                  for h = 1:size(freq_bands,1)
%                     in_range = (f_axis >= freq_bands(h,1)) & (f_axis <= freq_bands(h,2));
%                     [~, max_ind] = max(pxx(in_range));
%                     f_axis_range = f_axis(in_range);
%                     peak_freqs_pre(j,h) = f_axis_range(max_ind);
%                  end
% 
%                  % Peak Freq Post
%                  [pxx,f_axis] = pwelch(dat_post(j,:),[],[],[],fs_in);
%                  for h = 1:size(freq_bands,1)
%                     in_range = (f_axis >= freq_bands(h,1)) & (f_axis <= freq_bands(h,2));
%                     [~, max_ind] = max(pxx(in_range));
%                     f_axis_range = f_axis(in_range);
%                     peak_freqs_post(j,h) = f_axis_range(max_ind);
%                  end
%             end
%             peak_freq_pre_all = cat(3, peak_freq_pre_all, peak_freqs_pre);
%             peak_freq_post_all = cat(3, peak_freq_post_all, peak_freqs_post);
%         end
%     end
% 
% 
%     a = 1;
% end
% 
% 
% 
% %% After processes:
% %fp = round(median(fp));
% is_stim_on = logical(is_stim_on);
% 
% 
% 
% %% get event onsets:
% stim_duration_thresh = 2; % minimum time in seconds of stim
% [onsets, offsets] = find_boolean_on(is_stim_on, stim_duration_thresh*fs_in);               % stim:
% [onsets_nostim, offsets_nostim] = find_boolean_on(~is_stim_on, stim_duration_thresh*fs_in); % non-stim
% % [onsets, offsets] = find_boolean_on(is_stim_on, stim_duration_thresh*fs_out);               % stim:
% % [onsets_nostim, offsets_nostim] = find_boolean_on(~is_stim_on, stim_duration_thresh*fs_out); % non-stim
% if load_raw
%     [onsets_full, offsets_full] = find_boolean_on(is_stim_on_full, stim_duration_thresh*fs_in);               % stim:
%     [onsets_nostima_full, offsets_nostim_full] = find_boolean_on(~is_stim_on_full, stim_duration_thresh*fs_in); % non-stim
% end
% 
% 
% 
% 
% 
% %% Set Strucutred Array
% VNS_dat.dat_dir = rootdir; % the directory containing the raw data
% VNS_dat.sampFreq = fs_out; % sample frequency
% VNS_dat.stim_onsets_inds = onsets;
% VNS_dat.stim_offset_inds = offsets;
% VNS_dat.electrode_order = electrode_order;
% VNS_dat.is_stim_on = is_stim_on;
% VNS_dat.blk_index = blk_index; % index of list of blocks that gave rise to current data.
% VNS_dat.anatomy_elecs = anatomy_elecs;
% if load_envelope
%     VNS_dat.ecog_theta_env = ecog_theta_env;
%     VNS_dat.ecog_alpha_env = ecog_alpha_env;
%     VNS_dat.ecog_beta_env = ecog_beta_env;
%     VNS_dat.ecog_lg_env = ecog_lg_env;
%     VNS_dat.ecog_hg_env = ecog_hg_env;
% end
% if load_phase
%     VNS_dat.ecog_theta_ang = ecog_theta_ang;
%     VNS_dat.ecog_alpha_ang = ecog_alpha_ang;
%     VNS_dat.ecog_beta_ang = ecog_beta_ang;
%     VNS_dat.ecog_lg_ang = ecog_lg_ang;
%     VNS_dat.ecog_hg_ang = ecog_hg_ang;
% end
% if load_vns_band_env
%     VNS_dat.ecog_vns_env = ecog_vns_env;
%     VNS_dat.ecog_vns_raw = ecog_vns_raw;
% end
% if load_vns_phase
%     VNS_dat.vns_phase_ds = vns_phase_ds
% end
% if load_raw
%     VNS_dat.all_data = all_data;
%     VNS_dat.sampFreq_raw = fs_in;
%     VNS_dat.vns_raw = vns_raw;
%     VNS_dat.is_stim_on_full = is_stim_on_full;
%     VNS_dat.onsets_raw = onsets_raw;
% end
% if load_raw_erps
%     VNS_dat.raw_erps = raw_erps;
% end
% 
% 
% if load_peak_freqs
%     VNS_dat.peak_freq_pre_all = peak_freq_pre_all;
%     VNS_dat.peak_freq_post_all = peak_freq_post_all;
% end
% 
% 
% 
% %% Delete Bad Channels:
% % bad_chan_list = [];
% % %bad_chan_list = [1:19 41 43 60 71];
% % %bad_chan_list = [1:19 41 43 60 70 71]; % used in first set of recordings
% % bad_chan_list = [1:19 41 43 60 70 71 89 120]; % used in 9/4 recordings
% % good_channels = true(length(anatomy_elecs),1); % boolean reports true for bad channels
% % good_channels(bad_chan_list) = false;
% % VNS_dat.good_channels = good_channels;
% 
% 
% %% DELETE BAD TRIALS;
% %bad_trials_list = [16, 27, 65, 96];
% % % % bad_trials_list = [16, 27, 65, 96]; % list of stimulation onsets with bad stuff happening:
% % % % if ~load_full_dir
% % % %     bad_trials_list(bad_trials_list > length(VNS_dat.stim_onsets_inds)) = [];
% % % % end
% % % % VNS_dat.stim_onsets_inds(bad_trials_list) = [];
% % % % VNS_dat.stim_offset_inds(bad_trials_list) = [];
% 
% 
% 
% end
