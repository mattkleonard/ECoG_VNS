subj = 'EC131';
rootdir = '/Users/mattleonard/Documents/Research/data';

recType = 'clinical';       % 'clinical' or 'TDT'
ch = 1:136;                 % indices of ECoG channels to use
EKG_ch = [137 138];         % indices of EKG channels
load_range = [];            % if empty, load all files in h5, otherwise load specified files
notch_60Hz_flag = 1;        % whether to notch 60Hz + harmonics
notch_VNS_flag = 1;         % whether to notch VNS artifact
CARflag = 0;                % whether to CAR the data
fsDs = 400;                 % desired downsampled fs
ERP_times = [-10 10];       % window for ERPs

VNS_duty_cycle = [2,26,2];  % timing of VNS duty cycle (ramp,duration,ramp)

%% SET UP DIRECTORIES

dat_dir = [rootdir filesep subj filesep 'VNS/h5'];
brain_dir = [rootdir filesep '..' filesep 'pia/data_store2/imaging/subjects'];
out_dir = [rootdir filesep subj filesep 'VNS'];

load([brain_dir '/' subj '/elecs/' recType '_elecs_all.mat']);

ch = [ch EKG_ch];


%% RUN VNS_process_raw_v2

VNS_dat = VNS_process_raw_v2(subj, dat_dir, brain_dir, EKG_ch, ch,...
    'load_range',load_range,...
    'fsDs',fsDs,...
    'ERP_times',ERP_times',...
    'recType',recType,...
    'VNS_duty_cycle', VNS_duty_cycle,...
    'notch_60Hz_flag',true,...
    'notch_VNS_flag',true,...
    'CARflag',CARflag,...
    'save_raw_erps',true,...
    'load_env',false,'notch_data',false);