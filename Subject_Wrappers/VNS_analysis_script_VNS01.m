% VNS process data script:
root_dir = '/Users/changlab/Documents/data/VNS/VNS_Initial_Test';
dat_dir = [root_dir '/data'];
brain_dir = [root_dir '/BrainPlot'];

%% Find bad channels:
VNS_raw = VNS_process_raw(dat_dir, 127 ,1:120, [],'load_RAW_erps',true,'load_env',false,'load_full_dir',false,'notch_data',false,'VNS_duty_cycle', [21 3]);
[is_bad_chan] = plot_spectra_debug(VNS_raw, 1000, 0.6);
clear VNS_raw
%% Find bad times:
VNS_dat = VNS_process_raw(dat_dir, 127 ,1:120, [],'load_env',true, 'VNS_duty_cycle', [21 3]);
VNS_dat.good_channels = ~is_bad_chan;

VNS_dat = artifact_detector_3(VNS_dat);
save([root_dir '/VNS_dat_env.mat'],'VNS_dat' ,'-v7.3')

out_table_mean = generate_block_analysis(VNS_dat,20);

out_table_lump = generate_block_analysis_lump(VNS_dat,20,brain_dir,'VNS01','both');

[kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10);




%% Phase Stuff:
clear VNS_dat
VNS_dat_phase = VNS_process_raw(dat_dir, 127 ,1:120, [],'load_env',false,'load_phase',true,'VNS_duty_cycle', [21 3]);
VNS_dat_phase.good_channels = ~is_bad_chan;
save([root_dir '/VNS_dat_phase.mat'],'VNS_dat_phase' ,'-v7.3')

out_table_phase = generate_block_analysis_phase(VNS_dat_phase,20);


%% Interference:
VNS_dat_artifact_env = VNS_process_raw(dat_dir, 127, 1:127,[], 'load_env', false, 'notch_data', false,'load_vns_band',true, 'VNS_duty_cycle', [21 3]); % must mod VNS env stuff, to get things to 10hz
VNS_dat_artifact_env.good_channels = ~is_bad_chan; 
save([root_dir '/VNS_dat_artifact_env.mat'],'VNS_dat_artifact_env' ,'-v7.3')
% ekg_channel is ch 7, 1:6 are eeg
x = corrcoef(VNS_dat_artifact_env.ecog_vns_env');
vns_corr_env = abs(x(7,8:127)); % ekg channel is actually first due to sort
x = corrcoef(VNS_dat_artifact_env.ecog_vns_raw(:,VNS_dat_artifact_env.is_stim_on)');
vns_corr_raw = abs(x(7,8:127));

figure; 
plot(1:length(vns_corr_env), vns_corr_env, 'Marker','x', 'LineStyle','none', 'Color','r');
xlabel('Channel Number'); ylabel('Correlation to VNS Intensity')
title({'Effect of Stimulation Settings on the Correlation Between the Intensity of'; 'VNS Harmonics on Electrodes and EKG Channel'})
hold on;
plot(1:length(vns_corr_raw), vns_corr_raw, 'Marker','x', 'LineStyle','none', 'Color','b');
legend('VNS Frequency Envelope', 'Raw VNS Frequency')


vns_corr_env(~VNS_dat_artifact_env.good_channels) = 0; % assumes that eeg & ekg were dropped
plot_brain_elecs_dat(vns_corr_env, brain_dir,'VNS01', 'both')
title({'Correlation Between the ECoG and EKG Power in VNS Harmonics';''})
view([0 30]);
