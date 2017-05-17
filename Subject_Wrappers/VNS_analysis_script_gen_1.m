%% General VNS Processing Script.
% as per the usual, you will need to make ammendments as necessary
% for an individual subject
% in principle this script could be automated into a function,
% however given the long run times for generating data structures, and the
% importance of Visually inspecting Figures as they are generated,
% it seems premature to carry out such a design

%% Subject parameters:
root_dir = '/Users/changlab/Documents/data/VNS/EC131';
dat_dir = [root_dir '/EC131_759989b0']; % Contains a Block of .h5 files
brain_dir = [root_dir '/BrainPlot']; % Contains clinical electrode map and mesh

fs_data = 1024;
duty_cycle = [30 2];
analysis_win = duty_cycle(1);
kl_win = [-duty_cycle(1) 2*duty_cycle(1)];

subj = 'EC131';
hem = 'lh';
p_thresh = 0.05;

ekg_chan = 137;
ecog_elecs = 1:136;



%% Find bad channels:
VNS_raw = VNS_process_raw(dat_dir, ekg_chan ,ecog_elecs, [],'load_RAW_erps',true,'load_env',false,'load_full_dir',false,'notch_data',false,'VNS_duty_cycle', duty_cycle);
[is_bad_chan] = plot_spectra_debug(VNS_raw, fs_data, 0.6);
clear VNS_raw
%% Find bad times:
VNS_dat = VNS_process_raw(dat_dir, ekg_chan, ecog_elecs, [],'load_env',true,'VNS_duty_cycle', duty_cycle);
VNS_dat.good_channels = ~is_bad_chan;

VNS_dat = artifact_detector_3(VNS_dat);
save([root_dir '/VNS_dat_env.mat'],'VNS_dat' ,'-v7.3')

%% Get summary stats on full on/off
out_table_mean = generate_block_analysis(VNS_dat,analysis_win);
out_table_mean = generate_block_analysis(VNS_dat,analysis_win,[],[],true); % plot only significant erps...


[out_table_lump, pvals_onoff] = generate_block_analysis_lump(VNS_dat,analysis_win,brain_dir,subj,hem);

is_sig_diff = pvals_onoff<(p_thresh/size(pvals_onoff,1));
[kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10,kl_win); % plot all channels
[kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, true, VNS_dat.anatomy_elecs,[],10,kl_win, [],[],is_sig_diff);

clear VNS_dat
%% Phase Stuff:
VNS_dat_phase = VNS_process_raw(dat_dir, ekg_chan, ecog_elecs, [],'load_env',false,'load_phase',true,'VNS_duty_cycle', duty_cycle);
VNS_dat.good_channels = ~is_bad_chan;
save([root_dir '/VNS_dat_phase.mat'],'VNS_dat_phase' ,'-v7.3')

%out_table_phase = generate_block_analysis_phase(VNS_dat_phase,analysis_win);
out_table_phase = generate_block_analysis_phase(VNS_dat_phase,analysis_win, true, p_thresh,true);

out_table_phase_lump = generate_block_analysis_lump_phase(VNS_dat_phase, analysis_win,brain_dir,'EC131','lh');



%% Peak Frequency Bands:
VNS_dat_peak_band = VNS_process_raw(dat_dir, ekg_chan, ecog_elecs, [],'load_env',false, 'load_peak_freqs', true,'VNS_duty_cycle', duty_cycle);
out_table_peak = zeros(5,2);
for k = 1:5
    dat_pre = squeeze(VNS_dat_peak_band.peak_freq_pre_all(VNS_dat_peak_band.good_channels,k,:));
    dat_post = squeeze(VNS_dat_peak_band.peak_freq_post_all(VNS_dat_peak_band.good_channels,k,:));
%    [~,pvals] = ttest2(dat_pre,dat_post,'Dim',2);
    for j = 1:size(dat_pre,1)
        [~,pvals(j)] = kstest2(dat_pre(j,:),dat_post(j,:));
    end
    is_sig = pvals <p_thresh;
    out_table_peak(k,1) = sum(is_sig)/length(pvals);
    out_table_peak(k,2) = mean((mean(dat_post(is_sig,:),2) - mean(dat_pre(is_sig,:),2))./std([dat_post(is_sig,:) dat_pre(is_sig,:)],[],2));
end




%% Interference:
VNS_dat_artifact_env = VNS_process_raw(dat_dir, ekg_chan, [ecog_elecs ekg_chan],[], 'load_env', false, 'notch_data', false,'load_vns_band',true, 'VNS_duty_cycle', duty_cycle); % must mod VNS env stuff, to get things to 10hz
VNS_dat_artifact_env.good_channels = ~is_bad_chan; 
save([root_dir '/VNS_dat_artifact_env.mat'],'VNS_dat_artifact_env' ,'-v7.3')

]x,p] = corrcoef(VNS_dat_artifact_env.ecog_vns_env');
vns_corr_env = abs(x(1,2:size(x,2))); % ekg channel is actually first due to sort
vns_corr_env_p = p(1,2:size(x,2)) < p_thresh;
[x,p] = corrcoef(VNS_dat_artifact_env.ecog_vns_raw(:,VNS_dat_artifact_env.is_stim_on)');
vns_corr_raw = abs(x(1,2:137));
vns_corr_raw_p = p(1,2:size(x,2)) < p_thresh;

figure; 
plot(1:length(vns_corr_env), vns_corr_env, 'Marker','x', 'LineStyle','none', 'Color','r');
xlabel('Channel Number'); ylabel('Correlation to VNS Intensity')
title({'Effect of Stimulation Settings on the Correlation Between the Intensity of'; 'VNS Harmonics on Electrodes and EKG Channel'})
hold on;
plot(1:length(vns_corr_raw), vns_corr_raw, 'Marker','x', 'LineStyle','none', 'Color','b');
legend('VNS Frequency Envelope', 'Raw VNS Frequency')


figure;
plot(vns_corr_raw, vns_corr_env, 'rx')
xlabel('Correlation Raw VNS Signal'); ylabel('Correlation to VNS Envelope');
title('Comparison of Measures of VNS Interference')
annotation('textbox',[0.7 0.8 0.1 0.04],'String',['r = ' num2str(corr(vns_corr_raw', vns_corr_env'),3)])



vns_corr_env(~vns_corr_env_p) = 0; % assumes that eeg & ekg were dropped
plot_brain_elecs_dat(vns_corr_env, brain_dir,subj, hem)
title({'Correlation Between the ECoG and EKG Power in VNS Harmonics'})
view(230,0)

%vns_corr_raw(~vns_corr_raw_p) = 0;
%plot_brain_elecs_dat(vns_corr_env, brain_dir,w, 'lh')
%view(230,0)
%title({'Correlation Between the ECoG and EKG Signal in VNS Harmonics'})
