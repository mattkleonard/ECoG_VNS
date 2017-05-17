figure; subplot(2,1,1);
time_axis = linspace(0,length(data)/(60*fs_dat),length(data));
plot(time_axis,data)
xlabel('Time (min)'); ylabel('Intensity'); title('RAW EKG Signal')
ylim([-300 600]); axis tight
legend('Raw Data')

subplot(2,1,2);
plot(time_axis,dat_filt)
hold on; plot(time_axis, vns_stim_env, 'r','LineWidth',2)
xlabel('Time (min)'); ylabel('Intensity'); title('Power EKG Channel Filtered for VNS Harmonics')
legend({'Filtered EKG Signal','Power in VNS Harmonics'})

% get time locking
figure;
subplot(1,2,1)
plot(time_axis,dat_filt);
hold on; plot(time_axis, vns_stim_env, 'r')
thresh = mean([median(vns_stim_env(vns_stim_env > mean(vns_stim_env))) median(vns_stim_env(vns_stim_env < mean(vns_stim_env)))]);
plot([0 60], [thresh thresh], 'k--') 
plot(time_axis, 180*(vns_stim_env > thresh), 'k', 'LineWidth',2)
title('Estimation of Stimulation Onset Times')
legend({'Filtered EKG Data', 'VNS Power', 'Threshold', 'Estimated ON Times'})

subplot(1,2,2);
vns_pulse_ideal = 0.5*(generate_ideal_vns_pulse(100)+1);
plot(linspace(0,44,44*100), vns_pulse_ideal,'k', 'LineWidth',2);
hold on;
plot(time_axis, vns_stim_env/max(vns_stim_env),'r', 'LineWidth',2);
xlim([0, 50]); ylim([0 1.2])
xlabel('Time (s)'); ylabel('Intensity'); 
title({'Compare Time Courses of Theoretical and'; 'Observed VNS Duty Cycle in Selected ROI'})
legend('Idealized VNS Cycle', 'Measured VNS Power')


%final out
figure; 
plot(linspace(0,length(ekg_data)/(60*fs_in),length(ekg_data)), ekg_data/max(ekg_data));
hold on;
(ecg_data_vns + min(ecg_data_vns))/max(ecg_data_vns + min(ecg_data_vns));
plot(linspace(0,length(ecg_data_vns)/(60*fs_in),length(ecg_data_vns)), (ecg_data_vns - min(ecg_data_vns))/max(ecg_data_vns - min(ecg_data_vns)))
plot(linspace(0,length(is_stim_on_xcorr)/(60*fs_out), length(is_stim_on_xcorr)), is_stim_on_xcorr, 'k', 'LineWidth',2)
axis tight; ylim([-0.5 1.2])
legend('RAW Data', 'VNS Power', 'VNS ON')
xlabel('Time (min)'); ylabel('Normalized Intensity')
title('Visualization of Algorithm for finding VNS ON Times')