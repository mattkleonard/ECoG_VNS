function [kl_spectrograms, time_axis, f_axis] = kl_get_spectrograms(VNS_dat, plot_out, varargin);
% function makes a quick approximate of a spectrogram by taking Analytic
% Envelope for several data
% 
% Inputs:
% VNS_dat - VNS output structure must contain envelope data
% plot_out - generate a channel wise spectrogram
% 
% Variable inputs:
% 1 - anatomical label (n_chans cell array)
% 2 - electrode colors - useful for comparing to brain anatomy\
% 3 - window length (default 10)
% 4 - kl_window (default [-30, 60] (seconds about onsets to study)
% 5 - image save_directory
% 6 - figure title appelation (only possible if Exist variable input 4
% 7 - is_sig (channels x frequency band matrix of boolean for significant
% electrodes will only plot electrodes in the spectrogram which are significant
% in some band).
% outputs:
% kl_spectrograms (n_chan x time x frequency matrix)
% time_axis = time axis for kl_spectrograms
% f_axis = frequency axis for kl_spectrograms
%
is_good = VNS_dat.good_channels;
bands = {'theta', 'alpha', 'beta', 'lg', 'hg'};
color_chans = false;
label_anatomy = false;
window_len = 10;
plot_spectrograms = true;
plot_time_series = true; % additional plot flag to plot the time series for each band
if length(varargin) > 0
    if ~isempty(varargin{1})
        label_anatomy = true;
        anatomy = varargin{1};
        anatomy = anatomy(is_good);
    end
 end
if length(varargin) > 1
    if ~isempty(varargin{2})
        color_chans = true;
        elec_colors = varargin{2};
    end
end
if length(varargin) > 2
    if ~isempty(varargin{3})
        window_len = varargin{3};
    end
end
kl_win = [-30 60];
if length(varargin) > 3
    if ~isempty(varargin{4})
        kl_win = varargin{4};
    end
end

save_figs = false;
if length(varargin) > 4
    if ~isempty(varargin{5})
        save_figs = true;
        fig_dir = varargin{5};
        figure_appellation = '';
    end
end
if save_figs & (length(varargin) > 5)
    if ~isempty(varargin{6})
        figure_appellation = varargin{6};
    end
end
is_sig = true(sum(is_good), length(bands));
if length(varargin) > 6
    if ~isempty(varargin{7})
        is_sig = varargin{7};
    end
end



% NaN bad times:
nan_bad_times = false;
if sum(strncmpi(fieldnames(VNS_dat),'is_bad_time',11)) > 0;
    nan_bad_times = true;
end
% % % 
% % % if sum(strcmpi(fieldnames(VNS_dat),'is_bad_timept')) == 1;
% % %     is_good_time = ~VNS_dat.is_bad_timept;
% % % else
% % %     is_good_time = true(size(VNS_dat.is_stim_on));
% % % end

%[kl_div, time_axis] = tstat_time_course(VNS_dat.ecog_theta_env(is_good,:), VNS_dat.stim_onsets_inds, VNS_dat.sampFreq, window_len);
ecog_dat = VNS_dat.ecog_theta_env(is_good,:);
%ecog_dat(:,~is_good_time) = NaN;
[kl_div, time_axis] = kl_divergence_time_course_2(ecog_dat, VNS_dat.stim_onsets_inds, VNS_dat.sampFreq, window_len);
f_axis = mean([4,7; 8,15; 18, 30; 33 55; 70 150],2);  % mean frequencies of each band


kl_spectrograms = zeros(size(kl_div,1), size(kl_div,2), length(f_axis));


%[kl_div, time_axis] = tstat_time_course(VNS_dat.ecog_theta_env(is_good,:), VNS_dat.stim_onsets_inds, VNS_dat.sampFreq, window_len);
for b = 1:length(bands)
    ecog_dat = getfield(VNS_dat, ['ecog_' bands{b} '_env']);
    ecog_dat(~is_good,:) = [];
    if nan_bad_times
        is_bad_timept = getfield(VNS_dat, ['is_bad_time_' bands{b}]);
        ecog(:,is_bad_timept) = NaN;
    end
    [kl_div, time_axis] = kl_divergence_time_course_2(ecog_dat, VNS_dat.stim_onsets_inds, VNS_dat.sampFreq, window_len, kl_win);
    kl_spectrograms(:,:,b) = kl_div;
%[kl_div, time_axis] = tstat_time_course(VNS_dat.ecog_alpha_env(is_good,:), VNS_dat.stim_onsets_inds, VNS_dat.sampFreq, window_len);
end


cmax = max(kl_spectrograms(:));
cmin = min(kl_spectrograms(:));
if plot_out
    is_sig_ch = sum(is_sig,2)>3;
    kl_spectrograms_sig = kl_spectrograms(is_sig_ch,:,:);
    anatomy_sig = anatomy(is_sig_ch);
    if plot_spectrograms
        [sub_plot_coords] = position_subplot_grid(size(kl_spectrograms_sig,1),0.25);
        figure('units','normalized','outerposition',[0 0 1 1]);
        for n = 1:size(kl_spectrograms_sig,1)
            dat = squeeze(kl_spectrograms_sig(n,:,:));
            subplot('position', sub_plot_coords(n,:)); hold on;
            imagesc(time_axis, f_axis, dat')
            caxis([cmin cmax])
            set(gca,'YDir', 'normal')
            set(gca, 'XTick', [])
            %set(gca, 'YTick', [])
            axis tight;
            hold on;
            plot([0 0], get(gca,'YLim'), 'r')
            plot([34 34], get(gca,'YLim'), 'r')
            ax = gca;
            if color_chans
                ax.LineWidth = 4;
                ax.XColor = elec_colors(n,:);
                ax.YColor = elec_colors(n,:);
            end
            if label_anatomy
                title(anatomy_sig{n}, 'FontSize',8)
            end
        end
        if save_figs
            drawnow
            pause(5)
            set(gcf, 'PaperPositionMode', 'auto')
            saveas(gcf,[fig_dir filesep 'KL_Divergence_Spectrograms_' figure_appellation '.png']);
        end
    end
    if plot_time_series
        band_names = {'Theta', 'Alpha', 'Beta', 'Low Gamma', 'High Gamma'};
        figure('units','normalized','outerposition',[0 0 0.7 0.7]);
        for k = 1:length(band_names)
            subplot(2,3,k)
            kl_band = squeeze(kl_spectrograms(is_sig(:,k),:,k));
            plot(repmat(time_axis,size(kl_band,1),1)', kl_band')
            hold on; plot([0 0], get(gca, 'YLim'),'k'); plot([0 0], get(gca, 'YLim'),'k')
            xlabel('Time (s)'); ylabel('KL-Divergence'); legend('Ch 1','Ch 2','...')
            title(['Time Series of KL-Divergence of ' band_names{k} ' Envelope'])
            %title(['Time Series of T-statistic of ' band_names{k} ' Envelope'])
            %title(['Time Series of KS-Statistic of ' band_names{k} ' Envelope'])

            axis tight
            xlim([-15, 45])
        end
        if save_figs
            drawnow
            pause(5)
            set(gcf, 'PaperPositionMode', 'auto')
            saveas(gcf,[fig_dir filesep 'KL_Divergence_TimeCourse_' figure_appellation '.png']);
        end
    end
a=1;


end





