function [out_table, pvals] = generate_block_analysis_lump(VNS_dat, varargin)
% This function generates the analysis comparing ALL data points in the time window
% before stimulation to All data points in the time window after stimulation
% It retunrs a table showing the percentage of electrodes in each band that
% show a significant difference betwen these two populations and also list
% the average effect sizes among the significant electrodes
% as a variable option in will plot the magnitude of the signifcance
% (Cohen's D) on the brain.
%
% Inputs:
% VNS_dat - VNS data structure that includes envelope data
%
% Variable Inputs:
% 1 - window length(pre vs post stim)
% 2 - the directory of brain data (necessary for brain plots)
% 3 - subject name
% 4 - hemisphere to plot (default 'lh')
% 5 - p_thresh, the p value threshold for significance (default
% 0.05/nchans)
%
% Outputs: Data Table containing percentage of significant electrodes per
% band and the average effect size of among these electrodes.
% Also can plot brain plots of this effect.


window_size = 30; % time length of each block
if length(varargin)>0
    if ~isempty(varargin{1})
        window_size = varargin{1};
    end
end

%brain_dir = '/Users/changlab/Documents/data/VNS/EC131/BrainPlot';
if length(varargin)>1
    if ~isempty(varargin{2})
        brain_dir = varargin{2};
    end
end

include_subj_name = false;
if length(varargin)>2
    if ~isempty(varargin{3})
        include_subj_name = true;
        subj = varargin{3};
    end
end
include_hemisphere = false;
if length(varargin)>3
    if ~isempty(varargin{4})
        include_hemisphere = true;
        hem = varargin{4};
    end
end


% Analysis
table_cols = {'mean', 'mean_sig_effect_size'};
classes = {'theta', 'alpha', 'beta', 'lg', 'hg'};
data_table = zeros(length(classes), length(table_cols));


is_good_chan = VNS_dat.good_channels;
onsets = VNS_dat.stim_onsets_inds;
fs = VNS_dat.sampFreq;

p_thresh = 0.05/sum(is_good_chan);
if length(varargin)>4
    if ~isempty(varargin{5})
        p_thresh = varargin{5};
    end
end

nan_bad_times = false;
if sum(strncmpi(fieldnames(VNS_dat),'is_bad_time',11)) > 0;
    nan_bad_times = true;
end


effect_sizes = zeros(sum(is_good_chan),length(classes));
pvals = zeros(sum(is_good_chan),length(classes));

%[erps,time_axis] = make_vns_erps(VNS_dat.ecog_hg_env, onsets,fs, [-window_size window_size]);
for b = 1:length(classes)
    classes{b};
    ecog = getfield(VNS_dat, ['ecog_' classes{b} '_env']);
    ecog(~is_good_chan,:) = [];
    if nan_bad_times
        is_bad_timept = getfield(VNS_dat, ['is_bad_time_' classes{b}]);
        ecog(:,is_bad_timept) = NaN;
    end

    [erps,time_axis] = make_vns_erps(ecog, onsets,fs, [-window_size window_size]);

    is_post_onset = (time_axis > 0);
    is_pre_onset = (time_axis <0);
    % means
    erps_post = erps(:,is_post_onset,:); dims = size(erps_post);
    ecog_post = reshape(erps_post, dims(1), prod(dims(2:3)));
    
    erps_pre = erps(:,is_pre_onset,:); dims = size(erps_pre);
    ecog_pre = reshape(erps_pre, dims(1), prod(dims(2:3)));
    
    [~,pvals(:,b)] = ttest2(ecog_post,ecog_pre,'Dim',2);
    data_table(b,1) = 100*sum(pvals(:,b)<p_thresh)/sum(is_good_chan);
    effect_sizes(:,b) = (nanmean(ecog_post,2) - nanmean(ecog_pre,2))./nanstd([ecog_post ecog_pre],[],2);
    data_table(b,2) = mean(effect_sizes(pvals(:,b)<p_thresh,b));
end

out_table = table(classes','VariableNames',{'Band'});
out_table = [out_table array2table(data_table, 'VariableNames', table_cols)];

plot_erps = false;
if plot_erps
    twin = [-10 10];
    ecog = getfield(VNS_dat, ['ecog_hg_env']);
    ecog(~is_good_chan,:,:) = [];
    if nan_bad_times
        ecog(:,is_bad_timept) = NaN;
    end
    [erps,time_axis] = make_vns_erps(ecog, onsets,fs, twin);

    anatomy_elecs = VNS_dat.anatomy_elecs(is_good_chan);

    is_twin = (time_axis >= twin(1)) & (time_axis<=twin(2));
    time_axis(~is_twin) = [];
    erps(~is_good_chan,:,:) = [];
    erps(:,~is_twin,:) = [];
    sub_plot_coords = position_subplot_grid(size(erps,1),0.2);
    figure;
    for i = 1:size(erps,1);
        subplot('position', sub_plot_coords(i,:)); hold on;
        shadedErrorBar(time_axis, mean(squeeze(erps(i,:,:)),2), nansem(squeeze(erps(i,:,:)),2),'r',1);
        axis tight;
        ax = gca;
        plot(ax.XLim, [0 0],'k'); plot([0 0], ax.YLim, 'k'); plot([2 2], ax.YLim,'m');
        ax.XTick = [];
        ax.YTick = [];
        axis tight;
        title(anatomy_elecs{i}, 'FontSize',8)
    end
end


plot_brain = false;
if plot_brain
    %brain_dir = '/Users/changlab/Documents/data/VNS/EC131/BrainPlot';
    for b = 1:length(classes)
        dat = effect_sizes(:,b);
        dat(pvals(:,b)>=p_thresh) = 0;
        dat_full = zeros(length(is_good_chan),1);
        dat_full(is_good_chan) = dat;
        if include_subj_name
            if include_hemisphere
                plot_brain_elecs_dat(dat_full, brain_dir,subj,hem);
            else
                plot_brain_elecs_dat(dat_full, brain_dir,subj);
            end
        else
            plot_brain_elecs_dat(dat_full, brain_dir);
        end
        view(240, 0);
        title(['Effect Size of Change in ' classes{b} ' during VNS']);
    end
    a = 1;
end
a = 1;

end

    
