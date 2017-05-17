function [out_table] = generate_block_analysis_phase(VNS_dat, varargin)
%% This function generates the standard set of comparisons of the circular 
% means variances and pairwise correlations on the PHASE data
% between the time blocks when stimulation is present
% and the time blocks when stimulation is absent.
%
% Inputs:
% VNS_dat - standard VNS data structure with phase data included
%
% variable inputs:
%
% 1 - window_size (the size the analysis window for blocks relative to
% onset (default 30 - means compare means of data 30 after onsets to means
% 30 seconds before onsets) - is a variable with VNS duty cycle
%
% 2 - plot_mean_dist (boolean, default true), plot the histogram of the
% circular means for each band
%
% 3 - p_thresh, the p value threshold for significance
%
% 4 - plot_sig_only - (default false) contingent on plot_mean_dist (only plot electrodes that show significance) 
%
% Here significance is determined by the Kuiper Test, except the cicular
% correlations which use the KS-test.
%
% Output: out_table, table containing the percentage of of electrodes that
% show significant difference between stim on and stim off blocks.
% also makes some graphs...

%% Variable Inputs:
window_size = 30; % time length of each block
if length(varargin)>0
    if ~isempty(varargin{1})
        window_size = varargin{1};
    end
end
plot_mean_dist = true;
if length(varargin)>1
    if ~isempty(varargin{2})
        plot_mean_dist = varargin{2};
    end
end
p_thresh = 0.05; %/sum(is_good_chan);
if length(varargin)>2
    if ~isempty(varargin{3})
        p_thresh = varargin{3};
    end
end
plot_only_sig = false;
if length(varargin)>3
    if ~isempty(varargin{4})
        plot_only_sig = varargin{4};
    end
end



% Analysis
table_cols = {'mean', 'var', 'pw_corr'};
classes = {'theta', 'alpha', 'beta', 'lg', 'hg'};
data_table = zeros(length(classes), length(table_cols));


is_good_chan = VNS_dat.good_channels;
onsets = VNS_dat.stim_onsets_inds;
fs = VNS_dat.sampFreq;
%p_thresh = 0.05; %/sum(is_good_chan);

nan_bad_times = false;
if sum(strncmpi(fieldnames(VNS_dat),'is_bad_time',11)) > 0;
    nan_bad_times = true;
end






pvals_mean = zeros(sum(is_good_chan),length(classes));
for b = 1:length(classes)
    ecog = getfield(VNS_dat, ['ecog_' classes{b} '_ang']);
    ecog = ecog(is_good_chan,:);
    if nan_bad_times
        ecog(:,is_bad_timept) = NaN;
    end

    [erps,time_axis] = make_vns_erps(ecog, onsets,fs, [-window_size window_size]);
    %erps(:,:, wayy_too_high_trial) = [];
    is_post_onset = (time_axis > 0);
    is_pre_onset = (time_axis <0);
    % means
    means_post = squeeze(circ_mean(erps(:,is_post_onset,:),[],2));
    means_pre = squeeze(circ_mean(erps(:,is_pre_onset,:),[],2));
    pvals = zeros(size(means_post,1),1);
    for k = 1:size(means_post,1)
        pvals(k) = circ_kuiper_eq(means_post(k,:), means_pre(k,:));
    end
    pvals_mean(:,b) = pvals;
    %[~,pvals] = ttest2(means_post,means_pre,'Dim',2);
    data_table(b,1) = 100*sum(pvals<p_thresh)/sum(is_good_chan);
    % vars:
    vars_post = squeeze(circ_var(erps(:,is_post_onset,:),[],[],2));
    vars_pre = squeeze(circ_var(erps(:,is_pre_onset,:),[],[],2));
    %[~,pvals] = ttest2(vars_post,vars_pre,'Dim',2);
    pvals = zeros(size(vars_post,1),1);
    for k = 1:size(vars_post,1)
        pvals(k) = circ_kuiper_eq(vars_post(k,:), vars_pre(k,:));
    end
    data_table(b,2) = 100*sum(pvals<p_thresh)/sum(is_good_chan);
    
    % pair wise corr:
    pw_corr_trial_pre = zeros(sum(is_good_chan),sum(is_good_chan),size(erps,3));
    pw_corr_trial_post = zeros(sum(is_good_chan),sum(is_good_chan),size(erps,3));
    for k = 1:size(erps,3)
        pw_corr_trial_pre(:,:,k) = circ_corcoeff(squeeze(erps(:,is_pre_onset,k))');
        pw_corr_trial_post(:,:,k) = circ_corcoeff(squeeze(erps(:,is_post_onset,k))');
    end
    
    for k = 1:size(pw_corr_trial_pre,1);
         pw_corr_trial_pre(k,k,:) = 1;
    end
    %
    %[~,pvals] = ttest2(pw_corr_trial_pre,pw_corr_trial_post, 'Dim',3);
    for i = 1:size(pw_corr_trial_pre,1)
        for j = 1:size(pw_corr_trial_pre,2)
            [~,pvals(i,j)] = kstest2(squeeze(pw_corr_trial_pre(i,j,:)),squeeze(pw_corr_trial_post(i,j,:)));
        end
    end
    data_table(b,3) = 100*sum(pvals(:)<p_thresh)/(size(pvals,1)^2-size(pvals,1));
end

out_table = table(classes','VariableNames',{'Band'});
out_table = [out_table array2table(data_table, 'VariableNames', table_cols)];


if plot_mean_dist
    for b = 1:length(classes)
        %is_good_chan
        ecog = getfield(VNS_dat, ['ecog_' classes{b} '_ang']);
        ecog = ecog(is_good_chan,:);
        
        anatomy_elecs = VNS_dat.anatomy_elecs(is_good_chan);
        if plot_only_sig
            is_sig = (pvals_mean(:,b) < p_thresh);
            ecog = ecog(is_sig,:);
            anatomy_elecs = anatomy_elecs(is_sig);
        end

        if nan_bad_times
            ecog(:,is_bad_timept) = NaN;
        end

        [erps,time_axis] = make_vns_erps(ecog, onsets,fs, [-window_size window_size]);
        is_post_onset = (time_axis > 0);
        is_pre_onset = (time_axis <0);

        means_post = squeeze(circ_mean(erps(:,is_post_onset,:),[],2));
        means_pre = squeeze(circ_mean(erps(:,is_pre_onset,:),[],2));
        pvals = zeros(size(means_post,1),1);
        for k = 1:size(means_post,1)
            pvals(k) = circ_kuiper_eq(means_post(k,:), means_pre(k,:));
        end
    
        
        sub_plot_coords = position_subplot_grid(size(erps,1),0.2);
        figure;
        for i = 1:size(erps,1)        
            subplot('position', sub_plot_coords(i,:)); hold on;
            max_val = max([means_pre(i,:) means_post(i,:)]);
            min_val = min([means_pre(i,:) means_post(i,:)]);
            nbins = ceil(sqrt(size(means_pre,2)));
            histogram(means_pre(i,:), linspace(min_val, max_val, nbins),'EdgeColor','none')
            histogram(means_post(i,:), linspace(min_val, max_val, nbins), 'EdgeColor','none')
            axis tight;
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            title(anatomy_elecs{i}, 'FontSize',8)
            if ~plot_only_sig & (pvals(i) < p_thresh)
                set(gca, 'XColor','r'); set(gca, 'YColor', 'r');
                set(gca, 'LineWidth', 3);
            end
        end
    end
end
a = 1;

end

    
