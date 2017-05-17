function [out_table] = generate_block_analysis(VNS_dat, varargin)
%% This function generates the standard set of comparisons of the means variances and
% pairwise correlations between the time blocks when stimulation is present
% and the time blocks when stimulation is absent.
% for the high gamma response.
%
% Inputs:
% VNS_dat - standard VNS data structure with envelope data included
%
% variable inputs:
%
% 1 - window_size (the size the analysis window for blocks relative to
% onset (default 30 - means compare means of data 30 after onsets to means
% 30 seconds before onsets) - is a variable with VNS duty cycle
%
% 2 - plot_erps (boolean, default true), plot the ERPs for each band time
% locked to stimulation onset
%
% 3 - p_thresh, the p value threshold for significance
%
% 4 - plot_only_sig (boolean default = false) only active if plot erps is
% active (will only plot significant ERPs)
% 
% Here significance is determined by the T-test Test, except the pairwise
% correlations which use the KS-test.
%
% Output: out_table, table containing the percentage of of electrodes that
% show significant difference between stim on and stim off blocks.
% also makes some graphs...

window_size = 30; % time length of each block
if length(varargin)>0
    if ~isempty(varargin{1})
        window_size = varargin{1};
    end
end
plot_erps = true;
if length(varargin)>1
    if ~isempty(varargin{2})
        plot_erps = varargin{2};
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

nan_bad_times = false;
if sum(strncmpi(fieldnames(VNS_dat),'is_bad_time',11)) > 0;
    nan_bad_times = true;
end





%[erps,time_axis] = make_vns_erps(VNS_dat.ecog_hg_env, onsets,fs, [-window_size window_size]);
%wayy_too_high = erps>12;
%wayy_too_high_trial = (sum(squeeze(sum(wayy_too_high(is_good_chan,:,:),1)),1)>1);
pvals_means = zeros(sum(is_good_chan),length(classes));
for b = 1:length(classes)
    classes{b};
    ecog = getfield(VNS_dat, ['ecog_' classes{b} '_env']);
    ecog = ecog(is_good_chan,:);
    if nan_bad_times
        is_bad_timept = getfield(VNS_dat, ['is_bad_time_' classes{b}]);
        ecog(:,is_bad_timept) = NaN;
    end

    [erps,time_axis] = make_vns_erps(ecog, onsets,fs, [-window_size window_size]);
    %erps(:,:, wayy_too_high_trial) = [];
    is_post_onset = (time_axis > 0);
    is_pre_onset = (time_axis <0);
    % means
    pvals = zeros(size(erps,1),1); % channels, 0
    means_post = squeeze(nanmean(erps(:,is_post_onset,:),2));
    means_pre = squeeze(nanmean(erps(:,is_pre_onset,:),2));
    [~,pvals] = ttest2(means_post,means_pre,'Dim',2);
    pvals_means(:,b) = pvals; % p value for each channel 
    data_table(b,1) = 100*sum(pvals<p_thresh)/sum(is_good_chan); % why this way?
    % vars:
    pvals = zeros(size(erps,1),1);
    vars_post = squeeze(nanvar(erps(:,is_post_onset,:),[],2));
    vars_pre = squeeze(nanvar(erps(:,is_pre_onset,:),[],2));
    [~,pvals] = ttest2(vars_post,vars_pre,'Dim',2);
    data_table(b,2) = 100*sum(pvals<p_thresh)/sum(is_good_chan);
    
    % pair wise corr:
    pw_corr_trial_pre = zeros(sum(is_good_chan),sum(is_good_chan),size(erps,3)); % 
    pw_corr_trial_post = zeros(sum(is_good_chan),sum(is_good_chan),size(erps,3));
    for k = 1:size(erps,3)
        pw_corr_trial_pre(:,:,k) = corrcoef(squeeze(erps(:,is_pre_onset,k))', 'rows', 'pairwise'); % these are just generating NANs for beta
        pw_corr_trial_post(:,:,k) = corrcoef(squeeze(erps(:,is_post_onset,k))', 'rows', 'pairwise'); % why?
    end
    for k = 1:size(pw_corr_trial_pre,1);
        pw_corr_trial_pre(k,k,:) = 1;
    end
    %[~,pvals] = ttest2(pw_corr_trial_pre,pw_corr_trial_post, 'Dim',3);

    pvals = zeros(size(pw_corr_trial_pre,1),size(pw_corr_trial_pre,1));
    for i = 1:size(pw_corr_trial_pre,1)
        for j = 1:size(pw_corr_trial_pre,2)
            [~,pvals(i,j)] = kstest2(squeeze(pw_corr_trial_pre(i,j,:)), squeeze(pw_corr_trial_post(i,j,:)));
        end
    end
    data_table(b,3) = 100*sum(pvals(:)<p_thresh)/(size(pvals,1)^2-size(pvals,1));
end

out_table = table(classes','VariableNames',{'Band'});
out_table = [out_table array2table(data_table, 'VariableNames', table_cols)];


if plot_erps
    for b = 1:length(classes)
        anatomy_elecs = VNS_dat.anatomy_elecs(is_good_chan);
        twin = [-30 30];
        ecog = getfield(VNS_dat, ['ecog_' classes{b} '_env']);
        ecog = ecog(is_good_chan,:);
        if nan_bad_times
            is_bad_timept = getfield(VNS_dat, ['is_bad_time_' classes{b}]);
            ecog(:,is_bad_timept) = NaN;
        end
        if plot_only_sig
            is_sig = pvals_means(:,b) < p_thresh;
            ecog = ecog(is_sig,:);
            anatomy_elecs = anatomy_elecs(is_sig);
        end
        [erps,time_axis] = make_vns_erps(ecog, onsets,fs, twin); 

        

        is_twin = (time_axis >= twin(1)) & (time_axis<=twin(2));
        time_axis(~is_twin) = [];
        erps(:,~is_twin,:) = [];
        sub_plot_coords = position_subplot_grid(size(erps,1),0.2);
        figure;
        for i = 1:size(erps,1);
            subplot('position', sub_plot_coords(i,:)); hold on;
            shadedErrorBar(time_axis, nanmean(squeeze(erps(i,:,:)),2), nansem(squeeze(erps(i,:,:)),2),'r',1);
            axis tight;
            ax = gca;
            plot(ax.XLim, [0 0],'k'); plot([0 0], ax.YLim, 'k'); plot([2 2], ax.YLim,'m'); %why???
            ax.XTick = [];
            ax.YTick = [];
            axis tight;
            title(anatomy_elecs{i}, 'FontSize',8)
            if ~plot_only_sig & (pvals_means(i,b) < p_thresh)
                ax.LineWidth = 3;
                ax.XColor = 'b';
                ax.YColor = 'b';
            end
        end
    end
end
a = 1;

end

    
