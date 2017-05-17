function visually_inspect_channels(VNS_dat)

if sum(strcmp(fieldnames(VNS_dat), 'ecog_theta_env')) == 1

    bands = {'ecog_theta_env','ecog_alpha_env','ecog_beta_env','ecog_lg_env','ecog_hg_env'};
    blank_edge_art_flag = 1;
    plot_evnt_flag = 0;

    nSec = 2;

    %%

    figure;
    for j = 1:length(bands)
        subplot(length(bands),1,j);
        %     for i = 1:length(VNS_dat.stim_onsets_inds)
        %         fprintf('%d\n',i);
        %         plot(VNS_dat.(bands{j})(:,(VNS_dat.stim_onsets_inds(i) - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds(i) + VNS_dat.sampFreq*nSec))); %VNS_dat.stim_offset_inds(i)));
        plot(VNS_dat.(bands{j})');
        hold on;
        %     end
        set(gca,'XLim',[-1000 26000]);
    end

    %% FIND BAD CHANNELS

    % figure;
    for j = 1:length(bands)
        fprintf('Plotting band %s\n',bands{j});
        figure('Name',bands{j});
        for i = 1:size(VNS_dat.(bands{j}),1)
            p = plotGridPosition(i,size(VNS_dat.(bands{j}),1),ceil(sqrt(size(VNS_dat.(bands{j}),1))));
            subplot('Position',p);

            plotDat = VNS_dat.(bands{j})(i,:);
            if blank_edge_art_flag
                plotDat(1:20) = NaN;
                plotDat(end-20:end) = NaN;
            end
            plot(plotDat,'Color',[0.8 0.1 0.1]);
            axis tight;
            set(gca,'YLim',[min(min(plotDat)) max(max(plotDat))],...
                'XTickLabel',[],'YTickLabel',[]);
            text(0,(max(get(gca,'YLim')) - ((max(get(gca,'YLim')))*0.1)),num2str(i));
            line(get(gca,'XLim'),[0 0],'Color','k');
            %     shadedErrorBar(1:length((VNS_dat.stim_onsets_inds(i) - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds(i) + VNS_dat.sampFreq*nSec)),...
            %         mean(VNS_dat.(bands{j})(:,VNS_dat.stim_onsets_inds - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds + VNS_dat.sampFreq*nSec)),...
            %         std(VNS_dat.(bands{j})(:,VNS_dat.stim_onsets_inds - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds + VNS_dat.sampFreq*nSec)));
            hold on;

            if plot_evnt_flag
                for k = 1:length(VNS_dat.stim_onsets_inds)
                    line([VNS_dat.stim_onsets_inds(k) VNS_dat.stim_onsets_inds(k)],get(gca,'YLim'),'Color','b')
                end
            end
        end
    end
else
    
    blank_edge_art_flag = 1;
    plot_evnt_flag = 1;

    nSec = 2;

    %%

    figure;
        %     for i = 1:length(VNS_dat.stim_onsets_inds)
        %         fprintf('%d\n',i);
        %         plot(VNS_dat.(bands{j})(:,(VNS_dat.stim_onsets_inds(i) - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds(i) + VNS_dat.sampFreq*nSec))); %VNS_dat.stim_offset_inds(i)));
    plot(VNS_dat.raw_erps(:));

    %% FIND BAD CHANNELS

    figure;
    for i = 1:size(VNS_dat.raw_erps,1)
        p = plotGridPosition(i,size(VNS_dat.raw_erps,1),ceil(sqrt(size(VNS_dat.raw_erps,1))));
        subplot('Position',p);

        plotDat = VNS_dat.raw_erps(i,:);
        if blank_edge_art_flag
            plotDat(1:20) = NaN;
            plotDat(end-20:end) = NaN;
        end
        plot(plotDat,'Color',[0.8 0.1 0.1]);
        axis tight;
        set(gca,'YLim',[min(min(plotDat)) max(max(plotDat))],...
            'XTickLabel',[],'YTickLabel',[]);
        text(0,(max(get(gca,'YLim')) - ((max(get(gca,'YLim')))*0.1)),num2str(i));
        line(get(gca,'XLim'),[0 0],'Color','k');
            %     shadedErrorBar(1:length((VNS_dat.stim_onsets_inds(i) - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds(i) + VNS_dat.sampFreq*nSec)),...
            %         mean(VNS_dat.(bands{j})(:,VNS_dat.stim_onsets_inds - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds + VNS_dat.sampFreq*nSec)),...
            %         std(VNS_dat.(bands{j})(:,VNS_dat.stim_onsets_inds - VNS_dat.sampFreq*nSec):(VNS_dat.stim_onsets_inds + VNS_dat.sampFreq*nSec)));
        hold on;

        if plot_evnt_flag
            for k = 1:length(VNS_dat.stim_onsets_inds)
                %line([VNS_dat.stim_onsets_inds(k) VNS_dat.stim_onsets_inds(k)],get(gca,'YLim'),'Color','b')
            end
        end
    end
    
end

end
