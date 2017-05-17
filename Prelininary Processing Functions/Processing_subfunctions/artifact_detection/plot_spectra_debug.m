function [is_bad_chan] = plot_spectra_debug(VNS_raw, varargin)
%% This function is a plotting script to assist with debugging clinical data.
% it looks at the spectrograms of the raw data and uses the shape of a
% channel's spectrogram to determine whether a channel is measuring real
% neurological signal or not. Additionally it generates a plot of the
% spectrograms of each channel with problematic ones outlined in red so
% that a human can review the selection of bad channels.
%
% Inputs: VNS_raw - VNS data structure with Raw Data ERPs included
%
% Variable Inputs:
% var1 - sample rate (default 1024, has been seen at 1000 and 512)
%
% var2 - correlation threshold (default 0.6 - raise to reject more bad
% chans)
%
% var3 - VNS frequency (default 25hz, has been seen at 10 hz before)
% this and the wall frequency are exlucded in the estimate of pinknoise.
%
% var4 - clear white matter electrodes (electrodes labeled as white matter
% will be cleared automatically as bad channels
% 
% var5 - plot_out (default true) generates a plot of the spectrogram of
% each channel 
%
% outputs:
% is_bad_chan - 1d boolean array (ON if channel does not have a pinknoise
% distribution, or if it is labeled as whitematter
%

%% Load Variable Inputs:
fs_dat = 1024; % magic constant, may need to be revised later.
if length(varargin)>0
    if ~isempty(varargin{1});
        fs_dat = varargin{1};
    end
end

pinkNoise_thresh = 0.6; % pinkness must be along these lines to be seen as good
if length(varargin)>1
    if ~isempty(varargin{2});
        pinkNoise_thresh = varargin{2};
    end
end

vns_freq = 25;
if length(varargin)>2
    if ~isempty(varargin{3});
        vns_freq = varargin{3};
    end
end

clear_white_matter = true;
if length(varargin)>3
    if ~isempty(varargin{4});
        clear_white_matter = varargin{4};
    end
end


plot_spectra = true;
if length(varargin)>4
    if ~isempty(varargin{5});
        plot_spectra = varargin{5};
    end
end

notch_interenrence = true; % ignore frequencies with known artifacts:
if notch_interenrence
    wall_freq = 60;
    
    wall_harms = wall_freq:wall_freq:(fs_dat/2);
    vns_harms = vns_freq:vns_freq:(fs_dat/2);
    interfence_harms = [wall_harms vns_harms];
    interference_spread = 1;
end

is_bad_chan = strcmpi(VNS_raw.anatomy_elecs,'Left-Cerebral-White-Matter');
flatten = @(x) x(:);


if plot_spectra
    [sub_plot_coords] = position_subplot_grid(size(VNS_raw.raw_erps,1),0.1);
    figure;
    for n = 1:size(VNS_raw.raw_erps,1)
        ch_raw = flatten(VNS_raw.raw_erps(n,:,:));
        subplot('position', sub_plot_coords(n,:));
        [pxx,f] =pwelch(ch_raw, [],[],[],fs_dat);
        if notch_interenrence
            is_good_freq = f>30 & f<(fs_dat*0.33); % set maximum value to ignore anti aliasing
            for k = 1:length(interfence_harms)
                is_int = (f >= (interfence_harms(k)-interference_spread)) & (f <= (interfence_harms(k)+interference_spread));
                is_good_freq(is_int) = false;
            end
            is_bad_chan(n) = (corr(1./smooth(pxx(is_good_freq),10),f(is_good_freq))<pinkNoise_thresh);
            plot(f(is_good_freq),1./smooth(pxx(is_good_freq),10),'r.')
        else
            is_bad_chan(n) = (corr(1./smooth(pxx,10),f)<pinkNoise_thresh);
        end
        semilogy(f,pxx)
        axis tight
        %ylim([10^(-6) 0.01])
        if is_bad_chan(n)
            set(gca, 'XColor','r'); set(gca, 'YColor', 'r');
            set(gca, 'LineWidth', 3);
        end
        set(gca,'YTick', []); % set(gca, 'XTick',[]); 
        text((max(get(gca,'XLim')) - (max(get(gca,'XLim')) * 0.25)),(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.25)),num2str(n));
    end
end

%% Clear white matter electrodes
if clear_white_matter
    is_bad_chan = is_bad_chan | strcmpi(strrep(VNS_raw.anatomy_elecs,'Right','Left'),'Left-Cerebral-White-Matter');
end

end
