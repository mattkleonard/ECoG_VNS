function kl_div = kl_div(X,Y,varargin)
% This function takes two 1d arrays of data X and Y 
% and finds the KL-divergence
% The binning follows a modification of the freedman-diaconis rule
%
% Variable inputs
% data_range - [roi_min, roi_max]

X = X(:);
Y = Y(:); % make nx1 arrays

%% Set bin parameters:
data_ampl = iqr([X; Y]); % amplitude of variation (set range to be 2*IQR)
data_center = median([X; Y]);
data_range = [data_center - data_ampl, data_center + data_ampl];
if length(varargin)>0
    if ~isempty(varargin{1})
        data_range = varargin{1};
    end
end

%% Set bins
n_bins = ceil(mean([length(X), length(Y)])^(1/3));
bin_edges = linspace(data_range(1),data_range(2), n_bins);

p = histcounts(X,bin_edges)/length(X);
q = histcounts(Y,bin_edges)/length(Y);

% Clear zeros
is_zero_val = (p==0) | (q==0);
p = p(~is_zero_val);
q = q(~is_zero_val);
kl_div = sum(p.*log2(p./q));

end