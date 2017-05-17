function tstat = t_welch(X1,X2,varargin)
%% Calculates welch's estimate of the t-statistic for two data sets x,y
% if 'dim' is specified, it will calculate the tstatistic along that
% dimension.
% Inputs:
% X1,X2 - two matricies of data (may be unequal along the dimension tstat
% is calcuated upon
%
% Variable input:
% dimension - n<ndims(X1, X2) - dimension to take the tstat along
xdims = size(X1);
if (length(xdims) == 2) & (xdims(1) == 1 | xdims(2)==1)
    % 1 D data:
    tstat = (nanmean(X1) - nanmean(X2))/sqrt(nanvar(X1)/length(X1)+nanvar(X2)/length(X2));
else
    % Determine the direction 
    dim = 1;
    if ~isempty(varargin);
        if ~isempty(varargin{1});
            dim = varargin{1};
        end
    end
    tstat = (nanmean(X1,dim) - nanmean(X2,dim))./sqrt(nanvar(X1,[],dim)./size(X1,dim)+nanvar(X2,[],dim)./size(X2,dim));
end



end

