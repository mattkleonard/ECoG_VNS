function corcoeffs_X = circ_corcoeff(X)
% this function takes an nxm matrix X of circular data (ranging [0, 2pi),
% and returns an mxm matrix where the element i,j is equal to the circular 
% correlation coefficient between X(:,i) and X(:,j)
%   Input:
%     X             Matrix of circular data
%
%   Output:
%   corcoeffs_X     Pairwise Circular-Correlation Matrix for columns of X
%
%   This function was inspired by
%   Circular Satsitics Toolbox - Philipp Berens, 2009
%
%   Created by Ben Lucas - September 23, 2016

X_bar = angle(mean(exp(sqrt(-1)*X),1));     % Circular mean of each column
X_circ = sin(gsubtract(X,X_bar));   % sine transform of X - X_bar

num = X_circ'*X_circ;
X_circ_sqr = sum(X_circ.^2,1);
den = sqrt(X_circ_sqr'*X_circ_sqr);

corcoeffs_X = num./den;

end