function pval = circ_kuiper_eq(alpha1,alpha2)
%% This function calculates pvalue of kuiper statistic of the empiircal
% cdf distributions derevied from input data alpha1 and alpha2
%
% This equation is derived from the work of M.A. Stephens
%   Stephens, M.A., (1970) "Use of the Kolmogorov-Smirnov, Cramer-Von Mises and
%         Related Statistics Without Extensive Tables", Journal of the Royal
%         Statistical Society. Series B, 32(1):115-122.
% and took some infrastructure from Matlab's kstest2
alpha1 = alpha1(:);
alpha2 = alpha2(:);

binEdges    =  [-inf ; sort([alpha1;alpha2]) ; inf];

binCounts1  =  histc (alpha1 , binEdges, 1);
binCounts2  =  histc (alpha2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

%
% Compute the test statistic of interest.
%

deltaCDF  =  sampleCDF2 - sampleCDF1;

Kuiper_stat   =  max(deltaCDF) + max(-1*deltaCDF);

%
% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
%

n1     =  length(alpha1);
n2     =  length(alpha2);
n      =  n1 * n2 /(n1 + n2);
lambda =  max((sqrt(n) + 0.155 + 0.24/sqrt(n)) * Kuiper_stat , 0);

%
%  Use the asymptotic Q-function to approximate the 2-sided P-value.
%
%   j       =  (1:101)';
%   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
%   pValue  =  min(max(pValue, 0), 1);

pval = (8*lambda.^2-2).*exp(-2*lambda.^2);
pval  =  min(max(pval, 0), 1);


%H  =  (alpha >= pValue);
end