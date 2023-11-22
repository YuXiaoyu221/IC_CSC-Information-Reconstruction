function [dispersion,meanw] = covw(x,w)

% [covaw meanw] = covw(x,w,varargin)
%
%   COVW Covariance matrix with weights.
%
%   COVW(W,X), if X is a vector, returns the variance.  For matrices,
%   where each row is an observation, and each column a variable,
%   COVW(X,W) is the variance-covariance matrix.  DIAG(COVW(X,W)) is
%   a vector of variances for each column, and SQRT(DIAG(COVW(X,W)))
%   is a vector of standard deviations. 
%
%   Observations are weighted with vector of weights, W, which has all
%   nonnegative values.  COVW(X,W) gives the weighted estimate of the
%   variance-covariance matrix.
%   
%   COVW(X,W) normalizes by (N-1) where N is the number of
%   observations.
%
%   COVW(X,W,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.
%   COVW(X,W,0) is the same as COVW(X,W).
%
%   The weighted mean is removed from each column before calculating the result.

if(size(w,2)~=1)
    w = w(:); % w is now a column vector
end
[m,n] = size(x);

if m==1  % Handle special case
    dispersion = 0;
else
    sumw = sum(w);
    meanw = sum(repmat(w,1,n).*x) / sumw;%B = repmat(A,[m n p...])
    %xc = x - repmat(sum(repmat(w,1,n).*x),m,1) / sumw;  % Remove weighted mean
    xc = x - repmat(meanw,m,1);  % Remove weighted mean
    clear x
    % direct method uses way too much memory:
    %dispersion = xc' * diag(w) * xc / sumw;
    c = repmat(sqrt(w),1,n) .* xc;
    dispersion = xc' * xc / sumw;
    clear xc;    
    dispersion = m/(m-1) * dispersion;
end