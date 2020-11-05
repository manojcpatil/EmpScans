function [Stats]=describe(N)
%
% [Stats]=describe(N)
%
% This function returns descriptive Statistics and empirical CDF
%      for given Scan Statistics
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.

[row,col]=size(N);
Stats.Mean=mean(N);
Stats.Mode=mode(N);
Stats.Median=median(N);
Stats.Min=min(N);
Stats.Max=max(N);
Stats.Range=range(N);
Stats.Stdev=std(N);
Stats.Variance=var(N);
if min(row,col)==1
    t=tabulate(N);
    Stats.Empirical_CDF=[t(:,1) t(:,2)/sum(t(:,2)) cumsum(t(:,3))./100];
    hold on
    Stats.Diag=cdfplot(N); grid  off; xlabel('x values'); ylabel('F(x)'); title('Empirical Cumulative Probability distribution Function'); 
    hold off
else
    Stats.Variance=cov(N);
    Stats.Correlation=corr(N);
    Stats.Mat=plotmatrix(N); title('Matrix Plot for the Scan Statistics under study')
end