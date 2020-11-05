% EmpScans Toolbox
%
%   Toolbox for Scan statistics with MATLAB
%
%   Scans_N1 	-#windows of length at most k containing r successes (i.e. occurrence of state `ch') where recounting starts from scratch each time  r  successes encountered.
%   Scans_N2	-#windows of length k containing r successes where recounting starts only after the completion of k-tuple of successive trials with at least r successes.
%   Scans_M     -#overlapping k-tuples which contain at least r successes.
%   Scans_ML	-#l-overlapping ( 0<=l<=k-1 ) k-tuples containing at least r successes.
%   Scans_E     -#k-tuples containing exactly r successes where recounting starts only after the completion of k-tuple of successive trials containing exactly r successes.
%   Scans_L     -Maximum number of successes in a scanning window of length k moving along the sequence of n
%   Scans_WN1	-Waiting number of trials for which Scans_N1(s,ch,n,k,r) = u
%   Scans_WN2	-Waiting number of trials for which Scans_N2(s,ch,n,k,r) = u
%   Scans_WM	-Waiting number of trials for which Scans_M(s,ch,n,k,r) = u
%   Scans_WML	-Waiting number of trials for which Scans_ML(s,ch,n,k,L,r) = u
%   Scans_WE	-Waiting number of trials for which Scans_E(s,ch,n,k,r) = u
%   describe    -Descriptive statistics of Scan statistics in vector/matrix N
%
% Authors:   
%   Patil Manoj C. and Shinde Ramakrishna L., 2019.