%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Example MATLAB scripts for using EmpScans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------
%% Non-overlapping counting first kind
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_N1(sample_seq,6,16,4,2)
stat_wait  = Scans_WN1(sample_seq,6,16,4,2,3)
% 
% Output
%statistic  = 4
%stat_wait  = 11
% 
%------------------------------------------------------------
%% Non-overlapping counting Second kind
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_N2(sample_seq,6,16,4,2)
stat_wait  = Scans_WN2(sample_seq,6,16,4,2,3)
% 
% Output
% statistic  = 3
% stat_wait  = 14
% 
%------------------------------------------------------------
%% Overlapping counting
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_M(sample_seq,6,16,4,2)
stat_wait  = Scans_WM(sample_seq,6,16,4,2,7)
% 
% Output
% statistic  = 10
% stat_wait  = 11
% 
%------------------------------------------------------------

%% l-Overlapping counting
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_ML(sample_seq,6,16,4,1,2)
stat_wait  = Scans_WML(sample_seq,6,16,4,1,2,3)
% 
% Output
% statistic  = 4
% stat_wait  = 10
% 
%------------------------------------------------------------
%% Exact counting
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_E(sample_seq,6,16,4,2)
stat_wait  = Scans_WE(sample_seq,6,16,4,2,3)
% 
% Output
% statistic  = 3
% stat_wait  = 14
% 
%------------------------------------------------------------
%% Other Scan Statistic (L_{n,k}^{(i)})}
% Input 
sample_seq = [6,5,6,9,1,6,6,6,8,4,6,3,2,6,6,2];
statistic  = Scans_L(sample_seq,6,16,4)
% 
% Output
% statistic  = 3
% 
%------------------------------------------------------------

%% Joint Distributions of Scan Statistics}
%% Same counting schemes and number of successes with different scanning window sizes (k_1,k_2,..., k_s)}
% Input 
sample_seq = [6,5,6,9,6,6,6,6,8,4,6,3,6,6,6,2,3,6,6,5];
K = [3,5,7];
D = Scans_M(sample_seq,6,20,K,3)
% 
% Output
% D =
%  3 13 14
% 
%------------------------------------------------------------
% Five sample stochastic sequences of length n=20
samples=[2,1,3,3,2,1,2,2,2,2,2,2,2,1,2,2,2,2,1,2
3,2,2,2,2,2,3,2,2,2,2,2,2,2,2,3,1,3,3,2
2,2,1,2,1,2,2,2,2,3,3,2,2,1,1,2,2,2,3,2
3,2,1,2,3,1,2,1,3,2,2,2,1,2,1,1,2,1,3,3
2,2,3,2,3,2,2,3,1,1,1,2,3,3,1,2,2,2,2,2];
%------------------------------------------------------------
% Calculation of (M^(2)_{20,3,3},M^(2)_{20,4,3}, M^(2)_{20,5,3}) for samples
% Input 
k = [3,4,5];
D = Scans_M(samples,2,20,k,3)
% 
% Output
% D =
%  7 13 13
%  9 13 13
%  3 8 14
%  1 3 5
%  3 5 7
% 
%------------------------------------------------------------
%% Same counting schemes and window lengths but with different number of successes ($r_1, r_2,. . ., r_s $)

% Calculation of (M^(2)_{20,4,2},M^(2)_{20,4,3}, M^(2)_{20,4,4}) for samples
% Input 
 r = [2,3,4];
 D = Scans_M(samples,2,20,4,r)
% 
% Output
% D =
% 14 13 5
% 14 13 7
% 17 8 1
% 10 3 0
% 10 5 2
% 
%------------------------------------------------------------

%% Different counting schemes, window lengths and number of successes

% Calculation of (N^(1)(2)_{20,3,2},N^(2)(2)_{20,4,3}, E^(2)_{20,5,4}) for samples
% Input 
D(:,1) = Scans_N1(samples,2,20,3,2);
D(:,2) = Scans_N2(samples,2,20,4,3);
D(:,3) = Scans_E(samples,2,20,5,4);
D
% 
% Output
% D =
%  6 4 3
%  6 4 3
%  6 3 2
%  3 1 1
%  4 2 1
% 
%------------------------------------------------------------
%% Joint distribution of scan statistics for different states

% Calculation of (N^(1)(1)_{20,5,3},N^(1)(2)_{20,5,3}) for samples
% Input 
D(:,1) = Scans_N1(samples,1,20,5,3);
D(:,2) = Scans_N1(samples,2,20,5,3);
D
% 
% Output
%  0 4
%  0 4
%  0 4
%  1 1
%  1 2
% 
%------------------------------------------------------------
%% Descriptive Statistics

% Input 
samples = binornd(1,0.6,1000,20); % Simulate sample sequences
N = Scans_N1(samples,1,20,5,3);   % Calculate N^(1)_{20,5,3}
D = describe(N)                   % Descriptive Statistics  
% 
%------------------------------------------------------------
%------------------------------------------------------------
% Output
%  Mean: 3.4300
%  Mode: 4
%Median: 4
%   Min: 0
%   Max: 6
% Range: 6
% Stdev: 0.9829
%  Variance: 0.9661
% Empirical_CDF: [7x3 double]
%  Diag: 159.0021
% 
%------------------------------------------------------------

%------------------------------------------------------------
% Input 
D.Empirical_CDF
% 
% Output
%  00.00100.0010
% 1.00000.02900.0300
% 2.00000.13600.1660
% 3.00000.33300.4990
% 4.00000.38000.8790
% 5.00000.11600.9950
% 6.00000.00501.0000
% 
%------------------------------------------------------------

%------------------------------------------------------------
% Input 
E = Scans_E(samples,1,20,5,3);
V = var(E)
% 
% Output
% V =
% 0.5364
% 
%------------------------------------------------------------

%------------------------------------------------------------
% Input 
t = tabulate(E);
Emp_PMF=[t(:,1) t(:,2)./sum(t(:,2))] % Estimation of PMF
Emp_CDF=[t(:,1) cumsum(t(:,3))./100] % Estimation of CDF
% 
% Output
% Emp_PMF =
%  0        0.0120
% 1.0000    0.1110
% 2.0000    0.4460
% 3.0000    0.4170
% 4.0000    0.0140
% Emp_CDF =
%  0        0.0120
% 1.0000    0.1230
% 2.0000    0.5690
% 3.0000    0.9860
% 4.0000    1.0000
% 
%------------------------------------------------------------
%% Estimation of Joint Probabilities
% Input 
size = 100000; p = 0.7; n = 25; k = [3 5]; r = 2;
samples = binornd(1,p,size,n);    % Generating Sample stochastic sequences
E = Scans_E(samples,1,n,k,r);    
[TABLE, Chi2, P, RowCol_Labels] = crosstab(E(:,1),E(:,2));
RowCol_Labels    % Row and Column labels for Joint Probability table
pE=TABLE./size   % Empirical joint probabilities
plotmatrix(E,'o')
title('Matrix plot for (E^{(1)}_{25,3,2},E^{(1)}_{25,5,2})')
%------------------------------------------------------------
%% Application
% Input 
% X=(Y<8); if daily minimum temperature data is available and stored in Y.
X = binornd(1,0.3,50,120);
D = Scans_N2(X,1,120,10,7);
tabulate(D)
% 
% Output
%  ValueCount   Percent
% 0       31     61.00%
% 1       16     32.00%
% 2        2      4.00%
% 3        1      2.00%% 
%------------------------------------------------------------
%% Exact probability distribution of Scan Statistics in one order Markov
%% chains

%% Non-overlapping counting first kind
% Input 
n=20; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_N1(ch,n,k,r,m,P)
%
% Output
% Val_Prob =
%          0    0.0000
%     1.0000    0.0001
%     2.0000    0.0010
%     3.0000    0.0064
%     4.0000    0.0276
%     5.0000    0.0857
%     6.0000    0.1900
%     7.0000    0.2900
%     8.0000    0.2748
%     9.0000    0.1158
%    10.0000    0.0086
%------------------------------------------------------------
%% Non-overlapping counting Second kind
% Input 
n=20; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_N2(ch,n,k,r,m,P)
%
% Output
% Val_Prob =
%          0    0.0000
%     1.0000    0.0002
%     2.0000    0.0032
%     3.0000    0.0419
%     4.0000    0.3412
%     5.0000    0.6135
%------------------------------------------------------------
%% Overlapping counting
% Input 
n=20; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_M(ch,n,k,r,m,P)
%
% Output
% Val_Prob =
%          0    0.0000
%     1.0000    0.0000
%     2.0000    0.0000
%     3.0000    0.0001
%     4.0000    0.0002
%     5.0000    0.0004
%     6.0000    0.0009
%     7.0000    0.0018
%     8.0000    0.0037
%     9.0000    0.0072
%    10.0000    0.0132
%    11.0000    0.0235
%    12.0000    0.0393
%    13.0000    0.0674
%    14.0000    0.1010
%    15.0000    0.1464
%    16.0000    0.1161
%    17.0000    0.4790
%% L-Overlapping counting
% Input 
n=20; k=4; r=2; m=1; ch=1; L=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_ML(ch,n,k,L,r,m,P)
%
% Output
% Val_Prob =
%          0    0.0000
%     1.0000    0.0001
%     2.0000    0.0008
%     3.0000    0.0067
%     4.0000    0.0449
%     5.0000    0.2335
%     6.0000    0.7139
%% Exact Counting
% Input 
n=20; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_E(ch,n,k,r,m,P)
%
% Output
% Val_Prob =
%          0    0.1322
%     1.0000    0.3190
%     2.0000    0.3437
%     3.0000    0.1747
%     4.0000    0.0301
%     5.0000    0.0004



%% Waiting time distribution in Exact Counting
% Input 
u=2; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_WE(ch,u,k,r,m,P)
% 
% Val_Prob =
% 
%     8.0000    0.0449
%     9.0000    0.0449
%    10.0000    0.0486
%    11.0000    0.0488
%    12.0000    0.0481
%    13.0000    0.0461
%    14.0000    0.0441
%    15.0000    0.0421
%    16.0000    0.0401
%    17.0000    0.0381
%    18.0000    0.0362
%    19.0000    0.0344
%    20.0000    0.0325
%    21.0000    0.0307
%    22.0000    0.0290
%    23.0000    0.0273
%    24.0000    0.0256
%    25.0000    0.0241
%    26.0000    0.0226
%    27.0000    0.0211
%    28.0000    0.0198
%    29.0000    0.0185
%    30.0000    0.0172
%    31.0000    0.0161
%    32.0000    0.0150
%    33.0000    0.0139
%    34.0000    0.0130
%    35.0000    0.0120
%    36.0000    0.0112
%    37.0000    0.0104
%    38.0000    0.0096
%    39.0000    0.0089
%    40.0000    0.0083
%    41.0000    0.0076
%    42.0000    0.0071
%    43.0000    0.0065
%    44.0000    0.0060
%    45.0000    0.0056
%    46.0000    0.0051
%    47.0000    0.0048
%    48.0000    0.0044
%    49.0000    0.0040
%    50.0000    0.0037
%    51.0000    0.0034
%    52.0000    0.0032
%    53.0000    0.0029
%    54.0000    0.0027
%    55.0000    0.0025
%    56.0000    0.0023
%    57.0000    0.0021
%    58.0000    0.0019
%    59.0000    0.0018
%    60.0000    0.0016
%    61.0000    0.0015
%    62.0000    0.0014
%    63.0000    0.0013
%    64.0000    0.0012
%    65.0000    0.0011
%    66.0000    0.0010
%    67.0000    0.0009
%% Waiting time distribution in Nonoverlaping Counting first kind
% Input 
u=2; k=4; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_WN1(ch,u,k,r,m,P)
% 
% Val_Prob =
% 
%     4.0000    0.3072
%     5.0000    0.2611
%     6.0000    0.1805
%     7.0000    0.0972
%     8.0000    0.0604
%     9.0000    0.0378
%    10.0000    0.0235
%    11.0000    0.0137
%    12.0000    0.0079
%    13.0000    0.0046
%    14.0000    0.0026
%    15.0000    0.0015
%    16.0000    0.0009
%    17.0000    0.0005
%    18.0000    0.0003
%    19.0000    0.0002
%    20.0000    0.0001
%    21.0000    0.0000
%    22.0000    0.0000
%    23.0000    0.0000
%    24.0000    0.0000
%    25.0000    0.0000
%    26.0000    0.0000
%    27.0000    0.0000
%    28.0000    0.0000
%    29.0000    0.0000
%    30.0000    0.0000
%    31.0000    0.0000
%    32.0000    0.0000
%% Waiting time distribution in Nonoverlaping Counting Second kind
% Input 
u=2; k=4; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_WN2(ch,u,k,r,m,P)
% Output
% Val_Prob =
% 
%     8.0000    0.8034
%     9.0000    0.0778
%    10.0000    0.0513
%    11.0000    0.0302
%    12.0000    0.0172
%    13.0000    0.0090
%    14.0000    0.0050
%    15.0000    0.0028
%    16.0000    0.0015
%    17.0000    0.0008
%    18.0000    0.0005
%    19.0000    0.0003
%    20.0000    0.0001
%    21.0000    0.0001
   
%% Waiting time distribution in Overlaping Counting Scheme
% Input 
u=2; k=4; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_WM(ch,u,k,r,m,P)
%
% Output
% Val_Prob =
% 
%     5.0000    0.8515
%     6.0000    0.0584
%     7.0000    0.0399
%     8.0000    0.0233
%     9.0000    0.0128
%    10.0000    0.0065
%    11.0000    0.0035
%    12.0000    0.0019
%    13.0000    0.0010
%    14.0000    0.0006
%    15.0000    0.0003
%    16.0000    0.0002
%    17.0000    0.0001
%% Waiting time distribution in L-Overlaping Counting Scheme
% Input 
u=2; k=4; L=1; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_WML(ch,u,k,L,r,m,P)
%
%     7.0000    0.8127
%     8.0000    0.0746
%     9.0000    0.0497
%    10.0000    0.0288
%    11.0000    0.0161
%    12.0000    0.0083
%    13.0000    0.0045
%    14.0000    0.0025
%    15.0000    0.0014
%    16.0000    0.0007
%    17.0000    0.0004
%    18.0000    0.0002
%    19.0000    0.0001
%    20.0000    0.0001

n=500; k=4; r=2; m=1; ch=1;
P=[0.4,0.6;0.2,0.8];
Val_Prob=PrScan_N1(ch,n,k,r,m,P)










% clear all; clc;
% n=16;   states=[0 1];
% alpha=[1 0];
% TPM=[0.2,0.8;0.4,0.6];
% samples=10000;
% 
% 
% tic
% y=GenSeq(n+1,states,alpha,TPM,samples);
% y1=y(:,[2:n+1]);
% toc
% save data1 y1;
% % elapsed time around 675 sec
% 
% load data1.mat
% % tic
% %  dN1=Scans_N1(y1,1,10,4,2);
% %  tN1=tabulate(dN1);
% %  pN1=[tN1(:,1) tN1(:,3)./100]
% % toc
% % tic
% %  dN2=Scans_N2(y1,1,10,4,2);
% %  tN2=tabulate(dN2);
% %  pN2=[tN2(:,1) tN2(:,3)./100]
% % toc
% % tic
% %  dM=Scans_M(y1,1,10,4,2);
% %  tM=tabulate(dM);
% %  pM=[tM(:,1) tM(:,3)./100]
% % toc
% % tic
% %  dML=Scans_ML(y1,1,10,4,2,2);
% %  tML=tabulate(dML);
% %  pML=[tML(:,1) tML(:,3)./100]
% % toc
% % tic
% %  dE=Scans_E(y1,1,10,4,2);
% %  tE=tabulate(dE);
% %  pE=[tE(:,1) tE(:,3)./100]
% % toc
% 
% 
% 
% % Elapsed time is 437.218040 seconds.
% % 
% % pN1 =
% % 
% %          0    0.0003
% %     1.0000    0.0135
% %     2.0000    0.1695
% %     3.0000    0.5724
% %     4.0000    0.2380
% %     5.0000    0.0064
% % 
% % Elapsed time is 21.453988 seconds.
% % 
% % pN2 =
% % 
% %          0    0.0003
% %     1.0000    0.0270
% %     2.0000    0.9727
% % 
% % Elapsed time is 18.245398 seconds.
% % 
% % pM =
% % 
% %          0    0.0003
% %     1.0000    0.0013
% %     2.0000    0.0041
% %     3.0000    0.0130
% %     4.0000    0.0345
% %     5.0000    0.0669
% %     6.0000    0.1475
% %     7.0000    0.7324
% % 
% % Elapsed time is 17.985588 seconds.
% % 
% % pML =
% % 
% %          0    0.0003
% %     1.0000    0.0037
% %     2.0000    0.0333
% %     3.0000    0.1678
% %     4.0000    0.7950
% % 
% % Elapsed time is 18.515276 seconds.
% % 
% % pE =
% % 
% %          0    0.1680
% %     1.0000    0.5197
% %     2.0000    0.3122
% % 
% % Elapsed time is 18.177430 seconds.
% 
% 
% tic
% dN1=Scans_N1(y1,1,10,4,[2,3]);
% pN1=crosstab(dN1(:,1),dN1(:,2))./samples
% toc
% 
% tic
% dN2=Scans_N2(y1,1,10,4,[2,3]);
% pN2=crosstab(dN2(:,1),dN2(:,2))./samples
% toc
% 
% tic
% dM=Scans_M(y1,1,10,4,[2,3]);
% pM=crosstab(dM(:,1),dM(:,2))./samples
% toc
% 
% tic
% dML=Scans_ML(y1,1,10,4,2,[2,3]);
% pML=crosstab(dML(:,1),dML(:,2))./samples
% toc
% 
% tic
% dE=Scans_E(y1,1,10,4,[2,3]);
% pE=crosstab(dE(:,1),dE(:,2))./samples
% toc
% 
% 
% % OUTPUT
% % pN1 =
% % 
% %     0.0003         0         0         0
% %     0.0079    0.0056         0         0
% %     0.0393    0.1265    0.0037         0
% %     0.0005    0.1977    0.3741         0
% %          0         0    0.1832    0.0547
% %          0         0         0    0.0064
% % 
% % Elapsed time is 28.038621 seconds.
% % 
% % pN2 =
% % 
% %     0.0003         0         0
% %     0.0079    0.0191         0
% %     0.0398    0.3342    0.5987
% % 
% % Elapsed time is 21.459502 seconds.
% % 
% % pM =
% % 
% %     0.0003         0         0         0         0         0         0         0
% %     0.0013         0         0         0         0         0         0         0
% %     0.0031    0.0010         0         0         0         0         0         0
% %     0.0060    0.0046    0.0024         0         0         0         0         0
% %     0.0089    0.0097    0.0097    0.0062         0         0         0         0
% %     0.0108    0.0153    0.0170    0.0159    0.0079         0         0         0
% %     0.0108    0.0231    0.0318    0.0419    0.0256    0.0144         0         0
% %     0.0069    0.0235    0.0469    0.0871    0.1220    0.1437    0.1346    0.1678
% % 
% % Elapsed time is 22.458893 seconds.
% % 
% % pML =
% % 
% %     0.0003         0         0         0         0
% %     0.0028    0.0010         0         0         0
% %     0.0094    0.0148    0.0091         0         0
% %     0.0205    0.0508    0.0796    0.0170         0
% %     0.0151    0.0629    0.1995    0.3175    0.2000
% % 
% % Elapsed time is 20.617286 seconds.
% % 
% % pE =
% % 
% %     0.0066    0.0547    0.1066
% %     0.0079    0.1857    0.3261
% %     0.0398    0.2079    0.0645
% % 
% % Elapsed time is 20.226440 seconds.
% 
% 
% dN2=Scans_N2(y1,1,10,4,[2,3]);
% [TABLE, CHI2, P, LABELS]=crosstab(dN2(:,1),dN2(:,2));
% TABLE
% LABELS
% % OUTPUT
% % TABLE =
% %            1           1           0           0           0
% %            4          42          42           0           0
% %           53         406        1132         790           0
% %           12         261        1793        4399        1064
% % LABELS = 
% %     '1'    '0'
% %     '2'    '1'
% %     '3'    '2'
% %     '4'    '3'
% %      []    '4'