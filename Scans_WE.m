function no=Scans_WE(samples,ch,n,k,r,u)
% -------------------------------------------------------------------------
%
%     WE = Scans_WE(samples,ch,n,k,r,u)
%           This function returns the number of trials required to get u
%           scanning windows of length k containing exactly r successes
%           (occurrence of 'ch') where recounting starts only after the
%           completion of k-tuple of successive trials containing exactly
%           r occurrences of state ‘ch’.
%
%       Example:
%           Input:
%               sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%               statistic = Scans_E(sample_seq,1,16,4,2)
%               stat_wait = Scans_WE(sample_seq,1,16,4,2,3)
%
%           Output:
%               statistic = 3
%               stat_wait = 14
%
%       Reference:
%       Mood, A. M. (1940) The Distribution Theory of Runs, Annals of
%       Mathematical Statistics, 11, 367–392.
%
%       See also:
%           Scans_N1, Scans_N2, Scans_M, Scans_E and Scans_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
% -------------------------------------------------------------------------
if nargin<6
    error('Scans_WE:TooFewInputs','Requires at least six input arguments.');
else
    [m,n]=size(samples);
    f=length(r);
    for h=1:f
        samples1=(samples==ch(h));
        for i=1:m
            temp=samples1(i,:);
            j=k-1;    const=0;
            while const<u(h)&&j<n
                j=j+1;
                temp1=sum(temp(j-k(h)+1:j));
                const=const+(temp1==r(h)); no(i,h)=j;
                if temp1==r(h)    % exactly r successes in k tuples
                    j=j+k(h)-1;
                end
            end
            if const<u(h)
                no(i,h)=(n+1);
            end
        end
    end
end