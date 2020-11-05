function no=Scans_WN2(samples,ch,n,k,r,u)
% -------------------------------------------------------------------------
%     WN = Scans_WN2(samples,ch,n,k,r,u)
%           This function returns the number of trials required to get u
%           scanning windows of length k containing r successes (occurence
%           of 'ch') where recounting starts only after the completion of
%           k-tuple of successive trials with at least r successes.
%
%       Example:
%       Input
%         sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%         statistic  = Scans_N2(sample_seq,1,16,4,2)
%         stat_wait  = Scans_WN2(sample_seq,1,16,4,2,3)
%
%       Output
%         statistic  = 3
%         stat_wait  = 14
%
%       Reference:
%           Balakrishnan, N. and Koutras M. V. (2002) Runs and scans with
%           applications, A Wiley-Interscience Publication, John Wiley &
%           Sons, Inc.
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
if nargin<6
    error('Scans_WN2:TooFewInputs','Requires at least six input arguments.');
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
                const=const+(temp1>=r(h)); no(i,h)=j;
                j=j+(k(h)-1)*(temp1>=r(h));
            end
            if const<u(h)
                no(i,h)=(n+1);
            end
        end
    end
end