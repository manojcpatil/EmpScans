function no=Scans_WN1(samples,ch,n,k,r,u)
% -------------------------------------------------------------------------
%     WN = Scans_WN1(samples,ch,n,k,r,u)
%           This function returns the number of trials required to get u
%           scanning windows of length at most k containing r successes
%           (occurence of 'ch') where recounting starts from scratch each
%           time r successes encountered.
%
%       Example:
%       Input:
%           sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%           statistic  = Scans_N1(sample_seq,1,16,4,2)
%           stat_wait  = Scans_WN1(sample_seq,1,16,4,2,3)
%
%       Output
%           statistic  = 4
%           stat_wait  = 11
%
%       Reference:
%           Balakrishnan, N. and Koutras M. V. (2002) Runs and scans with
%           applications, A Wiley-Interscience Publication, John Wiley &
%           Sons, Inc.
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
if nargin<6
    error('Scans_WN1:TooFewInputs','Requires at least six input arguments.');
else
    [m,n]=size(samples);
    f=length(r);
    for h=1:f
        samples1=(samples==ch(h));
        for i=1:m
            temp=[zeros(k(h)-r(h),1);samples1(i,:)'];
            const=0;j=0;
            while const<u(h)&&j<=(n-k(h))
                j=j+1;
                temp1=sum(temp((j:j+k(h)-1)));
                const=const+(temp1==r(h));no(i,h)=j+k(h)-r(h)-1;
                if temp1==r(h)
                    temp(j:j+k(h)-1)=0;    j=j+k(h)-r(h)-1;
                    % As recounting starts only after completion of r success in
                    % atmost k tuples
                end
            end
            if const<u(h)
                no(i,h)=(n+1);
            end
        end
    end
end