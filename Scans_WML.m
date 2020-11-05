function no=Scans_WML(samples,ch,n,k,L,r,u)
% -------------------------------------------------------------------------
%
%     WML = Scans_WML(samples,ch,n,k,L,r,u)
%           This function returns the number of trials required to get u scanning
%           windows of length k containing r successes (occurence of 'ch')
%           which may have overlapping part of length L with the previous
%           scanning window that has been counted.
%
%     Example
%       Input:
%           sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%           statistic = Scans_ML(sample_seq,1,16,4,1,2)
%           stat_wait = Scans_WML(sample_seq,1,16,4,1,2,3)
%
%       Output:
%           statistic = 4
%           stat_wait = 10
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.

if nargin<6
    error('Scans_WML:TooFewInputs','Requires at least six input arguments.');
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
                j=j+(k(h)-1-L)*(temp1>=r(h));
            end
            if const<u(h)
                no(i,h)=(n+1);
            end
        end
    end
end