function no=Scans_L(samples,ch,n,k)
% -------------------------------------------------------------------------
%
%       L = Scans_L(samples,ch,n,k)
%           This function returns the maximum number of successes
%           (occurence of 'ch') in scanning windows of length k
%
%       Example
%       Input:
%           sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%           statistic = Scans_L(sample_seq,1,16,4)
%
%       Output:
%           statistic = 3
%
%       Reference:
%       Glaz, J. and Naus, J. I. (1991) Tight bounds and Approximations for
%       scan statistic probabilities for discrete data, Annals of Applied
%       Probability, 1, 306–318.
%
%       See also:
%           Scans_N1, Scans_N2, Scans_M, Scans_E and Scans_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
if nargin<4
    error('Scans_WN1:TooFewInputs','Requires at least four input arguments.');
else
    [m,n]=size(samples);
    f=length(k);
    for h=1:f
        samples1=(samples==ch(h));
        for i=1:m
            temp1=samples1(i,:)';
            for j=0:n-k(h)
                temp(j+1)=sum(double(temp1(j+1:j+k(h))==1));
            end
            no(i,h)=max(temp);clear temp
        end
    end
end