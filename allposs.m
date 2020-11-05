function mat=allposs(r, states)
% mat=allposs(r) this function returns the all possible permutations of
% 1,2,3,4 of length r+1
% which is useful for checking the rth order markov dependency
if nargin<2
    error('allposs:TooFewInputs','Requires at least two input arguments.');
elseif r<1
    error('allposs:LengthError','Length should be grater than 1.');
else
    [arow,bcol]=size(states);
    if arow==1
        states=states';
    end
    mat=states;
    while r>1
        [m,n]=size(mat);
        diff=[mat([1:m],:) states(ones(m,1))];
        for i=2:length(states)
            diff([(i-1)*m+1:i*m],:)=[mat([1:m],:) states(i*ones(m,1))];
        end
        mat=diff;
        r=r-1;
    end
end