function D=Coef(n,i,A,B)
if n==1
    if i<=1
        D=A*(i==0)+B*(i==1);
    else
        D=A*0;
    end
else
    if (i>n)||(i<0)
        D=A*0;
    else
        D=Coef(n-1,i,A,B)*A+Coef(n-1,i-1,A,B)*B;
    end
end
