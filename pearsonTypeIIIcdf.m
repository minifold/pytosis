function y=pearsonTypeIIIcdf(x,q)

a = q(2);% shape parameter
b = q(1);% scale parameter
delta1 = q(3);%

y = gamcdf(x-delta1,a,b);

end %end function