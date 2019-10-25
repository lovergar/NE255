function [out] = atanshift(x)

out=2*pi*(x<=0)+atan(x);
 %out=2*pi*(atan(x)<=0)+atan(x)
end