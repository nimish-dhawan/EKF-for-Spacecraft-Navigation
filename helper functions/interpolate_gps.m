%
function [y,dy] = interpolate_gps(xa,ya,x,n)

% inputs xa,ya
% INPUT
% ya    x,y or z coordinate
% xa    time at each 15 minute interval
% x     time at which you want the position
% n     order of interpolation
%
% OUTPUT
% y     new coordinate at desired time
% dy    error estimate of this value


nmax = 24*60/5+1;     % clock corrections given every 5 min. 
ns = 1;
% where do you want the interpolation
diff = abs(x - xa(1));

for i = 1:n
    difft = abs(x-xa(i));
    if difft < diff
        ns = i;
        diff = difft;
    end
    c(i) = ya(i);
    d(i) = ya(i);
end

% initial approximation to y
y = ya(ns);
ns = ns-1;

for m = 1:n-1
    for i = 1:n-m
        ho = xa(i)-x;
        hp = xa(i+m)-x;
        w  = c(i+1)-d(i);
        den = ho-hp;
        if den == 0
            'failure'
        end
        den = w/den;
        d(i) = hp*den;
        c(i) = ho*den;
    end
    if 2*ns < n-m
        dy = c(ns+1);
    else
        dy = d(ns);
        ns = ns-1;
    end
    y = y+dy;
end
