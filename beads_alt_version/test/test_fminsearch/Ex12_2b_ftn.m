function f = Ex12_2b_ftn(Y, xa, ya, A)

% compute f
x = Y(1);
y = Y(2);
f = (x-xa).^2 + (y-ya).^2 + A*(sin(y)*cos(x));

return
