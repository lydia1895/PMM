

%{
x = linspace(1,2,20);
y = 2*x.^2+1;
polyfit(x,y,2)
%}

n= 140
XIn = linspace(1,2,n)
YIn = linspace(1,2,n)
nx = size(XIn)
ny=size(YIn)
ZIn = zeros(n,n);
nz = size(ZIn)

for i = 1:n
    for j = 1:n
        %ii = i
        %jj = j
        ZIn(i,j) = XIn(i)^2+3*YIn(j);
    end
end

[XOut, YOut, ZOut] = prepareSurfaceData(XIn, YIn, ZIn)

g = fit([XOut, YOut],ZOut,'poly33')
g.p00