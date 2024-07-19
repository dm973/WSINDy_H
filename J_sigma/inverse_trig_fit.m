f = @(x) 1./(1+0.9*cos(x));
M = 200;
x = linspace(0,4*pi,M)';
y = f(x);
plot(x,y)

fq = 1:16;
L = length(fq);
A = ones(M,L*2+1);
for j=2:L
    A(:,j) = cos(fq(j-1)*x);
    A(:,j+L) = sin(fq(j-1)*x);
end

w = A \ y