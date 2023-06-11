tic
n = 200;
A = 500;
a = zeros(1, n);
tic
for i = 1:n
    a(i)= max(abs(eig(rand(A))));
end
toc