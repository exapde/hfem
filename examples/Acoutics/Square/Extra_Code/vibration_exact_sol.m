num = 100;
x = linspace(0,1,num);
y = linspace(0,1,num);
t = 1;
u = zeros(num,num);
v = zeros(num,num);

for i = 1:num
    for j = 1:num
        u(i,j) = (1/(sqrt(2)*pi))*sin(pi*x(i))*sin(pi*y(j))*sin(sqrt(2)*pi*t);
        v(i,j) = sin(pi*x(i))*sin(pi*y(j))*cos(sqrt(2)*pi*t);
    end
end

figure(2)
surf(x,y,u)

