function an = LagrangePolynomialCoefficient(tn, t)

k = length(tn);
an = zeros(1,k);
for i = 1:k
    an(i) = 1;
    for j = 1:k
        if (i ~= j)
            an(i) = an(i)*((t-tn(j))/(tn(i)-tn(j)));
        end
    end
end

return;

t(1) = 0.435866521500000;
t(2) = 0.717933260750000;
t(3) = 1.000000000000000;
dt = ([t(1)-0 t(2:end)-t(1:end-1)]);

% Stage 1
tt = t(1);
% linear
tn1 = [0-dt(3) 0];
an1 = LagrangePolynomialCoefficient(tn1, tt);
% quadratic
tn2 = [0-dt(2)-dt(3) 0-dt(3) 0];
an2 = LagrangePolynomialCoefficient(tn2, tt);
% cubic
tn3 = [0-sum(dt) 0-dt(2)-dt(3) 0-dt(3) 0];
an3 = LagrangePolynomialCoefficient(tn3, tt);
% P4
tn4 = [-sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0];
an4 = LagrangePolynomialCoefficient(tn4, tt);
% P5
tn5 = [-sum(dt)-dt(2)-dt(3) -sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0];
an5 = LagrangePolynomialCoefficient(tn5, tt);
% P6
tn6 = [-2*sum(dt) -sum(dt)-dt(2)-dt(3) -sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0];
an6 = LagrangePolynomialCoefficient(tn6, tt);
An = [0 0 0 0 0 an1; 0 0 0 0 an2; 0 0 0 an3; 0 0 an4; 0 an5; an6];

% Stage 2
tt = t(2);
% linear
tn1 = [0 dt(1)];
an1 = LagrangePolynomialCoefficient(tn1, tt);
% quadratic
tn2 = [-dt(3) 0 dt(1)];
an2 = LagrangePolynomialCoefficient(tn2, tt);
% cubic
tn3 = [0-dt(2)-dt(3) 0-dt(3) 0 dt(1)];
an3 = LagrangePolynomialCoefficient(tn3, tt);
% P4
tn4 = [-sum(dt) -dt(2)-dt(3) -dt(3) 0 dt(1)];
an4 = LagrangePolynomialCoefficient(tn4, tt);
% P5
tn5 = [-sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0 dt(1)];
an5 = LagrangePolynomialCoefficient(tn5, tt);
% P6
tn6 = [-sum(dt)-dt(2)-dt(3) -sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0 dt(1)];
an6 = LagrangePolynomialCoefficient(tn6, tt);
An = [0 0 0 0 0 an1; 0 0 0 0 an2; 0 0 0 an3; 0 0 an4; 0 an5; an6];

% Stage 3
tt = t(3);
% linear
tn1 = [dt(1) dt(1)+dt(2)];
an1 = LagrangePolynomialCoefficient(tn1, tt);
% quadratic
tn2 = [0 dt(1) dt(1)+dt(2)];
an2 = LagrangePolynomialCoefficient(tn2, tt);
% cubic
tn3 = [0-dt(3) 0 dt(1) dt(1)+dt(2)];
an3 = LagrangePolynomialCoefficient(tn3, tt);
% P4
tn4 = [-dt(2)-dt(3) -dt(3) 0 dt(1) dt(1)+dt(2)];
an4 = LagrangePolynomialCoefficient(tn4, tt);
% P5
tn5 = [-sum(dt) -dt(2)-dt(3) -dt(3) 0 dt(1) dt(1)+dt(2)];
an5 = LagrangePolynomialCoefficient(tn5, tt);
% P6
tn6 = [-sum(dt)-dt(3) -sum(dt) -dt(2)-dt(3) -dt(3) 0 dt(1) dt(1)+dt(2)];
an6 = LagrangePolynomialCoefficient(tn6, tt);
An = [0 0 0 0 0 an1; 0 0 0 0 an2; 0 0 0 an3; 0 0 an4; 0 an5; an6];

