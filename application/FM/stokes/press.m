function [Rp Rp_q] = press(q)
% Pressure equation \epsilon p - Rp = 0

[ng ncq] = size(q);
nd = round(sqrt(ncq));

q  = reshape(q,[ng,nd,nd]);

Rp = q(:,1,1);
Rp_q = zeros(ng,nd,nd);
Rp_q(:,1,1) = 1;
for d=2:nd
    Rp = Rp + q(:,d,d);
    Rp_q(:,d,d) = 1;
end
Rp_q = reshape(Rp_q,[ng ncq]);

