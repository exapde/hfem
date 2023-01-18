function [Rp Rp_q] = press(q)
% Pressure equation \epsilon p - Rp = 0

[ng ncq] = size(q);
if ncq==3
    Rp = 0.5*(q(:,1) + q(:,3));

    Rp_q = zeros(ng,ncq);
    Rp_q(:,1) = 0.5;
    Rp_q(:,ncq) = 0.5;
else
    Rp = 0.5*(q(:,1) + q(:,4) + q(:,6));

    Rp_q = zeros(ng,ncq);
    Rp_q(:,1) = 0.5;
    Rp_q(:,4) = 0.5;
    Rp_q(:,ncq) = 0.5;    
end
