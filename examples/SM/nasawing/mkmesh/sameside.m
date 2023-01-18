function same = sameside(p1, p2, q1, q2)

x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);
a1 = q1(1); b1 = q1(2);
a2 = q2(1); b2 = q2(2);

if ((y1-y2)*(a1-x1)+(x2-x1)*(b1-y1))*((y1-y2)*(a2-x1)+(x2-x1)*(b2-y1))>0
    same = 1;
else
    same = 0;
end
