function [n1,t1,t2] = plane(p1,p2,p3)

t1 = p2 - p1; 
t1 = t1/norm(t1);
t2 = p3 - p1; 
t2 = t2/norm(t2);
n1 = cross(t1,t2);
n1 = n1/norm(n1);
t2 = cross(n1,t1);

pm = (p1+p2+p3)/3;
if norm(p1+t2-pm)>norm(p1-t2-pm)
    t2 = -t2;
end

return; 

p1 = rand(1,3);
p2 = rand(1,3);
p3 = rand(1,3);
[n1,t1,t2] = plane(p1,p2,p3);

x = [p1(1) p2(1) p3(1) p1(1)];
y = [p1(2) p2(2) p3(2) p1(2)];
z = [p1(3) p2(3) p3(3) p1(3)];
figure(1); clf; plot3(x,y,z,'-o');
hold on;
x = [p1(1) p1(1)+n1(1)];
y = [p1(2) p1(2)+n1(2)];
z = [p1(3) p1(3)+n1(3)];
plot3(x,y,z,'-rx');
x = [p1(1) p1(1)+t1(1)];
y = [p1(2) p1(2)+t1(2)];
z = [p1(3) p1(3)+t1(3)];
plot3(x,y,z,'-rx');
x = [p1(1) p1(1)+t2(1)];
y = [p1(2) p1(2)+t2(2)];
z = [p1(3) p1(3)+t2(3)];
plot3(x,y,z,'-rx');
axis equal;

