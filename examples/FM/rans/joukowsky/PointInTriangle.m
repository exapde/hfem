function in = PointInTriangle(pt, v1, v2, v3)

b1 = sign(pt, v1, v2) < 0;
b2 = sign(pt, v2, v3) < 0;
b3 = sign(pt, v3, v1) < 0;

in = ((b1 == b2) && (b2 == b3));

function b = sign(p1, p2, p3)

 b = (p1(1) - p3(1)) * (p2(2) - p3(2)) - (p2(1) - p3(1)) * (p1(2) - p3(2));

