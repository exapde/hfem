


param = {1.4,0,1e6,0.72,0.5,1.5};

UDG = rand(1,12);
V   = rand(2,1);
G   = [1+rand rand; rand 1+rand];
g   = det(G);
gb  = g;
dgbdx = rand(2,1);
p   = [V; G(:); g; gb; dgbdx];
p   = [rand(4,1); p];
p   = p';

[f,f_UDG] = flux(p,UDG,param,[]);

sz = size(f_UDG);n  = sz(end);m  = prod(sz)/n;
f_UDG = reshape(f_UDG,[m n]);

es = 1e-6;
for i = 1:length(UDG)
    UDG1 = UDG;
    UDG1(i) = UDG(i) + es;
    f1 = flux(p,UDG1,param,[]);

    UDG2 = UDG;
    UDG2(i) = UDG(i) - es;
    f2 = flux(p,UDG2,param,[]);
    
    df = (f1(:)-f2(:))/(2*es);
    e(i)  = max(abs(df-f_UDG(:,i)));        
end
e


nl = rand(1,2);
nl = nl/sqrt(nl(1)^2+nl(2)^2);
UH = rand(1,4);
[fh,fh_UDG,fh_UH] = fhat(nl,p,UDG,UH,param,[],[]);

sz = size(fh_UDG);n  = sz(end);m  = prod(sz)/n;
fh_UDG = reshape(fh_UDG,[m n]);
sz = size(fh_UH);n  = sz(end);m  = prod(sz)/n;
fh_UH = reshape(fh_UH,[m n]);

for i = 1:length(UDG)
    UDG1 = UDG;
    UDG1(i) = UDG(i) + es;
    fh1 = fhat(nl,p,UDG1,UH,param,[],[]);

    UDG2 = UDG;
    UDG2(i) = UDG(i) - es;
    fh2 = fhat(nl,p,UDG2,UH,param,[],[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    eh(i)  = max(abs(dfh-fh_UDG(:,i)));        
end
eh

for i = 1:length(UH)
    UH1 = UH;
    UH1(i) = UH(i) + es;
    fh1 = fhat(nl,p,UDG,UH1,param,[],[]);

    UH2 = UH;
    UH2(i) = UH(i) - es;
    fh2 = fhat(nl,p,UDG,UH2,param,[],[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    euh(i)  = max(abs(dfh-fh_UH(:,i)));        
end
euh



nl = rand(1,2);
nl = nl/sqrt(nl(1)^2+nl(2)^2);
UH = rand(1,4);
[fh,fh_UDG,fh_UH] = fhatwall(nl,p,UDG,UH,param,[],[]);

sz = size(fh_UDG);n  = sz(end);m  = prod(sz)/n;
fh_UDG = reshape(fh_UDG,[m n]);
sz = size(fh_UH);n  = sz(end);m  = prod(sz)/n;
fh_UH = reshape(fh_UH,[m n]);

for i = 1:length(UDG)
    UDG1 = UDG;
    UDG1(i) = UDG(i) + es;
    fh1 = fhatwall(nl,p,UDG1,UH,param,[],[]);

    UDG2 = UDG;
    UDG2(i) = UDG(i) - es;
    fh2 = fhatwall(nl,p,UDG2,UH,param,[],[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    ehw(i)  = max(abs(dfh-fh_UDG(:,i)));        
end
ehw

for i = 1:length(UH)
    UH1 = UH;
    UH1(i) = UH(i) + es;
    fh1 = fhatwall(nl,p,UDG,UH1,param,[],[]);

    UH2 = UH;
    UH2(i) = UH(i) - es;
    fh2 = fhatwall(nl,p,UDG,UH2,param,[],[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    euhw(i)  = max(abs(dfh-fh_UH(:,i)));        
end
euhw


ib=1;
nl = rand(1,2);
nl = nl/sqrt(nl(1)^2+nl(2)^2);
UH = rand(1,4);
ui = rand(1,4);
[fh,fh_UDG,fh_UH] = fbou(ib,ui,nl,p,UDG,UH,param,[]);

sz = size(fh_UDG);n  = sz(end);m  = prod(sz)/n;
fh_UDG = reshape(fh_UDG,[m n]);
sz = size(fh_UH);n  = sz(end);m  = prod(sz)/n;
fh_UH = reshape(fh_UH,[m n]);

for i = 1:length(UDG)
    UDG1 = UDG;
    UDG1(i) = UDG(i) + es;
    fh1 = fbou(ib,ui,nl,p,UDG1,UH,param,[]);

    UDG2 = UDG;
    UDG2(i) = UDG(i) - es;
    fh2 = fbou(ib,ui,nl,p,UDG2,UH,param,[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    eb(i)  = max(abs(dfh-fh_UDG(:,i)));        
end
eb

for i = 1:length(UH)
    UH1 = UH;
    UH1(i) = UH(i) + es;
    fh1 = fbou(ib,ui,nl,p,UDG,UH1,param,[]);

    UH2 = UH;
    UH2(i) = UH(i) - es;
    fh2 = fbou(ib,ui,nl,p,UDG,UH2,param,[]);
    
    dfh = (fh1(:)-fh2(:))/(2*es);
    ebh(i)  = max(abs(dfh-fh_UH(:,i)));        
end
ebh


absolute = 0;
[An,Anm] = getan(p,nl,UH,param,absolute);
sz = size(Anm);n  = sz(end);m  = prod(sz)/n;
Anm = reshape(Anm,[m n]);

for i = 1:length(UH)
    UH1 = UH;
    UH1(i) = UH(i) + es;
    An1 = getan(p,nl,UH1,param,absolute);

    UH2 = UH;
    UH2(i) = UH(i) - es;
    An2 = getan(p,nl,UH2,param,absolute);
    
    dfh = (An1(:)-An2(:))/(2*es);
    em(i)  = max(abs(dfh-Anm(:,i)));        
end
em

