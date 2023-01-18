function ent2cpu = edgpartition(elcon, f, t, t2f, re, ce, elem2cpu, face2cpu, nproc, overlappinglevel)

[nfe,ne] = size(t2f);
npf = numel(elcon)/(nfe*ne);
nf = size(f,2);

elcon = reshape(elcon,[npf,nfe ne]);
facecon = -ones(npf,nf);
for i = 1:nf
    fe = f(end-1:end,i)+1; % neighboring elements of face i
    if1 = t2f(:,fe(1))==i-1; % location of face i on element fe(1)
    facecon(:,i) = elcon(:,if1,fe(1));
end

entmax = max(facecon(:))+1;
ent2cpu = -ones(entmax,1);
nent = 0;
for i = 1:nproc
    disp(['processor: ' num2str(i)]);

    % list of  elements in nonoverlapping subdomain i
    intelem = find(elem2cpu==(i-1));

    % list of  faces in nonoverlapping subdomain i
    intface = find(face2cpu==(i-1));

    elem = intelem;
    for j=1:overlappinglevel
        % TO DO: C++ version of node2elem
        elem = node2elem(t(:,elem)+1,re,ce);
    end
    extelem = setdiff(elem,intelem);
    extintelem = [intelem; extelem];

    face = t2f(:,extintelem)+1;
    face = unique(face(:));
    extface = setdiff(face,intface);
    extintface = [intface; extface];

    intent = mkintent(facecon(:,extintface)+1,length(intface),face2cpu(extintface),i-1);
    ent2cpu(intent) = i-1;

    nent = nent + length(intent);
end

if nent~=entmax
    error('EDG partition is incorrect');
end
if min(ent2cpu(:))<0
    error('EDG partition is incorrect');
end
