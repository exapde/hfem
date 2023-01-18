
function setapplicationpath(application)

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'digaso') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'DIGASO'))
    up = up+1;
end
ss=strcat(ps(1:is(end-up)),'application');
warning off; rmpath(genpath(ss)); warning on;
addpath(genpath(strcat(ss,sslash,application)));
