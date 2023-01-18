function mkmovie(dir,sz,fn)

if nargin==1
    cd(dir);
end

if nargin<2
    !mogrify -crop 0x0 *.png
    !mogrify -geometry 100%x100% *.png
else
    str = ['!mogrify -crop ' num2str(sz(1)) 'x' num2str(sz(2)) '+' num2str(sz(3)) '+' num2str(sz(4)) ' *.png'];
    eval(str);
end

if nargin < 3
    !mencoder "mf://*.png" -mf fps=24 -o a.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq:vbitrate=1000:vpass=1
    !mencoder "mf://*.png" -mf fps=24 -o a.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq:vbitrate=2000:vpass=2
else
    str = ['!mencoder "mf://*.png" -mf fps=24 -o ' fn '.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq:vbitrate=2000:vpass=1'];
    eval(str);
    str = ['!mencoder "mf://*.png" -mf fps=24 -o ' fn '.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq:vbitrate=2000:vpass=2'];
    eval(str);
end

!mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mjpeg -oac copy -o a.avi
!mencoder mf://*.png -mf fps=24 -ovc lavc -lavcopts vcodec=mjpeg -o a.avi

