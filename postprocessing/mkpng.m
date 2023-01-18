function mkpng(fn,frame,res)

if nargin<3, res=100; end

colormap(jet(65536));
delete(findobj(gcf,'tag','Colorbar'));
fn=sprintf('%s%05d.png',fn,frame);
print('-dpng',sprintf('-r%d',res),fn);


