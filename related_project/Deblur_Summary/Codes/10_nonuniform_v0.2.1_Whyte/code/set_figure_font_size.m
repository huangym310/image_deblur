% function []=gg_set_figure_font_size(ffs,fsize)
%
% ffs ... figure handle
% fsize ... fontsize
%
function []=set_figure_font_size(ffs,fsize)


%    set(H(i,jj).hsht,'linewidth',2,'markersize',5)
%   set(ff(1))
%    ffs = ff(1);
ffa = findobj(ffs,'type','axes');
set(findobj(ffs,'type','axes'),'fontsize',fsize);
%delete(get(findobj(ffs,'type','axes'),'title')); % delete title
aa=get(findobj(ffs,'type','axes'),'xlabel');
if iscell(aa)
   for i = 1:length(aa)
       set(aa{i},'fontsize',fsize);
   end;
else
   set(aa,'fontsize',fsize);
end;

aa=get(findobj(ffs,'type','axes'),'ylabel');
if iscell(aa)
   for i = 1:length(aa)
       set(aa{i},'fontsize',fsize);
   end;
else
   set(aa,'fontsize',fsize);
end;

aa=get(findobj(ffs,'type','axes'),'zlabel');
if iscell(aa)
   for i = 1:length(aa)
       set(aa{i},'fontsize',fsize);
   end;
else
   set(aa,'fontsize',fsize);
end;
%   set(get(findobj(ffs,'type','axes'),'ylabel'),'fontsize',fsize);

aa = get(findobj(ffs,'type','axes'),'title');
if iscell(aa)
   for i = 1:length(aa)
       set(aa{i},'fontsize',fsize);
   end;
else
   set(aa,'fontsize',fsize);
end;

return;