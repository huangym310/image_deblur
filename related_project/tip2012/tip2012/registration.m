function [R, H] = registration(G,blocks)

% Registration
% Images is divided into nonoverlapping windows and in each a shift by normalized
% correlation is estimated. Affine transform paramaters are calculated from
% the local shifts.
%
% [R,H] = registration(G,blocks)
%
% G ... input images (cell array)
% blocks ... number of blocks [m,n]
% 
% R ... registered G images
% H ... calcualted homographies
%
% Reference image is determined by refind variable

refind = 1;% identify reference image
verbose = 2;
isize = size(G{1});
gsize = isize(1:2);
P = length(G);
if size(G{1},3) > 1
    cind = 2; %green channel
else
    cind = 1;
end
    

% number of blocks and block size is given
if iscell(blocks)
   noblocks = blocks{1};
   if length(blocks) < 2
      bsize = ceil(gsize./noblocks);
   else
      bsize = blocks{2};
   end
   if sum(bsize > gsize) 
      error('Block size is larger than the image size!');
   end
% only number of blocks is given; block size is calculated with no overlaps
else
   noblocks = blocks;
   bsize = ceil(gsize./noblocks);
end
disp(['Block size: ',num2str(bsize)]);

c = ceil(bsize/2);
sc = ceil(ceil(gsize./noblocks)/2);
scl = max(c,sc);
scr = max(bsize-c,sc-1);
posy = linspace(scl(1),gsize(1)-scr(1),noblocks(1));
posx = linspace(scl(2),gsize(2)-scr(2),noblocks(2));
locx = round(posx); locy = round(posy);
[X Y] = meshgrid(locx,locy);
loc = [Y(:).'; X(:).'];
bb = [loc(1,:)-c(1)+1; loc(1,:)+(bsize(1)-c(1));...
   loc(2,:)-c(2)+1; loc(2,:)+(bsize(2)-c(2)) ]; 
roi = [1+floor(bsize/4), bsize-floor(bsize/4)];
roi = [roi(1) roi(3) roi(2) roi(4)];
j = 1;
N = size(loc,2);
tv = zeros(2,N,P);
%sigma = zeros(N,1);
mask = zeros(size(tv));
%  estimate local shifts
for i = bb
   S = getBlock(G,i,cind); 
   tv(:,j,:) = preReg(S,roi,refind).';
   %sigma(j) = var(vec(S(:,:,refind)));
   j = j + 1;
end
% eliminate wrong shifts
mask = (abs(tv(1,:,:)) <= (bsize(1)-(roi(2)-roi(1)+1))/2) & ...
    (abs(tv(2,:,:)) <= (bsize(2)-(roi(4)-roi(3)+1))/2);
mask = squeeze(mask);
%mask = true(N,P);
% if verbose > 1
%    figure; dispIm(linscale(G{2}));
%    for i = bb
%       h = rectangle('Position',[i(3) i(1) i(4)-i(3) i(2)-i(1)],...
%          'EdgeColor',[1 0 0]);
%    end
%    line([loc(2,mask(:,2)); loc(2,mask(:,2))+tv(2,mask(:,2),2)],...
%        [loc(1,mask(:,2)); loc(1,mask(:,2))+tv(1,mask(:,2),2)],'Color',[1 0 0]);
%    %unvec(sigma,blocks)
% end 
% move the origin to the image center and put X before Y so it is compatible with my
% homography functions 
loc = loc-repmat(ceil(gsize/2).',[1 N]);
loc = loc([2 1],:);
tv = tv([2 1],:,:);
% calculate parameters of affine transform
a = zeros(6,P);
for j = 1:P
    if j == refind
        H{1} = eye(3);
        continue;
    end
    Z = zeros(1,nnz(mask(:,j)));
    O = ones(1,nnz(mask(:,j)));
    XT = [ loc(1,mask(:,j)), Z; ...
    Z, loc(1,mask(:,j)); ...
    loc(2,mask(:,j)), Z; ...
    Z, loc(2,mask(:,j)); ...
    O, Z;...
    Z, O];
    a(:,j) = (XT*XT.')\(XT*vec((loc(:,mask(:,j))+tv(:,mask(:,j),j)).'));
    H{j} = [reshape(a(1:4,j),[2 2]),a(5:6,j); [0 0 1]]; 
end   

% perform interpolation
%bb = getMaxInRec(H,gsize);
%[X Y] = meshgrid(bb(3):bb(4),bb(1):bb(2)); 
[X Y] = meshgrid([1:isize(2)]-ceil(isize(2)/2), ...
     [1:isize(1)]-ceil(isize(1)/2));
for j = 1:P
    R{j} = zeros(size(G{1}));
    for c = 1:size(G{1},3)
        R{j}(:,:,c) = homography(H{j},X,Y,G{j}(:,:,c));
    end
end
%crop images to get valid parts
nanRec = ~isnan(R{1}(:,:,cind));
for j = 2:P
    nanRec = nanRec.*~isnan(R{j}(:,:,cind));
end
i=0;
while ~isempty(find(nanRec(1+i:end-i,1+i:end-i) == 0, 1))
    i = i+1;
end
i = i+1;
for j=1:P
    R{j} = R{j}(1+i:end-i,1+i:end-i,:); 
end
end

function tvec = preReg(I,roi,refind)
%
% registration using normalized correlation
% local shift estimation
%

sa = size(I);
tvec = zeros(sa(3),2);
B = I(roi(1):roi(2),roi(3):roi(4),refind); 
sb = size(B);

CFB = conj(fft2(B,sa(1),sa(2)));
CFO = conj(fft2(ones(sb),sa(1),sa(2)));
C = real(ifft2(fft2(I).*repmat(CFB,[1 1 sa(3)])))./...
    sqrt(real(ifft2(fft2(I.^2).*repmat(CFO,[1 1 sa(3)]))));
for ii = 1:sa(3)
    [k l] = find(C(:,:,ii) == max(vec(C(:,:,ii))));
    tvec(ii,:) = [k(1) l(1)];
end
tvec = tvec - repmat([roi(1),roi(3)],sa(3),1);
end

function S=getBlock(G,bb,ind)
S = zeros(bb(2)-bb(1)+1,bb(4)-bb(3)+1,length(G));
   for i = 1:length(G)
      S(:,:,i) = G{i}(bb(1):bb(2),bb(3):bb(4),ind);
   end
end 