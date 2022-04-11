
%{ 

% out : Linear Contrast Stretched Image
% x : Input nodes 

y=[5 70 125 200 255];   % output nodes
out=linearCS(img,y) ;    % img ; monochromatic (single-band) image
imshow([ img out ]), shg



(c) P. Civicioglu, E. Besdok, 23.March.2022
%}


function [J,x]=linearCS(img,y)
[M,N]=size(img);
img=double(img);
Imin=min(img(:));
Imax=max(img(:));
x=linspace(Imin,Imax,numel(y));
fnc=fit(x',y','linear');
J=reshape(fnc(img),M,N);







