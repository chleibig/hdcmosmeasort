function showrois( CC, varargin)
%showrois( CC ) constructs a labelled matrix from CC and plots
%it in HSV colours. This function allows to visualize connected components
%or regions of interest without using labelmatrix and label2rgb which
%require the image processing toolbox; 1 indicates the background
%
%showrois( CC , whichones) shows only the ones given in the vector
%                                whichones
%
%Note that that overlapping regions may be hidden by others
%
% Input
% =====
%
% CC - struct as returned from bwconncomp with at least the fields
%    *.ImageSize
%    *.NumObjects
%    *.PixelIdxList
%
%
% christian.leibig@g-node.org, 05.09.2013
%


if nargin == 1
    %Show all.
    whichones = 1:CC.NumObjects;
end

if nargin == 2
    %Show only the given selection.
    whichones = varargin{1};
end

L = zeros(CC.ImageSize);
for i = 1:length(whichones)
    L(ind2sub(CC.ImageSize,CC.PixelIdxList{whichones(i)})) = whichones(i) + 1;
end
image(L);colormap([1 1 1;hsv(CC.NumObjects)]);
%title('ROIs');xlabel('sensor columns');ylabel('sensor rows');
axis square;
%colorbar;

end

