function showrois( CC, varargin)
%showrois( CC ) constructs a labelled matrix from CC and plots
%it in HSV colours. This function allows to visualize connected components
%or regions of interest without using labelmatrix and label2rgb which
%require the image processing toolbox; 1 indicates the background
%
%Note that overlapping regions may be hidden by others
%
% Input
% =====
%
% CC - struct as returned from bwconncomp with at least the fields
%    *.ImageSize
%    *.NumObjects
%    *.PixelIdxList
%
% Optional arguments are to be given as key/value pairs:
%
% 'selection',selection - vector specifying which regions of interest to
%                         show
% 'fillColumnFlag',fillColumnFlag - if true columns skipped (in recordings)
%                                   get filled for visualization purposes
%
% 'shuffle' - if true, colours get shuffled (default: false)
%
%
% christian.leibig@g-node.org, 05.09.2013
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Show all regions of interest.
selection = 1:CC.NumObjects;

fillColumnFlag = false;

shuffle = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch varargin{i}
        case 'selection'
            selection = (varargin{i+1});
        case 'fillColumnFlag'
            fillColumnFlag = (varargin{i+1});
        case 'shuffle'
            shuffle = (varargin{i+1});
        otherwise
            % Hmmm, something wrong with the parameter string
            error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize images
background = ones(CC.ImageSize);%white background
Lall = background;%overlay all in one without transparency
%Stack all ROI images on top of each other into the array layers 
%using 3 RGB channels for each ROI image; 
layers = NaN(CC.ImageSize(1),CC.ImageSize(2),3);

%construct (alpha-blended) colormap.
ALPHA = 0.5;
M = hsv(CC.NumObjects);
M = rgb2hsv(M);
M(:,2) = ALPHA * M(:,2);%just rescale saturation
M = hsv2rgb(M);
if shuffle
    M = M(randperm(size(M,1)),:);
end
M = [NaN NaN NaN;M];%[NaN NaN NaN] will be replaced by [1 1 1] 
%later on and is for the white background, all other colours are
%alpha-blended

for i = 1:length(selection)
    if fillColumnFlag
        [r,c] = ind2sub(CC.ImageSize,CC.PixelIdxList{selection(i)});
        %make use of MATLAB internal vector representation to fill skipped
        %columns:
        r = [r(:)';r(:)'];
        r = r(:)';
        c = [c(:)';c(:)'+1];
        c = c(:)';
        %overlapping ROI information is overwritten in Lall
        Lall(sub2ind(CC.ImageSize,r,c)) = selection(i) + 1;
        %overlapping ROI information is kept in Li
        Li = background;
        Li(sub2ind(CC.ImageSize,r,c)) = selection(i) + 1;
        layers = cat(4,layers,ind2rgb(Li,M));
        clear Li
    else
        %overlapping ROI information is overwritten in Lall
        Lall(CC.PixelIdxList{selection(i)}) = selection(i) + 1;
        %overlapping ROI information is kept in Li
        Li = background;
        Li(CC.PixelIdxList{selection(i)}) = selection(i) + 1;
        layers = cat(4,layers,ind2rgb(Li,M));
        clear Li;
    end
end

%without alpha blending:
%image(Lall);colormap(M);

%composite alpha blended image:
transparentROIs = nanmean(layers,4);%nanmean guarantees that we only take
%into account non-background pixels
transparentROIs(isnan(transparentROIs)) = 1;%now we assign the background
%to become white
image(transparentROIs);

% %with alpha blending via matlab (expensive and impractical for GUI):
% 
% %the following code would be faster if one would calculate the composite
% %image and then call image() only once.
% a = 0.5;%use the same alpha value for all layers
% for l=1:size(layers,4);
%     hold on;
%     image(layers(:,:,1:3,l),'AlphaData',a*layers(:,:,4,l));
% end

%title('ROIs');xlabel('sensor columns');ylabel('sensor rows');
axis square;
%set axes limits to visible area:
[r, c] = find(Lall);
ylim([min(r),max(r)]);
xlim([min(c),max(c)]);
%colorbar;

end



