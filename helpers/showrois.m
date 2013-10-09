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
%
% christian.leibig@g-node.org, 05.09.2013
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Show all regions of interest.
selection = 1:CC.NumObjects;

fillColumnFlag = false;

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
        otherwise
            % Hmmm, something wrong with the parameter string
            error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = zeros(CC.ImageSize);
for i = 1:length(selection)
    if fillColumnFlag
        [r,c] = ind2sub(CC.ImageSize,CC.PixelIdxList{selection(i)});
        %make use of MATLAB internal vector representation to fill skipped
        %columns:
        r = [r(:)';r(:)'];
        r = r(:)';
        c = [c(:)';c(:)'+1];
        c = c(:)';        
        L(sub2ind(CC.ImageSize,r,c)) = selection(i) + 1;
    else
        L(CC.PixelIdxList{selection(i)}) = selection(i) + 1;
    end
end
image(L);colormap([1 1 1;hsv(CC.NumObjects)]);
%title('ROIs');xlabel('sensor columns');ylabel('sensor rows');
axis square;
%set axes limits to visible area:
[r, c] = find(L);
ylim([min(r),max(r)]);
xlim([min(c),max(c)]);
%colorbar;

end

