function [im, cmap] = getColorCodedDepthProjection(zstack, varargin)
% DOCUMENTATION

% INPUT PARSER
p = inputParser;
p.addParameter('percentage', 0.5)
p.addParameter('cmap', 'rainbow')
p.parse( varargin{:} );

% paraemters needed for loop
[M,N,numz] = size(zstack);

% evaluate cmap
switch p.Results.cmap
    case 'rainbow'
        cmap = rainbowCMAP(numz);
    otherwise
        cmap = colormap(p.Results.cmap);
end

% sort stack
stack_sorted = sort(zstack, 3,'descend'); 
k = round(numz*p.Results.percentage);
weights = zstack > stack_sorted(:,:,k);
im = 0;
for k = 1:numz
    % im = im + gray2rgb(weights(:,:,k).*zstack(:,:,k),cmap(k,:));
    im = im + reshape( reshape(weights(:,:,k) .* zstack(:,:,k),[M*N,1]) * cmap(k,:), [M,N,3]);
end
im = im ./ max(im(:));

end