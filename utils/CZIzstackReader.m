function stack = CZIzstackReader(fileName)
% function to read multichannel czi data

h = MEXlibCZI('Open', fileName);
info = MEXlibCZI('GetInfo', h);
Z = info.dimBounds.Z(2);

temp = MEXlibCZI('GetSubBlockBitmap', h, 0)';
[M,N] = size(temp);
stack = zeros(M,N,Z,'single');
for loop = 0:Z-1
    % subblock_handle = MEXlibCZI('GetSubBlock', h, loop);
    % infoBlock = MEXlibCZI('GetInfoFromSubBlock', h, subblock_handle);
    stack(:,:,loop+1) = MEXlibCZI('GetSubBlockBitmap', h, loop)';
    disp(['read z-plane #:', num2str(loop)])
end
MEXlibCZI('Close',h)

end

% function [C] = splitChannelID(x)
% 
% if contains(x, 'Z') & contains(x, 'C') & contains(x, 'T') & contains(x, 'H')
%     s = split(x,'C');
%     s = split(s(2),'T');
%     C = str2double(cell2mat( s(1) ));
% else
%     C = split(x,'C');
%     C = str2double(C(2));
% end
% end


