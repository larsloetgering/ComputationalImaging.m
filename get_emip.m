function r = get_emip(zstack, percentage)
% computes extended maximum intensity projection
%
% r = get_emip(zstack, percentage)
% zstack: array with dimension (x,y,z)
% percentage: determines how many pixels along z are used
% zstack = image intensity (x,y,z)

x = sort(zstack, 3,'descend'); 
r = sum(x(:,:,1:round(percentage*size(x,3)) ), 3);
end