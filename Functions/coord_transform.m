function [x2,y2] = coord_transform(x1,y1,theta,rad)
if nargin < 4 %Default is deg
    rad = 0;
end

N = length(x1);
x2 = zeros(N,1);
y2 = zeros(N,1);
if rad == 0
    theta = theta*pi/180;
end
    
if length(theta) == 1;%just one angle to rotate
    theta = repmat(theta,N,1);
end

for i = 1:N
    x2(i) = cos(theta(i))*x1(i) + sin(theta(i))*y1(i);
    y2(i) = -sin(theta(i))*x1(i) + cos(theta(i))*y1(i);
end
return

