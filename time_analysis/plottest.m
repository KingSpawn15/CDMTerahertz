a=-8:.5:8;
[X,Y] = meshgrid(a); 

R = sqrt(X.^2 + Y.^2) + eps;

Z = sin(R)./R;
hold off
surf(X,Y,Z)
hold on
colormap hsv
alpha(.4)

zMaxY=max(Z);
zMaxX=max(Z,[],2)';
zMinY=min(Z);
zMinX=min(Z,[],2)';

surf([a;a]',8*ones(size([a;a]))',[zMaxY;zMinY]')
alpha(.4)
surf([a;a]',-8*ones(size([a;a]))',[zMaxY;zMinY]')
alpha(.4)
surf(8*ones(size([a;a]))',[a;a]',[zMaxX;zMinX]')
alpha(.4)
surf(-8*ones(size([a;a]))',[a;a]',[zMaxX;zMinX]')
alpha(.4)