function h = showDescriptor(X,T,desc)

h = figure;
patch('vertices',X,'Faces',T,'CData',desc,'FaceColor','interp','edgecolor','none');
axis equal;
axis off;
cameratoolbar;

end
