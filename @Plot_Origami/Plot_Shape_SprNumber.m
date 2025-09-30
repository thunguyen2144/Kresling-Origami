%% Plot the configuration of the model

function Plot_Shape_SprNumber(obj)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;

% Plot Dots
figure
view(View1,View2); 
set(gca,'DataAspectRatio',[1 1 1])
set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])

% The software support two ways to set up the plotting range
A=size(Vsize);
if A(1)==1    
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
else
    axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
end

barNum=size(assembly.bar.A_vec);
barNum=barNum(1);
barConnect=assembly.bar.node_ij_mat;

for j=1:barNum
    node1=assembly.node.coordinates_mat(barConnect(j,1),:);
    node2=assembly.node.coordinates_mat(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k');
end

% Number Dots
node0=assembly.node.coordinates_mat;
A=size(assembly.node.coordinates_mat);
N=A(1);

sprIJKL=obj.assembly.rotSpr.node_ijkl_mat;
sprNum=size(sprIJKL);
sprNum=sprNum(1);


for i=1:sprNum
    x=0.5*(node0(sprIJKL(i,2),1)+...
        node0(sprIJKL(i,3),1));
    y=0.5*(node0(sprIJKL(i,2),2)+...
        node0(sprIJKL(i,3),2));
    z=0.5*(node0(sprIJKL(i,2),3)+...
        node0(sprIJKL(i,3),3));
    text(x,y,z,num2str(i),'Color','blue');
end

