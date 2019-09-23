clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Ellipsoid RBF information files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID = 'ADP';
Weight_filename  = [ID,'_Weight.txt'];
DecayRate_filename = [ID,'_DecayRate.txt'];
Center_filename = [ID,'_Center.txt'];
Angle_filename = [ID,'_Angle.txt'];

cci = load(Weight_filename);
ddi = load(DecayRate_filename);
gp = load(Center_filename);
th = load(Angle_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained points initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pqrfile = [ID,'.pqr'];
fileID = fopen(pqrfile);
C = textscan(fileID,'%s %s %s %s %s %f %f %f %s %f');
fclose(fileID);
x0=[C{6},C{7},C{8},C{10}];
th1=0;
th2=0;
th3=1;

T1= [[cos(th1),-sin(th1),0];[sin(th1),cos(th1),0];[0,0,1]];
T2= [[cos(th2),0,-sin(th2)];[0,1,0];[sin(th2),0,cos(th2)]];
T3=[[1,0,0];[0,cos(th3),-sin(th3)];[0,sin(th3),cos(th3)]];         
                
TTT=T3'*T2'*T1';
x0(:,1:3)=   (T1*T2*T3*x0(:,1:3)')';
gsrange=[min(x0(:,1))-3,max(x0(:,1))+3,min(x0(:,2))-3,max(x0(:,2))+3,min(x0(:,3))-3,max(x0(:,3))+3];
di=0.5*ones(size(x0,1),1);
ci=exp(di.*x0(:,end).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained points initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th1 = th(:,1);
th2 = th(:,2);
th3 = th(:,3);
[xx,yy,zz]=meshgrid(gsrange(1):0.5:gsrange(2),gsrange(3):0.5:gsrange(4),gsrange(5):0.5:gsrange(6));
pp=0;
for k=1:size(gp,1)
    T1= [[cos(th1(k)),-sin(th1(k)),0];[sin(th1(k)),cos(th1(k)),0];[0,0,1]];
    T2= [[cos(th2(k)),0,-sin(th2(k))];[0,1,0];[sin(th2(k)),0,cos(th2(k))]];
    T3=[[1,0,0];[0,cos(th3(k)),-sin(th3(k))];[0,sin(th3(k)),cos(th3(k))]];
    XX= (T1*T2*T3*[xx(:)-gp(k,1),yy(:)-gp(k,2),zz(:)-gp(k,3)]')';
    t1=(exp(-(ddi(k,1)^2)*(XX(:,1)).^2-(ddi(k,2)^2)*(XX(:,2)).^2-(ddi(k,3)^2)*(XX(:,3)).^2));
    t=(cci(k)^2)*t1;
    pp=pp+t;
end
pp=reshape(pp,size(xx,1),size(xx,2),size(xx,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
zpp=0;
for k=1:size(x0,1)
    zpp=zpp+ci(k)*exp(-di(k)*((xx-x0(k,1)).^2+(yy-x0(k,2)).^2+(zz-x0(k,3)).^2));
end
[faces,verts]=isosurface(xx,yy,zz,zpp,1);
verts=(TTT*verts')';
p1 = patch('Faces',faces,'Vertices',verts);
p1.FaceColor = 'green';
p1.EdgeColor = 'none';
title('Original surface','interpreter','latex','FontSize',15);
daspect([1 1 1])
view(10,15);
axis tight
camlight
lighting gouraud
alpha(p1,0.5)
axis off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
pp=0;
for k=1:size(gp,1)
    T1= [[cos(th1(k)),-sin(th1(k)),0];[sin(th1(k)),cos(th1(k)),0];[0,0,1]];
    T2= [[cos(th2(k)),0,-sin(th2(k))];[0,1,0];[sin(th2(k)),0,cos(th2(k))]];
    T3=[[1,0,0];[0,cos(th3(k)),-sin(th3(k))];[0,sin(th3(k)),cos(th3(k))]];
    XX= (T1*T2*T3*[xx(:)-gp(k,1),yy(:)-gp(k,2),zz(:)-gp(k,3)]')';
    t1=(exp(-(ddi(k,1)^2)*(XX(:,1)).^2-(ddi(k,2)^2)*(XX(:,2)).^2-(ddi(k,3)^2)*(XX(:,3)).^2));
    t=(cci(k)^2)*t1;
    pp=pp+t;
end
pp=reshape(pp,size(xx,1),size(xx,2),size(xx,3));
[faces,verts]=isosurface(xx,yy,zz,pp,1);
verts=(TTT*verts')';     
p = patch('Faces',faces,'Vertices',verts);
p.FaceColor = 'red';
p.EdgeColor = 'none';
title('Final surface','interpreter','latex','FontSize',15);
daspect([1 1 1])
view(10,15);
axis tight
camlight
lighting gouraud
alpha(0.5)
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original surface overlapped with Final surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
zpp=0;
for k=1:size(x0,1)
    zpp=zpp+ci(k)*exp(-di(k)*((xx-x0(k,1)).^2+(yy-x0(k,2)).^2+(zz-x0(k,3)).^2));
end
[faces,verts]=isosurface(xx,yy,zz,zpp,1);
verts=(TTT*verts')';
p1 = patch('Faces',faces,'Vertices',verts);
p1.FaceColor = 'green';
p1.EdgeColor = 'none';
daspect([1 1 1])
view(10,15);
axis tight
camlight
lighting gouraud
alpha(p1,0.5)
axis off

hold on
pp=0;
for k=1:size(gp,1)
    T1= [[cos(th1(k)),-sin(th1(k)),0];[sin(th1(k)),cos(th1(k)),0];[0,0,1]];
    T2= [[cos(th2(k)),0,-sin(th2(k))];[0,1,0];[sin(th2(k)),0,cos(th2(k))]];
    T3=[[1,0,0];[0,cos(th3(k)),-sin(th3(k))];[0,sin(th3(k)),cos(th3(k))]];
    XX= (T1*T2*T3*[xx(:)-gp(k,1),yy(:)-gp(k,2),zz(:)-gp(k,3)]')';
    t1=(exp(-(ddi(k,1)^2)*(XX(:,1)).^2-(ddi(k,2)^2)*(XX(:,2)).^2-(ddi(k,3)^2)*(XX(:,3)).^2));
    t=(cci(k)^2)*t1;
    pp=pp+t;
end
pp=reshape(pp,size(xx,1),size(xx,2),size(xx,3));
[faces,verts]=isosurface(xx,yy,zz,pp,1);
verts=   (TTT*verts')';     
p = patch('Faces',faces,'Vertices',verts);
p.FaceColor = 'red';
p.EdgeColor = 'none';
title('Original surface overlapped with Final surface',...
    'interpreter','latex','FontSize',15);
daspect([1 1 1])
view(10,15);
axis tight
camlight
lighting gouraud
alpha(0.5)
axis off
