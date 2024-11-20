% clear; clc;

clear all;

datadir = '../data';
procs = 0:15;
n = 1;
klevel = 10;
particle_type=2;

%%
% cellcentered data contains the voronoi points and the depths
% at those points.
N=length(procs);

tri = [];
x = [];
y = [];
T = [];

for proc=procs
    %
%     filename = [datadir,'/T.dat.'];
%     fid = fopen(filename,'rb');
    pd = load([datadir,'/points.dat']);
    cp = load([datadir,'/cells.dat.',num2str(proc)]);
    cdp = load([datadir,'/celldata.dat.',num2str(proc)]);

    xp = pd(:,1);
    yp = pd(:,2);
    tri_proc = cp(:,3:5);
    d = cdp(:,4);
    xv = cdp(:,1);
    yv = cdp(:,2);
    dz = load([datadir,'/vertspace.dat']);

    % Total number of cells in the horizontal Nc and vertical Nk
    Nc = length(xv);
    Nk = length(dz);
    z = -(sum(dz(1:klevel))-dz(klevel)/2)*ones(size(xv));

%     arraysize = Nc*Nk;
%     fseek(fid,8*(arraysize*(n-1)+(klevel-1)*Nc),0);
%
%     phi = fread(fid,Nc,'float64');
%
%     phi(find(z<-d))=nan;
%
    tri = vertcat(tri, tri_proc);
    x = vertcat(x,xp);
    y = vertcat(y,yp);
%     T = vertcat(T,phi);

end

% unsurf(tri,x,y,0,'FaceColor','none','EdgeColor',[0.2,0.2,0.2],'LineWidth',0.4);

%%
lag=[];


lag_start=1;
lag_interval=10;
lag_end=25000;


for proc=procs
    fname=[datadir,'/lagout.dat.',num2str(proc)];
    %lagout=load("-ascii", fname);
    fid = fopen(fname, 'r');
    [lagout, count] = fscanf(fid, '%d %d %d %f %f %f');
    lagout = reshape(lagout, 6, size(lagout)/6)';
    fclose(fid);
    if (~isempty(lagout))
        lag=vertcat(lagout,lag);
    end
end

par=load('../data/2d_particles.dat');

xp=zeros(length(par),length(lag_start:lag_interval:lag_end)+1);
yp=zeros(length(par),length(lag_start:lag_interval:lag_end)+1);
zp=zeros(length(par),length(lag_start:lag_interval:lag_end)+1);

for n=1:length(lag)
    tp=floor((lag(n,1)-lag_start)/lag_interval)+2;
    pid=lag(n,2)+1;
    xp(pid,tp)=lag(n,4);
    yp(pid,tp)=lag(n,5);
    zp(pid,tp)=lag(n,6);
end

% for n=lag_start:lag_interval:lag_end
%     xp(:,p)=x(lag(:,1)==n);
%     yp(:,p)=y(lag(:,1)==n);
%     if particle_type==3
%         wp(:,p)=z(lag(:,1)==n);
%     end
%     p=p+1;
% end

%%

for p=1:length(par)
%p = floor(length(par)/2)
    scatter(xp(p,:),zp(p,:),'.');
    hold on;
end

grid on;
xlabel('x','fontname','Arial','fontsize',18);
ylabel('y','fontname','Arial','fontsize',18);
set(gca, 'YDir','reverse')
print(gcf,'-dtiff','-r300','parts');
% close;
