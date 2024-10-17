clear all; clc; close all;

%--- change from here

pbcmap1.a =   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 6;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 5;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 4;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
          
pbcmap1.b =   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            6 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            5 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            4 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            3 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            2 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

pbcmap2.a =   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 20 21;
            0 0 0 0 0 0 0 0 0 0 0 0 18 19 0;
            0 0 0 0 0 0 0 0 0 0 0 16 17 0 0;
            0 0 0 7 8 9 10 11 12 13 14 15 0 0 0;
            0 0 5 6 0 0 0 0 0 0 0 0 0 0 0;
            0 3 4 0 0 0 0 0 0 0 0 0 0 0 0;
            1 2 0 0 0 0 0 0 0 0 0 0 0 0 0];
          
pbcmap2.b =   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 20 21;
            0 0 0 0 0 0 0 0 0 0 0 0 18 19 0;
            0 0 0 0 0 0 0 0 0 0 0 16 17 0 0;
            0 0 0 7 8 9 10 11 12 13 14 15 0 0 0;
            0 0 5 6 0 0 0 0 0 0 0 0 0 0 0;
            0 3 4 0 0 0 0 0 0 0 0 0 0 0 0;
            1 2 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


exist =   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
            0 0 0 0 0 0 0 0 0 0 0 0 1 1 1;
            0 0 0 0 0 0 0 0 0 0 0 1 1 0 1;
            0 0 0 1 1 1 1 1 1 1 1 1 0 0 1;
            0 0 1 1 1 0 0 0 0 0 0 0 0 0 1;
            0 1 1 0 1 1 1 1 0 0 0 0 0 1 1;
            1 1 0 0 0 0 0 1 1 1 1 0 1 1 0;
            1 0 0 0 0 0 0 0 0 0 1 1 1 0 0;
            1 0 0 1 1 1 1 1 1 1 1 1 0 0 0;
            1 0 1 1 0 0 0 0 0 0 0 0 0 0 0;
            1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
            1 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
            
pbcmap1.elts = 6;
pbcmap2.elts = 21;

%---- to here

pbcmap.pbcmap1 = pbcmap1;
pbcmap.pbcmap2 = pbcmap2;
% pbcmap.pbcmap3 = pbcmap3;
        
[X,T]=genMesh(exist);
X = X-mean(X);
mean(X)
PBC = genPBC(pbcmap,X,T,exist);

dlmwrite("./X.csv",X,'delimiter',',','precision',5)
dlmwrite("./T.csv",T-1,'delimiter',',','precision',5)
dlmwrite("./PBC.csv",PBC,'delimiter',',','precision',5)


function T = killVoid(exist,T)
  exist90 = rot90(exist,-1);
  Tp=[];
  for i=1:(size(exist,1)*size(exist,2))
    if exist90(i)==1
      Tp = [Tp; T(i,:)];
    end
  end
  T=Tp;
end

function [X,T] = createHEXmesh(dom,nx,ny)
  npx = nx+1; npy = ny+1; npt = npx * npy; x1 = dom(1); x2 = dom(2);
  y1 = dom(3); y2 = dom(4); x = linspace(x1,x2,npx);  y = linspace(y1,y2,npy); 
  [x,y] = meshgrid(x,y);  X = [reshape(x',npt,1), reshape(y',npt,1)]; 
  T = zeros(nx*ny,4);
  for i=1:ny
    for j=1:nx
      ielem = (i-1)*nx+j; inode = (i-1)*(npx)+j;
      T(ielem,:) = [inode   inode+1   inode+npx+1   inode+npx];
    end
  end
end

function genPLY(X,T)

  tri = [1 2 3;
         1 3 4;
         5 6 7;
         5 7 8;
         1 2 6;
         1 6 5;
         4 3 7;
         4 7 8;
         2 3 7;
         2 7 6;
         1 4 8;
         1 8 5];
  Ttri = [];
  for i=1:size(T,1)
    for j=1:size(tri,1)
      Ttri = [Ttri; T(i,tri(j,:))];
    end
  end

  fileID = fopen('essai.ply','w');
  fprintf(fileID,"ply\n")
  fprintf(fileID,"format ascii 1.0\n")
  fprintf(fileID,"comment VCGLIB generated\n")
  fprintf(fileID,"element vertex %i\n",size(X,1))
  fprintf(fileID,"property float x\n")
  fprintf(fileID,"property float y\n")
  fprintf(fileID,"property float z\n")
  fprintf(fileID,"element face %i\n",size(Ttri,1))
  fprintf(fileID,"property list uchar int vertex_indices\n")
  fprintf(fileID,"end_header\n")

  strX = '%5f %5f %5f\n';
  for i=1:size(X,1)
    fprintf(fileID,strX,X(i,:))
  end

  strT = "3 %i %i %i\n";
  Ttri = Ttri-1;
  for i=1:size(Ttri,1)
    fprintf(fileID,strT,Ttri(i,:))
  end


  fprintf(fileID,"//+\n")
  fprintf(fileID,"Recombine Surface {1};\n")
  fprintf(fileID,"//+\n")
  fprintf(fileID,"Recombine Surface {1};\n")

end

function genMESH(X,T)


  strT = '%i %i %i %i %i %i %i %i %i\n';
  strX = '%5f %5f %5f %i\n';
  fileID = fopen('essai2.mesh','w');
  fprintf(fileID,"MeshVersionFormatted 1\n");
  fprintf(fileID,"Dimension 3\n");
  fprintf(fileID,"Vertices %i\n", size(X,1));
  for i=1:size(X,1)
    fprintf(fileID, strX, [X(i,:) 0]);
  end
  fprintf(fileID,"Hexahedra %i\n", size(T,1));
  for i=1:(size(T,1))
    fprintf(fileID, strT, [T(i,:) 0]);
  end
  fprintf(fileID,"End\n");
  fclose(fileID);

end

function checkObject(X,T)
  M = containers.Map;
  for i=1:size(T,1)
    x = X(T(i,:),:);
    centerX = floor(mean(x(:,1))*100000);
    centerY = floor(mean(x(:,2))*100000);
    centerZ = floor(mean(x(:,3))*100000);
    
    key = strcat(num2str(centerX),num2str(centerY),num2str(centerZ));
    M(key)=1;
  end
  if M.Count~=size(T,1)
    print("Duplicate element")
  end
  M = containers.Map;
  for i=1:size(X,1)
    x = X(i,:);
    centerX = floor(mean(x(:,1))*100000);
    centerY = floor(mean(x(:,2))*100000);
    centerZ = floor(mean(x(:,3))*100000);
    
    key = strcat(num2str(centerX),num2str(centerY),num2str(centerZ));
    M(key)=1;
  end
  if M.Count~=size(X,1)
    print("Duplicate node")
  end
end



function PBC = genPBC(pbcmap,X,T,exist)


pbcmap1 = pbcmap.pbcmap1;
pbcmap2 = pbcmap.pbcmap2;
% pbcmap3 = pbcmap.pbcmap3;

lelt=[];
for i=1:pbcmap1.elts
  idx1 = find(rot90(pbcmap1.a,-1)==i);
  idx2 = find(rot90(pbcmap1.b,-1)==i);
  lelt = [lelt; idx1 idx2];
end

for i=1:pbcmap2.elts
  idx1 = find(rot90(pbcmap2.a,-1)==i);
  idx2 = find(rot90(pbcmap2.b,-1)==i);
  lelt = [lelt; idx1 idx2];
end

% for i=1:pbcmap3.elts
%   idx1 = find(rot90(pbcmap3.a,-1)==i);
%   idx2 = find(rot90(pbcmap3.b,-1)==i);
%   lelt = [lelt; idx1 idx2];
% end

exist90 = rot90(exist,-1);
redmap1=[];
cpt=1;
for i=1:(size(exist90,1)*size(exist90,2))
  if exist90(i)==1
    redmap1=[redmap1 cpt];
    cpt=cpt+1;
  else
    redmap1=[redmap1 0];
  end
end
redmap=redmap1;



lnode = [];
for i=1:size(lelt,1)
  lnode = [lnode; T(redmap(lelt(i,1)),:)' T(redmap(lelt(i,2)),:)'];
end
%enlever doublons
dbls = [];
for i=1:size(lnode,1)
  dbls = [dbls 10000*lnode(i,1)+lnode(i,2)];
end
[~,idu] = unique(dbls);
lnode = lnode(idu,:);



PBC = [3*(lnode-1);3*(lnode-1)+1;3*(lnode-1)+2];


% vis periodic nodes here
%   plotHexMesh(X,T);
%   for i=1:size(lnode,1)
%     x1 = X(lnode(i,1),:);
%     x2 = X(lnode(i,2),:);
%     scatter3((x1(1)),(x1(2)),(x1(3)),"xr");
%     scatter3((x2(1)),(x2(2)),(x2(3)),"xr");
% 
%     pause
%   end



end


function [X5,T5]=genMesh(exist)


  [X,T] = createHEXmesh([-1 1 -1 1],size(exist,1),size(exist,2));
  % figure
  % scatter(X(:,1),X(:,2))
  % for i=1:size(X,1)
  %   text(X(i,1),X(i,2),num2str( i+size(X,1) )); hold on;
  % end

  X1 = [X zeros(size(X,1),1)];
  X2 = [X ones(size(X,1),1)*0.1];
  T1 = T;
  T2 = T + size(X,1);
  
  T = killVoid(exist,T);
  T1 = killVoid(exist,T1);
  T2 = killVoid(exist,T2);


  X5 = [X1;X2];
  T5 = [];
  for i=1:size(T,1)
  %   tmp = [];
    T5 = [T5; T1(i,:) T2(i,:)];
  end

  lexist = [];
  for i=1:size(X5,1)
    if ismember(i,T5)
      lexist=[lexist i];
    end
  end

  X6 = X5(lexist,:);

  T6 = T5;
  for i=1:( size(T6,1)*size(T6,2) )
    idx = find( T6(i)==lexist );
    T6(i)=idx;
  end



  X5=X6;
  T5=T6;
  X5(:,3) = X5(:,3)-0.1/2;

  checkObject(X5,T5)

  % dlmwrite("./X.csv",X5,'delimiter',',','precision',5)
  % dlmwrite("./T.csv",T5-1,'delimiter',',','precision',5)


  % genPLY(X5,T5)
  % genMESH(X5,T5)

end

function plotHexMesh(X,T)
  for i=1:size(T,1)
    x = X(T(i,:),:);
    fill3(x(:,1),x(:,2),x(:,3),'b'); hold on;
  end
  axis equal
end