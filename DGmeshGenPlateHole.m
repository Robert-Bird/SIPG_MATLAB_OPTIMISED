 function [coord,etpl_face,etpl,ed,BCx,BCy,diag]=DGmeshGenPlateHole(nelsR,nelsC)
% % clear
% nelsR = 2;
% nelsC = 2;

Ro = 10;           % Outer radius
Ri = 1;            % Inner Radius
nD = 2;
sfR=1;           % radial scalling factor
th = (pi/2)/nelsC;
nen=4;                       % Number of element nodes
tol=1e-6;
nels=nelsR*nelsC;            % Number of elements
nodes = nels*nen;            % Number of element nodes
coord=zeros(nodes,2);        % zero coordinates


rs=zeros(nelsR+1,1);
for i=1:(nelsR)   % Radius of the nodes along the edges
  rs(i+1)=rs(i)+sfR^(i-1);                                                  
end
rs=Ri+(Ro-Ri)*rs/rs(end);


k=1;
for nelR=1:nelsR;
    for nelC=1:nelsC
        th_1=th*(nelC-1);
        
        if abs(th_1)>tol && abs(th_1)<=(pi/4)
            Rmax_1=Ro/cos(th_1);
        elseif abs(th_1)>(pi/4) && abs(th_1-pi/2)>tol
            Rmax_1=Ro/cos(pi/2-th_1);
        else
            Rmax_1=Ro;
        end
        r_1=Ri+(Rmax_1-Ri)*(rs(nelR)-Ri)/(Ro-Ri);
        r_4=Ri+(Rmax_1-Ri)*(rs(nelR+1)-Ri)/(Ro-Ri);
        
        th_2=th*(nelC);
        if abs(th_2)>tol && abs(th_2)<=(pi/4)
            Rmax_2=Ro/cos(th_2);
        elseif abs(th_2)>(pi/4) && abs(th_2-pi/2)>tol
            Rmax_2=Ro/cos(pi/2-th_2);
        else
            Rmax_2=Ro;
        end
        r_2=Ri+(Rmax_2-Ri)*(rs(nelR)-Ri)/(Ro-Ri);
        r_3=Ri+(Rmax_2-Ri)*(rs(nelR+1)-Ri)/(Ro-Ri);
        coord(k  ,:)=r_1*[cos(th_1) sin(th_1)];                             
        coord(k+1,:)=r_2*[cos(th_2) sin(th_2)];                                  
        coord(k+2,:)=r_3*[cos(th_2) sin(th_2)];                                  
        coord(k+3,:)=r_4*[cos(th_1) sin(th_1)];                                  
        k=k+4;
    end
end
% coord(:,2)=coord(:,2).*(5/10);

%element topology
etpl= 1:nodes;
etpl= reshape(etpl,nen,nels)';
ed = 1:(nodes*nD);
ed = reshape(ed,nen*nD,nels)';


% Face information 
% postive element numbers for 2.
el_2_positive = reshape((repmat((1:(nelsC-1)),nelsR,1) +repmat(nelsC*[0:(nelsR-1)]',1,(nelsC-1)))',((nelsC-1)*nelsR),1);
el_4_negative = reshape((repmat((2:(nelsC)),nelsR,1)   +repmat(nelsC*[0:(nelsR-1)]',1,(nelsC-1)))',((nelsC-1)*nelsR),1);
el_3_positive = reshape((repmat(1:(nelsC),(nelsR-1),1) + repmat([nelsC:nelsC:((nelsR-1)*(nelsC))]'-nelsC,1,nelsC))',((nelsR-1)*nelsC),1);
el_1_negative = el_3_positive+nelsC;

% nodes necessary for postive face 2 connectivity
nodes_pve_3 = etpl(el_3_positive,[3,4])'; coord_pve_3 = coord(nodes_pve_3,:);
nodes_pve_2 = etpl(el_2_positive,[2,3])'; coord_pve_2 = coord(nodes_pve_2,:);

normalised = (((coord_pve_3(1:2:end,1)-coord_pve_3(2:2:end,1)).^2 )+ ((coord_pve_3(1:2:end,2)-coord_pve_3(2:2:end,2)).^2)).^(0.5);
coord_3_vet=[(coord_pve_3(1:2:end,1)-coord_pve_3(2:2:end,1))./normalised,(coord_pve_3(1:2:end,2)-coord_pve_3(2:2:end,2))./normalised];
normal_3 = [0 1;-1 0]*coord_3_vet';
h_3=normalised;

normalised = (((coord_pve_2(1:2:end,1)-coord_pve_2(2:2:end,1)).^2 )+ ((coord_pve_2(1:2:end,2)-coord_pve_2(2:2:end,2)).^2)).^(0.5);
coord_2_vet=[(coord_pve_2(1:2:end,1)-coord_pve_2(2:2:end,1))./normalised,(coord_pve_2(1:2:end,2)-coord_pve_2(2:2:end,2))./normalised];
normal_2 = [0 1;-1 0]*coord_2_vet';
h_2=normalised;
% Boundary conditions only apply to the edges ones.
ext = [(1:nelsC:(nelsR*nelsC)),((1:nelsC:(nelsR*nelsC)) + nelsC -1)];
ext_normal = [repmat([0 ,-1],nelsR,1);repmat([-1 ,0],nelsR,1)];

positive_elements = [el_2_positive; el_3_positive; ext']; 
positive_face     = [2*ones(size(el_2_positive,1),1);3*ones(size(el_3_positive,1),1);4*ones(nelsR,1);2*ones(nelsR,1)];
negative_elements = [el_4_negative; el_1_negative; zeros(size(ext,2),1)];
negative_face     = [4*ones(size(el_2_positive,1),1);1*ones(size(el_3_positive,1),1);zeros(nelsR,1);zeros(nelsR,1)];
etpl_face = [positive_elements,negative_elements,positive_face, negative_face,[normal_2';normal_3';ext_normal],[h_2;h_3;(rs(2:end)-rs(1:(end-1)));(rs(2:end)-rs(1:(end-1)))]];

