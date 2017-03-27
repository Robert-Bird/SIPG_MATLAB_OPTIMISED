function [coord,etpl_face,etpl,ed]=DGmeshGenPlateHole_8node(nelsR,nelsC)
% nelsR= 2; nelsC=2;
%% parameters controling the mesh
Ro = 10;                                                                   % Outer radius
Ri = 1;                                                                    % Inner Radius
nD = 2;                                                                   % No. of dimensions
sfR = 1;                                                                   % Radial scalling factor

%% initial information for mesh generation
nen = 8;                                                                   % Number of element nodes
nels = nelsR*nelsC;                                                        % Number of elements
nodes = nels*nen;                                                          % Number of element nodes
coord = zeros(nodes,2);                                                    % Zero coordinates
th = (pi/2)/nelsC;                                                         % Angle around the circumference
tol=1e-6;
rs=zeros(nelsR+1,1);
for i=1:(nelsR)
    rs(i+1)=rs(i)+sfR^(i-1);                                                 % radius scaling factor
end
rs=Ri+(Ro-Ri)*rs/rs(end);                                                  % nodal radii


k=1;
%% coordinate generation
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
        r_7=Ri+(Rmax_1-Ri)*(rs(nelR+1)-Ri)/(Ro-Ri);
        
        th_2=th*(nelC);
        if abs(th_2)>tol && abs(th_2)<=(pi/4)
            Rmax_2=Ro/cos(th_2);
        elseif abs(th_2)>(pi/4) && abs(th_2-pi/2)>tol
            Rmax_2=Ro/cos(pi/2-th_2);
        else
            Rmax_2=Ro;
        end
        r_3=Ri+(Rmax_2-Ri)*(rs(nelR)-Ri)/(Ro-Ri);
        r_5=Ri+(Rmax_2-Ri)*(rs(nelR+1)-Ri)/(Ro-Ri);
        
        coord(k  ,:)=r_1*[cos(th_1) sin(th_1)];%1
        coord(k+2,:)=r_3*[cos(th_2) sin(th_2)];%3
        coord(k+4,:)=r_5*[cos(th_2) sin(th_2)];%5
        coord(k+6,:)=r_7*[cos(th_1) sin(th_1)];%7
        
        coord(k+1,:)=(coord(k+2,:)+coord(k  ,:))./2; %2
        coord(k+3,:)=(coord(k+2,:)+coord(k+4,:))./2;%4
        coord(k+5,:)=(coord(k+4,:)+coord(k+6,:))./2;
        coord(k+7,:)=(coord(k  ,:)+coord(k+6,:))./2;
        
        k=k+8;
    end
end

%% Element topology
etpl= 1:nodes;
etpl= reshape(etpl,nen,nels)';
ed = 1:(nodes*nD);
ed = reshape(ed,nen*nD,nels)';
% hold on
% for nel=1:nelsC*nelsR
%     plot(coord(etpl(nel,:),1),coord(etpl(nel,:),2),'bx-');                  % Plotting original element map
%     plot(x(nel,:),y(nel,:),'-r');                                          % Plotting resultant element map
% end
%% Face connectivity
el_2_positive = reshape((repmat((1:(nelsC-1)),nelsR,1) +repmat(nelsC*[0:(nelsR-1)]',1,(nelsC-1)))',((nelsC-1)*nelsR),1);
el_4_negative = reshape((repmat((2:(nelsC)),nelsR,1)   +repmat(nelsC*[0:(nelsR-1)]',1,(nelsC-1)))',((nelsC-1)*nelsR),1);
el_3_positive = reshape((repmat(1:(nelsC),(nelsR-1),1) + repmat([nelsC:nelsC:((nelsR-1)*(nelsC))]'-nelsC,1,nelsC))',((nelsR-1)*nelsC),1);
el_1_negative = el_3_positive+nelsC;

%% Face perpendicular normal and heights
nodes_pve_3 = etpl(el_3_positive,[5,7])'; coord_pve_3 = coord(nodes_pve_3,:);
nodes_pve_2 = etpl(el_2_positive,[3,5])'; coord_pve_2 = coord(nodes_pve_2,:);

height = (((coord_pve_3(1:2:end,1)-coord_pve_3(2:2:end,1)).^2 )+ ((coord_pve_3(1:2:end,2)-coord_pve_3(2:2:end,2)).^2)).^(0.5);
coord_3_vet=[(coord_pve_3(1:2:end,1)-coord_pve_3(2:2:end,1))./height,(coord_pve_3(1:2:end,2)-coord_pve_3(2:2:end,2))./height];
normal_3 = [0 1; -1 0]*coord_3_vet';                                       % +ve element face 3 normal
h_3=height;                                                                % +ve element face 3 height

height = (((coord_pve_2(1:2:end,1)-coord_pve_2(2:2:end,1)).^2 )+ ((coord_pve_2(1:2:end,2)-coord_pve_2(2:2:end,2)).^2)).^(0.5);
coord_2_vet=[(coord_pve_2(1:2:end,1)-coord_pve_2(2:2:end,1))./height,(coord_pve_2(1:2:end,2)-coord_pve_2(2:2:end,2))./height];
normal_2 = [0 1; -1 0]*coord_2_vet';                                       % +ve element face 2 normal
h_2=height;                                                                % +ve element face 2 height

%% +ve and -ve element numbers
positive_elements = [el_2_positive; el_3_positive];
positive_face     = [2*ones(size(el_2_positive,1),1);3*ones(size(el_3_positive,1),1)];
negative_elements = [el_4_negative; el_1_negative];
negative_face     = [4*ones(size(el_2_positive,1),1);1*ones(size(el_3_positive,1),1)];
etpl_face = [positive_elements,negative_elements,positive_face, negative_face,[normal_2';normal_3'],[h_2;h_3]];

