function [k]=Linear_2D(nen,nels,num_faces,nD,etpl,etpl_face,ed,coord,E,v,nelblo_size)


%% ==================== Problem Defining Parameters ==================== %%
ndof=nen*nD; nndof=ndof^2; tndof=ndof*nels;                                % Setup of problem parameters
[vwp,vdNr,ngp,swp,sdNr,nsgp,sNr,~]=dershapefunc2D(nen,0);                    % Shape function production
if nen == 4
    p=1;
elseif nen == 8
    p=2;
end

[no_blocks,nelblo_all]=block_size(nels,nelblo_size);                       % Block size determination
K_all=zeros(nels,ndof,ndof);                                               % Assigning memory for global stiffness matrix
K_eig_all=zeros(nels,ndof,ndof);                                           % Global eigenvalue store

const = E/(1-(v^2)); B = const*(v); A = const*(1); C = E/(2*(1+v));        % 2D plane stress linear elastic stiffness matrix values
pen = 10*(p^2)*E/((1-(v^2)));                                              % dG penalty term for 2D linear elasticity 


%% ========================= Volume Integration ======================== %%
for bloc_count = 1:ceil(no_blocks)                                         % Block loop over all blocks
    nelblo = nelblo_all(bloc_count);                                       % Defining a variable nelblo to number of elements within a block (for readability)
    invJx = zeros(nelblo,nD); invJy=invJx;                                 % Assigning memory for variables within the loop.
    K_block = zeros(nelblo,ndof,ndof);                                     % Temporary storage variable for elements in the current block loop
    K_eig = zeros(nelblo,ndof,ndof);
    indx_coord = (1:nelblo)+((bloc_count-1)*nelblo_all(1));                % Element numbers in the current block loop
    
    Ecoord_x = reshape(coord(etpl(indx_coord,:),1),nelblo,nen);            % x coordinates of element node
    Ecoord_y = reshape(coord(etpl(indx_coord,:),2),nelblo,nen);            % y coordinates of element node
    
    for gp = 1:ngp                                                         % Start of Gauss point intergration loop for specific block
        indx=nD*(gp-1)+(1:nD);                                             % indx for the specific guass point loop
        Jx = vdNr(indx,:)*Ecoord_x';                                       % Matricies of Jacbian terms    | Jx(1,el_no), Jy(1,el_no) |  (el_no = element number)
        Jy = vdNr(indx,:)*Ecoord_y';                                       %                               | Jx(2,el_no), Jy(2,el_no) |
        Jdet = 1./((Jx(1,:).*Jy(2,:))-(Jx(2,:).*Jy(1,:)));                 % Element determinant
        invJx(:,1) = Jdet.* Jy(2,:);     invJx(:,2) = Jdet.*-Jy(1,:);      % Matrices of inverse Jacobian terms | invJx(el_no,1), invJx(el_no,2)|
        invJy(:,1) = Jdet.*-Jx(2,:);     invJy(:,2) = Jdet.* Jx(1,:);      %                                    | invJy(el_no,1), invJy(el_no,2)|
        dNx = (invJx * vdNr(indx,:));                                      % Matrices of global shape function derivatives   | dNx(el_no,1),... ,dNx(el_no,nen) |
        dNy = (invJy * vdNr(indx,:));                                      %                                                 | dNy(el_no,1),... ,dNy(el_no,nen) |
        w=((1./Jdet).*vwp(gp))';                                           % Specific weight calcuation for each element within the block.
        
        % Local stiffness matrix calculation
        for i = 1:nen % Row index
            Ai = (i-1)*nD;
            for j = i:nen % Column index
                Aj = (j-1)*nD;
                % Only upper triangle of matrix is calculated. The transpose (lower triangle) is later added, so leading diagonals need to be halved.
                if j==i; half=0.5; else half =1; end
                K_block(:,1+Ai,1+Aj)  = K_block(:,1+Ai,1+Aj) + (dNx(:,i).*dNx(:,j)*A + dNy(:,i).*dNy(:,j)*C).*w*half;
                K_block(:,1+Ai,2+Aj)  = K_block(:,1+Ai,2+Aj) + (dNx(:,i).*dNy(:,j)*B + dNy(:,i).*dNx(:,j)*C).*w;
                K_block(:,2+Ai,2+Aj)  = K_block(:,2+Ai,2+Aj) + (dNy(:,i).*dNy(:,j)*A + dNx(:,i).*dNx(:,j)*C).*w*half;
                if j>i;
                    K_block(:,2+Ai,1+Aj)  = K_block(:,2+Ai,1+Aj) + (dNy(:,i).*dNx(:,j)*B + dNx(:,i).*dNy(:,j)*C).*w;
                end
            end
        end
        
    end
    
    K_all(indx_coord,:,:)=K_block;
    K_eig_all(indx_coord,:,:)=K_eig;
    clear K_block K_eig
end

K_all = reshape(reshape(K_all,nels,(ndof)^2)',nels*(ndof)^2,1);
ed_i = reshape(repmat(ed,1,ndof)',1,ndof*ndof*nels);                       % A vector of the 'i' (row) locations of corressponding values in K_all vector
ed_j = reshape(ed',1,ndof*nels);
ed_j = reshape(repmat(ed_j,ndof,1),1,ndof*ndof*nels);                      % A vector of the 'j' (column) locations of corressponding values in K_all vector

clear Ai Aj Ecoord_x Ecoord_y Jdet Jx Jy a ans bloc_count const dNx dNy ... 
    gp half  i indx indx_coord invJx  invJy j nelblo nelblo_all no_blocks p w          


%% =================== Surface block loop setup =========================%%


etpl_iface = etpl_face((etpl_face(:,2)~=0),:);                             % Matrix of surface to surface connectivity information (internal)
tot_f=size(etpl_iface,1);                                                  % Number of inward surfaces
[no_blocks_i,nelblo_all_i]=block_size(tot_f,nelblo_size);                  % interal surface: no. of block loops and size of the corresponding face block (no_blocks_i,nelblo_all_i respectively)

glob_1=zeros(tot_f,ndof,ndof);                                             % Assigning memory for all surface integrals of the global stiffness matrices
glob_2=glob_1; glob_3=glob_1; glob_4=glob_1;

num_blocks = no_blocks_i;                                                  % Number of block loops

ed_p=ed(etpl_iface(:,1),:);                                                % +ve face steering matrix
ed_n=ed(etpl_iface(:,2),:);                                                % -ve face steering matrix

%% =================== Internal face Jacobian sort ======================%%

%% ================== Surface integral ================================= %%
for Block_n = 1:num_blocks
    %% ================== Surface integral setup (internal) ============ %%
    if Block_n <= no_blocks_i
        si = nelblo_all_i(Block_n);                                        % Number of surfaces in internal block loop
        block_index_i = (1:si)+(nelblo_all_i(1)*(Block_n-1));              % Surface numbers in internal block loop
        etpl_face_block = etpl_iface(block_index_i,:);                         % Matrix of, surface to surface element information in current block loop
        
        Ecoord_x_p = reshape(coord(etpl(etpl_face_block(:,1),:),1),si,nen);    % x coordinates of +ve element nodes
        Ecoord_y_p = reshape(coord(etpl(etpl_face_block(:,1),:),2),si,nen);    % y coordinates of +ve element nodes
        
        Ecoord_x_n = reshape(coord(etpl(etpl_face_block(:,2),:),1),si,nen);    % x coordinates of -ve element nodes
        Ecoord_y_n = reshape(coord(etpl(etpl_face_block(:,2),:),2),si,nen);    % y coordinates of -ve element nodes
        
        
        Jx_p = zeros(nD,si);    Jy_p = Jx_p;                               % +ve element Jacobian components
        Jx_n = Jx_p;            Jy_n = Jx_p;                               % -ve element Jacobian components
        invJx_p = Jx_p';        invJy_p = Jx_p';                           % +ve inverse Jacobian components
        invJx_n = Jx_p';        invJy_n = Jx_p';                           % -ve inverse Jacobian components
        dNx_p = zeros(si,nen);  dNy_p = dNx_p;                             % +ve global shape function derivative components
        dNx_n = dNx_p;          dNy_n = dNx_p;                             % -ve global shape function derivative components
        sNr_p = zeros(si,nen);  sNr_n = zeros(si,nen);                     % +ve/ive shape function components
        % Internal surface integral weight
        
        glob_pp = zeros(si,ndof,ndof); glob_pn = glob_pp;     glob_np = glob_pp;        glov_nn = glob_pp;    %
        
        nx_i = etpl_face_block(:,5);                                           % Unit normals for the external surface intergral
        ny_i = etpl_face_block(:,6);                                           %
    end
    
    
    
    %% =================== Surface integration Gauss point loop========= %%
    for gp = 1:nsgp
        %% =================== Surface number specific calculations ==== %%
        for fn = 1:num_faces
            indx_dNr =((gp-1)*nD)+((fn-1)*nD*nsgp)+(1:nD);                 % Referenece shape function derivative index for current gauss point and surface number
            indx_sNr = gp+((fn-1)*nsgp);                                   % Referenece shape function index for current gauss point and surface number
            %% =================== Internal dNr ======================== %%
            if Block_n <= no_blocks_i
                index_p=(etpl_face_block(:,3)==fn);                            % Logic statement to select +ve element surface integrals over the face number 'fn'
                index_n=(etpl_face_block(:,4)==fn);                            % Logic statement to select -ve element surface integrals over the face number 'fn'
                
                if sum(index_p)~=0
                    
                    Jx_p(:,index_p) = sdNr(indx_dNr,:)*Ecoord_x_p(index_p,:)'; % Jacobian calculated from non-reversed Gauss point ordering
                    Jy_p(:,index_p) = sdNr(indx_dNr,:)*Ecoord_y_p(index_p,:)'; %
                    
                    det = (Jx_p(1,index_p).*Jy_p(2,index_p))-(Jx_p(2,index_p).*Jy_p(1,index_p)); % Jacobian determinant calculation
                    det = 1./det';
                    
                    invJx_p(index_p,:) = [ det.*(Jy_p(2,index_p))',-det.*(Jy_p(1,index_p))']; % Inverse Jacobian  calculation
                    invJy_p(index_p,:) = [-det.*(Jx_p(2,index_p))', det.*(Jx_p(1,index_p))'];
                    
                    dNx_p(index_p,:) = invJx_p(index_p,:)*sdNr(indx_dNr,:); % Global shape function calculation
                    dNy_p(index_p,:) = invJy_p(index_p,:)*sdNr(indx_dNr,:); %
                    
                    sNr_p(index_p,:) = repmat(sNr(indx_sNr,:),sum(index_p),1);% Surface number specific shape function for +ve elements
                end
                
                if sum(index_n)~=0
                    indx_dNr = ((nsgp-gp)*nD)+((fn-1)*nD*nsgp)+(1:nD);     % Referenece shape function derivative index for current gauss point and surface number
                    indx_sNr = (nsgp-gp+1)+((fn-1)*nsgp);                  % Referenece shape function index for current gauss point and surface number
                    Jx_n(:,index_n) = sdNr(indx_dNr,:)*Ecoord_x_n(index_n,:)'; % Jacobian Calculation
                    Jy_n(:,index_n) = sdNr(indx_dNr,:)*Ecoord_y_n(index_n,:)'; %
                    det = (Jx_n(1,index_n).*Jy_n(2,index_n))-(Jx_n(2,index_n).*Jy_n(1,index_n));% Jacobian determinant calculation
                    det = 1./det';
                    
                    invJx_n(index_n,:) = [ det.*(Jy_n(2,index_n))',-det.*(Jy_n(1,index_n))']; % Inverse Jacobian  calculation
                    invJy_n(index_n,:) = [-det.*(Jx_n(2,index_n))', det.*(Jx_n(1,index_n))']; %
                    
                    dNx_n(index_n,:) = invJx_n(index_n,:)*sdNr(indx_dNr,:); % Global shape function calculation
                    dNy_n(index_n,:) = invJy_n(index_n,:)*sdNr(indx_dNr,:); %
                    
                    sNr_n(index_n,:) = repmat(sNr(indx_sNr,:),sum(index_n),1); % Surface number specific shape function for -ve elements
                    
                end
            end
        end
        
        weighti = repmat(pen,si,1)./((etpl_face_block(:,end)));                % Penalty internal integral weight (Penalty/element height)
        gwi = repmat(swp(gp)/2,si,1);                                      % Gauss point weight (internal)
        Jwp = etpl_face_block(:,end)/2;
        for i = 1:nen
            Ai = (i-1)*nD;
            
            % j2 is an index to count backwards through the local matrix
            % columns of D.
            j2=(nen+1)-i;
            Bj = (j2-1)*nD;
            
            %% =================== Calculation of repeated terms ======= %%
            % To reduce the number of total calculations, components of the
            % surface stiffness matrix have been taken out of the j=i:nen
            % loop.
            % =================== +ve components of C =================== %
            C2t_11 = ((dNx_p(:,i).*nx_i.*A) + (dNy_p(:,i).*ny_i.*C)).*gwi;
            C2t_12 = ((dNx_p(:,i).*ny_i.*B) + (dNy_p(:,i).*nx_i.*C)).*gwi;
            C2t_21 = ((dNy_p(:,i).*nx_i.*B) + (dNx_p(:,i).*ny_i.*C)).*gwi;
            C2t_22 = ((dNx_p(:,i).*nx_i.*C) + (dNy_p(:,i).*ny_i.*A)).*gwi;
            
            % =================== -ve components of C =================== %
            C2t_11n = ((dNx_n(:,i).*nx_i.*A) + (dNy_n(:,i).*ny_i.*C)).*gwi;
            C2t_12n = ((dNx_n(:,i).*ny_i.*B) + (dNy_n(:,i).*nx_i.*C)).*gwi;
            C2t_21n = ((dNy_n(:,i).*nx_i.*B) + (dNx_n(:,i).*ny_i.*C)).*gwi;
            C2t_22n = ((dNx_n(:,i).*nx_i.*C) + (dNy_n(:,i).*ny_i.*A)).*gwi;
            
            
            % =================== +ve components of D =================== %
            D2t_11 = ((dNx_p(:,j2).*nx_i.*A) + (dNy_p(:,j2).*ny_i.*C)).*gwi;
            D2t_12 = ((dNy_p(:,j2).*nx_i.*B) + (dNx_p(:,j2).*ny_i.*C)).*gwi;
            D2t_21 = ((dNx_p(:,j2).*ny_i.*B) + (dNy_p(:,j2).*nx_i.*C)).*gwi;
            D2t_22 = ((dNx_p(:,j2).*nx_i.*C) + (dNy_p(:,j2).*ny_i.*A)).*gwi;
            
            % =================== -ve components of D =================== %
            D2t_11n = ((dNx_n(:,j2).*nx_i.*A) + (dNy_n(:,j2).*ny_i.*C)).*gwi;
            D2t_12n = ((dNy_n(:,j2).*nx_i.*B) + (dNx_n(:,j2).*ny_i.*C)).*gwi;
            D2t_21n = ((dNx_n(:,j2).*ny_i.*B) + (dNy_n(:,j2).*nx_i.*C)).*gwi;
            D2t_22n = ((dNx_n(:,j2).*nx_i.*C) + (dNy_n(:,j2).*ny_i.*A)).*gwi;
            
            for j = i:nen
                Aj = (j-1)*nD;
                
                % i2 is index to count backwards through the local matrix
                % rows of D.
                i2=(nen+1)-j;
                Bi = (i2-1)*nD;
                if j==i; half=0.5; else half =1; end                       % As the transpose is added, the components on the leading diagonal are halved
                
                %% ================== Internal surface stiffness calcs = %%
                if Block_n <= no_blocks_i
                    
                    %========================================C1===============================%
                    glob_pp(:,1+Ai,1+Aj) = glob_pp(:,1+Ai,1+Aj) + C2t_11.*sNr_p(:,j).*Jwp*half;
                    glob_pp(:,1+Ai,2+Aj) = glob_pp(:,1+Ai,2+Aj) + C2t_12.*sNr_p(:,j).*Jwp;
                    glob_pp(:,2+Ai,2+Aj) = glob_pp(:,2+Ai,2+Aj) + C2t_22.*sNr_p(:,j).*Jwp*half;
                    
                    %========================================C2===============================%
                    glob_pn(:,1+Ai,1+Aj) = glob_pn(:,1+Ai,1+Aj) - C2t_11.*sNr_n(:,j).*Jwp*half;
                    glob_pn(:,1+Ai,2+Aj) = glob_pn(:,1+Ai,2+Aj) - C2t_12.*sNr_n(:,j).*Jwp;
                    glob_pn(:,2+Ai,2+Aj) = glob_pn(:,2+Ai,2+Aj) - C2t_22.*sNr_n(:,j).*Jwp*half;
                    
                    %========================================C3===============================%
                    glob_np(:,1+Ai,1+Aj) = glob_np(:,1+Ai,1+Aj) + C2t_11n.*sNr_p(:,j).*Jwp*half;
                    glob_np(:,1+Ai,2+Aj) = glob_np(:,1+Ai,2+Aj) + C2t_12n.*sNr_p(:,j).*Jwp;
                    glob_np(:,2+Ai,2+Aj) = glob_np(:,2+Ai,2+Aj) + C2t_22n.*sNr_p(:,j).*Jwp*half;
                    
                    %========================================C4===============================%
                    glov_nn(:,1+Ai,1+Aj) = glov_nn(:,1+Ai,1+Aj) - C2t_11n.*sNr_n(:,j).*Jwp*half;
                    glov_nn(:,1+Ai,2+Aj) = glov_nn(:,1+Ai,2+Aj) - C2t_12n.*sNr_n(:,j).*Jwp;
                    glov_nn(:,2+Ai,2+Aj) = glov_nn(:,2+Ai,2+Aj) - C2t_22n.*sNr_n(:,j).*Jwp*half;
                    
                    if j>i
                        glob_pp(:,2+Ai,1+Aj) = glob_pp(:,2+Ai,1+Aj) + C2t_21.*sNr_p(:,j).*Jwp;
                        glob_pn(:,2+Ai,1+Aj) = glob_pn(:,2+Ai,1+Aj) - C2t_21.*sNr_n(:,j).*Jwp;
                        glob_np(:,2+Ai,1+Aj) = glob_np(:,2+Ai,1+Aj) + C2t_21n.*sNr_p(:,j).*Jwp;
                        glov_nn(:,2+Ai,1+Aj) = glov_nn(:,2+Ai,1+Aj) - C2t_21n.*sNr_n(:,j).*Jwp;
                    end
                    
                    %========================================D1===============================%
                    glob_pp(:,1+Bi,1+Bj) = glob_pp(:,1+Bi,1+Bj) + D2t_11.*sNr_p(:,i2).*Jwp*half;
                    glob_pp(:,1+Bi,2+Bj) = glob_pp(:,1+Bi,2+Bj) + D2t_12.*sNr_p(:,i2).*Jwp;
                    glob_pp(:,2+Bi,2+Bj) = glob_pp(:,2+Bi,2+Bj) + D2t_22.*sNr_p(:,i2).*Jwp*half;
                    
                    %========================================D2===============================%
                    glob_pn(:,1+Bi,1+Bj) = glob_pn(:,1+Bi,1+Bj) + D2t_11n.*sNr_p(:,i2).*Jwp*half;
                    glob_pn(:,1+Bi,2+Bj) = glob_pn(:,1+Bi,2+Bj) + D2t_12n.*sNr_p(:,i2).*Jwp;
                    glob_pn(:,2+Bi,2+Bj) = glob_pn(:,2+Bi,2+Bj) + D2t_22n.*sNr_p(:,i2).*Jwp*half;
                    
                    %========================================D3===============================%
                    glob_np(:,1+Bi,1+Bj) = glob_np(:,1+Bi,1+Bj) - D2t_11.*sNr_n(:,i2).*Jwp*half;
                    glob_np(:,1+Bi,2+Bj) = glob_np(:,1+Bi,2+Bj) - D2t_12.*sNr_n(:,i2).*Jwp;
                    glob_np(:,2+Bi,2+Bj) = glob_np(:,2+Bi,2+Bj) - D2t_22.*sNr_n(:,i2).*Jwp*half;
                    
                    %========================================D4===============================%
                    glov_nn(:,1+Bi,1+Bj) = glov_nn(:,1+Bi,1+Bj) - D2t_11n.*sNr_n(:,i2).*Jwp*half;
                    glov_nn(:,1+Bi,2+Bj) = glov_nn(:,1+Bi,2+Bj) - D2t_12n.*sNr_n(:,i2).*Jwp;
                    glov_nn(:,2+Bi,2+Bj) = glov_nn(:,2+Bi,2+Bj) - D2t_22n.*sNr_n(:,i2).*Jwp*half;
                    
                    if j2>i2
                        glob_pp(:,2+Bi,1+Bj) = glob_pp(:,2+Bi,1+Bj) + D2t_21.*sNr_p(:,i2).*Jwp;
                        glob_pn(:,2+Bi,1+Bj) = glob_pn(:,2+Bi,1+Bj) + D2t_21n.*sNr_p(:,i2).*Jwp;
                        glob_np(:,2+Bi,1+Bj) = glob_np(:,2+Bi,1+Bj) - D2t_21.*sNr_n(:,i2).*Jwp;
                        glov_nn(:,2+Bi,1+Bj) = glov_nn(:,2+Bi,1+Bj) - D2t_21n.*sNr_n(:,i2).*Jwp;
                    end
                    
                    %========================================E1===============================%
                    glob_pp(:,1+Ai,1+Aj) = glob_pp(:,1+Ai,1+Aj) - sNr_p(:,i).*sNr_p(:,j).*weighti.*Jwp.*gwi*half;
                    glob_pp(:,2+Ai,2+Aj) = glob_pp(:,2+Ai,2+Aj) - sNr_p(:,i).*sNr_p(:,j).*weighti.*Jwp.*gwi*half;
                    
                    %========================================E4===============================%
                    glob_pn(:,1+Ai,1+Aj) = glob_pn(:,1+Ai,1+Aj) + sNr_p(:,i).*sNr_n(:,j).*weighti.*Jwp.*gwi*half;
                    glob_pn(:,2+Ai,2+Aj) = glob_pn(:,2+Ai,2+Aj) + sNr_p(:,i).*sNr_n(:,j).*weighti.*Jwp.*gwi*half;
                    
                    %========================================E2===============================%
                    glob_np(:,1+Ai,1+Aj) = glob_np(:,1+Ai,1+Aj) + sNr_n(:,i).*sNr_p(:,j).*weighti.*Jwp.*gwi*half;
                    glob_np(:,2+Ai,2+Aj) = glob_np(:,2+Ai,2+Aj) + sNr_n(:,i).*sNr_p(:,j).*weighti.*Jwp.*gwi*half;
                    
                    %========================================E3===============================%
                    glov_nn(:,1+Ai,1+Aj) = glov_nn(:,1+Ai,1+Aj) - sNr_n(:,i).*sNr_n(:,j).*weighti.*Jwp.*gwi*half;
                    glov_nn(:,2+Ai,2+Aj) = glov_nn(:,2+Ai,2+Aj) - sNr_n(:,i).*sNr_n(:,j).*weighti.*Jwp.*gwi*half;
                    
                    
                end
            end
        end
        
    end
    
    
    %% ==================  Storage of temporary stiffness block ======== %%
    if Block_n <= no_blocks_i
        glob_1(block_index_i,:,:) = glob_pp;
        glob_2(block_index_i,:,:) = glob_pn;
        glob_3(block_index_i,:,:) = glob_np;
        glob_4(block_index_i,:,:) = glov_nn;
        clear glob_pp glob_pn glob_np glob_nn
    end
end

clear_many_vars






%% ================== Reshaping of stiffness arrays ==================== %%
% Reshaping the storage arrays into a vector form for storage into a
% sparse matrix.
glob_1_rs = reshape(reshape(glob_1,tot_f,nndof)',tot_f*(nndof),1);
glob_2_rs = reshape(reshape(glob_2,tot_f,nndof)',tot_f*(nndof),1);
glob_3_rs = reshape(reshape(glob_3,tot_f,nndof)',tot_f*(nndof),1);
glob_4_rs = reshape(reshape(glob_4,tot_f,nndof)',tot_f*(nndof),1);


%% ================== Producing steering vectors ======================= %%
pos_i = reshape(repmat(ed_p,1,ndof)',1,nndof*tot_f);                          %
ed_pve = reshape(ed_p',1,ndof*tot_f);                                      % +ve element steering vector for row 'j' and column 'i'
pos_j = reshape(repmat(ed_pve,ndof,1),1,tot_f*nndof);                         %

neg_i = reshape(repmat(ed_n,1,ndof)',1,nndof*tot_f);                          %
ed_nve = reshape(ed_n',1,ndof*tot_f);                                      % -ve element steering vector for row 'j' and column 'i'
neg_j = reshape(repmat(ed_nve,ndof,1),1,tot_f*nndof);                         %

%% ================== Sparse storage =================================== %%

%% ================== Determining block loop size ====================== %%
% No. of block loops and size of the corresponding element block (no_blocks,nelblo_all respectively), for:
k = sparse(ndof*nels,ndof*nels);

% =================== Volumetric sparse storage ========================= %
k = k + sparse(ed_i,ed_j,K_all,tndof,tndof);
% =================== Surface sparse storage ============================ %
k = k -  sparse(pos_i,pos_j,glob_1_rs,tndof,tndof);
k = k -  sparse(pos_i,neg_j,glob_2_rs,tndof,tndof);
k = k -  sparse(neg_i,pos_j,glob_3_rs,tndof,tndof);
k = k -  sparse(neg_i,neg_j,glob_4_rs,tndof,tndof);
k=k+k';
end


