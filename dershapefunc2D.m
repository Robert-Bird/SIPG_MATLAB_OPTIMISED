function [vwp,vdNr,ngp,swp,sdNr,nsgp,sNr,Nr]=dershapefunc2D(nen,stress_flag)
if nen == 4 && stress_flag == 0
    nf = 4;                                                                % Number of faces
    ngp=4;                                                                 % Number of area Gauss points
    nsgp = 3;                                                              % Number of surface Gauss points
    nD = 2;                                                                % Number of dimensions
    g2=1/sqrt(3);                                                          % Area Gauss point position (magnitude)
    
    %% ========= 4 node quad shape function and derivative values ========== %%
    gp = [-1 -1  1  1; -1  1  1 -1 ]'.*g2;                                 % Area Gauss point position
    vwp=ones(4,1);                                                         % Area Guass point weights
    r = (ngp*nD);
    xsi=gp(:,1);  eta=gp(:,2); 	                                           % Reference coodinate position
    vdNr(1:2:r,1)  = -(1 - eta)./4;                                        % Area shape function derivative values
    vdNr(1:2:r,2)  = -(1 + eta)./4;
    vdNr(1:2:r,3)  =  (1 + eta)./4;
    vdNr(1:2:r,4)  =  (1 - eta)./4;
    vdNr(2:2:r,1)  = -(1 - xsi)./4;
    vdNr(2:2:r,2)  =  (1 - xsi)./4;
    vdNr(2:2:r,3)  =  (1 + xsi)./4;
    vdNr(2:2:r,4)  = -(1 + xsi)./4;
    r = r/nD;
    Nr(1:r,1) = (1 - xsi).*(1 - eta)./4;                                   % Area shape function values
    Nr(1:r,2) = (1 - xsi).*(1 + eta)./4;
    Nr(1:r,3) = (1 + xsi).*(1 + eta)./4;
    Nr(1:r,4) = (1 + xsi).*(1 - eta)./4;
    
    g1=sqrt(3/5);                                                          % Surface Gauss point position (magnitude)
    gp = [-1  -1  -1 -g1  0  g1  1   1   1  g1  0  -g1 ;                   % Surface Gauss point position
        -g1  0   g1   1  1  1   g1  0  -g1 -1  -1  -1]';                  %
    w1 = (5/9); w2 = (8/9);                                                % Edge Guass point weight
    swp = [w1 w2 w1];
    r = nsgp*nf*nD;
    xsi = gp(:,1); eta=gp(:,2);                                            % Reference coodinate position
    sdNr(1:2:r,1)  = -(1 - eta)./4;                                        % Edge shape function derivatives
    sdNr(1:2:r,2)  = -(1 + eta)./4;
    sdNr(1:2:r,3)  =  (1 + eta)./4;
    sdNr(1:2:r,4)  =  (1 - eta)./4;
    sdNr(2:2:r,1)  = -(1 - xsi)./4;
    sdNr(2:2:r,2)  =  (1 - xsi)./4;
    sdNr(2:2:r,3)  =  (1 + xsi)./4;
    sdNr(2:2:r,4)  = -(1 + xsi)./4;
    r=r/nD;
    sNr(1:r,1) = (1-xsi).*(1-eta)./4;                                      % Edge shape functions
    sNr(1:r,2) = (1-xsi).*(1+eta)./4;
    sNr(1:r,3) = (1+xsi).*(1+eta)./4;
    sNr(1:r,4) = (1+xsi).*(1-eta)./4;
    
    
elseif nen==8 && stress_flag == 0
    nf = 4;                                                                % Number of faces
    ngp=9;                                                                 % Number of area Gauss points
    nsgp = 4;                                                              % Number of surface Gauss points
    nD = 2;                                                                % Number of dimensions
    g1=sqrt(3/5);                                                          % Area Gauss point position (magnitude)
    
    %% ========= 8 node quad shape function and derivative values ========== %%
    gp = [-1 -1 -1  0  1  1  1  0 0;                                       % Area Gaus point position
        -1  0  1  1  1  0 -1 -1 0]'.*g1;                                 %
    w1 = (5/9)^2; w2 = (8/9).*(5/9); w3=(8/9)^2;                           % Area Gauss point weights
    vwp=[w1 w2 w1 w2 w1 w2 w1 w2 w3];                                      %
    r = (ngp*nD);
    xsi=gp(:,1);  eta=gp(:,2);                      	                   % Reference coordinate position
    vdNr(1:2:r,1)  = ((1 - eta).*((2.*xsi)+eta))./4;                       % Area shape function derivative values
    vdNr(1:2:r,2)  = -(1 - eta.^2)./2;
    vdNr(1:2:r,3)  = ((1 + eta).*((2.*xsi)-eta))./4;
    vdNr(1:2:r,4)  = -xsi.*(1 + eta);
    vdNr(1:2:r,5)  = ((1 + eta).*((2.*xsi)+eta))./4;
    vdNr(1:2:r,6)  =  (1 - eta.^2)./2;
    vdNr(1:2:r,7)  = ((1 - eta).*((2.*xsi)-eta))./4;
    vdNr(1:2:r,8)  = -xsi.*(1 - eta);
    vdNr(2:2:r,1)  = ((1 - xsi).*((2*eta)+xsi))./4;
    vdNr(2:2:r,2)  = -eta.*(1-xsi);
    vdNr(2:2:r,3)  =  (1-xsi).*((2*eta) - xsi)./4;
    vdNr(2:2:r,4)  =  (1 - (xsi.^2))./2;
    vdNr(2:2:r,5)  =  (1 + xsi).*((2*eta)+xsi)./4;
    vdNr(2:2:r,6)  = -eta.*(1+xsi);
    vdNr(2:2:r,7)  =  (1 + xsi).*((2*eta)-xsi)./4;
    vdNr(2:2:r,8)  = -(1 - xsi.^2)./2;
    r = r/nD;
    Nr(1:r,1) = (1 - xsi).*(1 - eta).*(-xsi-eta-1)./4;                     % Area shape function values
    Nr(1:r,2) = (1 - xsi).*(1 - eta.^2)./2;
    Nr(1:r,3) = (1 - xsi).*(1 + eta).*(-xsi+eta-1)./4;
    Nr(1:r,4) = (1 - xsi.^2).*(1 + eta)./2;
    Nr(1:r,5) = (1 + xsi).*(1 + eta).*( xsi+eta-1)./4;
    Nr(1:r,6) = (1 + xsi).*(1 - eta.^2)./2;
    Nr(1:r,7) = (1 + xsi).*(1 - eta).*(xsi-eta-1)./4;
    Nr(1:r,8) = (1 - xsi.^2).*(1 - eta)./2;
    
    g1=sqrt((3/7)+(2/7)*sqrt(6/5)); g2=sqrt((3/7)-(2/7)*sqrt(6/5));        % Surface Gauss point position (magnitude)
    gp = [-1  -1  -1 -1 -g1 -g2 g2  g1  1   1  1    1  g1  g2  -g2 -g1;    % Surface Gauss point position
        -g1 -g2 g2 g1  1   1   1   1  g1 g2  -g2 -g1 -1  -1  -1  -1]';
    w1 = (18-sqrt(30))/36; w2 = (18+sqrt(30))/36;                          % Edge Guass point weight
    swp = [w1 w2 w2 w1];
    r = nsgp*nf*nD;
    xsi = gp(:,1); eta=gp(:,2);                                            % Reference coordinate position
    sdNr(1:2:r,1)  = ((1 - eta).*((2.*xsi)+eta))./4;                       % Surface shape function derivative values
    sdNr(1:2:r,2)  = -(1 - eta.^2)./2;
    sdNr(1:2:r,3)  = ((1 + eta).*((2.*xsi)-eta))./4;
    sdNr(1:2:r,4)  = -xsi.*(1 + eta);
    sdNr(1:2:r,5)  = ((1 + eta).*((2.*xsi)+eta))./4;
    sdNr(1:2:r,6)  =  (1 - eta.^2)./2;
    sdNr(1:2:r,7)  = ((1 - eta).*((2.*xsi)-eta))./4;
    sdNr(1:2:r,8)  = -xsi.*(1 - eta);
    sdNr(2:2:r,1)  = ((1 - xsi).*((2*eta)+xsi))./4;
    sdNr(2:2:r,2)  = -eta.*(1-xsi);
    sdNr(2:2:r,3)  =  (1-xsi).*((2*eta) - xsi)./4;
    sdNr(2:2:r,4)  =  (1 - (xsi.^2))./2;
    sdNr(2:2:r,5)  =  (1 + xsi).*((2*eta)+xsi)./4;
    sdNr(2:2:r,6)  = -eta.*(1+xsi);
    sdNr(2:2:r,7)  =  (1 + xsi).*((2*eta)-xsi)./4;
    sdNr(2:2:r,8)  = -(1 - xsi.^2)./2;
    r = r/nD;
    sNr(1:r,1) = (1 - xsi).*(1 - eta).*(-xsi-eta-1)./4;                    % Surface shape function values
    sNr(1:r,2) = (1 - xsi).*(1 - eta.^2)./2;
    sNr(1:r,3) = (1 - xsi).*(1 + eta).*(-xsi+eta-1)./4;
    sNr(1:r,4) = (1 - xsi.^2).*(1 + eta)./2;
    sNr(1:r,5) = (1 + xsi).*(1 + eta).*(xsi+eta-1)./4;
    sNr(1:r,6) = (1 + xsi).*(1 - eta.^2)./2;
    sNr(1:r,7) = (1 + xsi).*(1 - eta).*(xsi-eta-1)./4;
    sNr(1:r,8) = (1 - xsi.^2).*(1 - eta)./2;
    
elseif nen==4 && stress_flag == 1
    gp = [-1 -1  1  1
        -1  1  1 -1 ]';
    ngp=size(gp,1);
    r1 = (ngp*2);
    xsi=gp(:,1) ;
    eta=gp(:,2) ;
    r=r1;
    vdNr(1:2:r,1)  = -(1 - eta)./4;
    vdNr(1:2:r,2)  = -(1 + eta)./4;
    vdNr(1:2:r,3)  =  (1 + eta)./4;
    vdNr(1:2:r,4)  =  (1 - eta)./4;
    vdNr(2:2:r,1)  = -(1 - xsi)./4;
    vdNr(2:2:r,2)  =  (1 - xsi)./4;
    vdNr(2:2:r,3)  =  (1 + xsi)./4;
    vdNr(2:2:r,4)  = -(1 + xsi)./4;
    vwp=[];swp=[]; sdNr =[]; nsgp=[];sNr=[];Nr=[];
elseif nen==8 && stress_flag == 1
    gp = [-1 -1 -1  0  1  1   1  0;
        -1  0  1  1  1  0  -1 -1]';
    ngp=size(gp,1);
    r1 = (ngp*2);
    xsi=gp(:,1) ;
    eta=gp(:,2) ;
    r=r1;
    vdNr(1:2:r,1)  = ((1 - eta).*((2.*xsi)+eta))./4;
    vdNr(1:2:r,2)  = -(1 - eta.^2)./2;
    vdNr(1:2:r,3)  = ((1 + eta).*((2.*xsi)-eta))./4;
    vdNr(1:2:r,4)  = -xsi.*(1 + eta);
    vdNr(1:2:r,5)  = ((1 + eta).*((2.*xsi)+eta))./4;
    vdNr(1:2:r,6)  =  (1 - eta.^2)./2;
    vdNr(1:2:r,7)  = ((1 - eta).*((2.*xsi)-eta))./4;
    vdNr(1:2:r,8)  = -xsi.*(1 - eta);
    
    vdNr(2:2:r,1)  = ((1 - xsi).*((2*eta)+xsi))./4;
    vdNr(2:2:r,2)  = -eta.*(1-xsi);
    vdNr(2:2:r,3)  =  (1-xsi).*((2*eta) - xsi)./4;
    vdNr(2:2:r,4)  =  (1 - (xsi.^2))./2;
    vdNr(2:2:r,5)  =  (1 + xsi).*((2*eta)+xsi)./4;
    vdNr(2:2:r,6)  = -eta.*(1+xsi);
    vdNr(2:2:r,7)  =  (1 + xsi).*((2*eta)-xsi)./4;
    vdNr(2:2:r,8)  = -(1 - xsi.^2)./2;
    vwp=[];swp=[]; sdNr =[]; nsgp=[];sNr=[];Nr=[];
end

