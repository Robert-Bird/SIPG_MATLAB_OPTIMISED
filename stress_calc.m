function [z]=stress_calc(E,v,coord,etpl,nD,nen,d,nels)
[~,dNr,ngp,~,~,~,~,~]=dershapefunc2D(nen,1); 

const = E/(1-(v^2));
bta = const*(1);
alh = const*(v);% Values in stiffness [De] matrix
gma = E/(2*(1+v));
De = [bta alh 0;
    alh bta 0;
    0  0  gma];
B = zeros(3,nD*nen);
p_sigma=zeros(max(etpl(:)),1);
for nel=1:nels                                                             % START of element loop
    ed=ones(nD,1)*(etpl(nel,:)-1)*nD+(1:nD).'*ones(1,nen);
    ed=reshape(ed,1,nen*nD);                                               % element degrees of freedom
    JT=dNr*coord(etpl(nel,:),:);                                           % Jacobian matrix (all Gauss points)
    temp=zeros(1,3);
    for gp=1:ngp                                                           % Start of Gauss point loop
        indx=nD*(gp-1)+(1:nD);                                             % index for current Gp information
        dNx=JT(indx,:)\dNr(indx,:);                                        % global shape func. derivatives
        B([1 3],1:2:end)=dNx;                                              % strain displacement matrix
        B([3 2],2:2:end)=dNx;
        strain=B*d(ed);                                                    % Strain
        s=De*strain;                                                       % Cauchy stress
        s1 = ((s(1)+s(2))/2)+(((s(1)-s(2))/2)^2 + s(3)^2)^0.5;
        s2 = ((s(1)+s(2))/2)-(((s(1)-s(2))/2)^2 + s(3)^2)^0.5;
        node=etpl(nel,gp);
        p_sigma(node) =  (s1-s2);
    end                                                                     
    
end
z = p_sigma;
end


