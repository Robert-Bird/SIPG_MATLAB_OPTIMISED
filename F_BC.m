function [f,BC,d,c]=F_BC(nelsC,nelsR,etpl,ed,coord,pr,nen,nD)
f=zeros(nelsR*nelsC*nen*nD,1); th=(pi/2)/(nelsC*2);d=f; c=d;
if nen == 4
    node(1)=3; node(2)=4;
    el_BC_1 = 1:nelsC:(nelsC*nelsR);                                           % Element numbers to horiztontal face roller BC
    el_BC_2 = nelsC:nelsC:(nelsC*nelsR);                                       % Element numbers to vertical face roller BC
    BC=reshape([ed(el_BC_1,[2,8]);ed(el_BC_2,[3,5])],1,4*nelsR);               % Variable of BC dof.
   
elseif nen ==8
    node(1)=5; node(2)=7;
    el_BC_1 = 1:nelsC:(nelsC*nelsR);                                           % Element numbers to horiztontal face roller BC
    el_BC_2 = nelsC:nelsC:(nelsC*nelsR);                                       % Element numbers to vertical face roller BC
    BC=reshape([ed(el_BC_1,[2,14,16]);ed(el_BC_2,[5,7,9])],1,6*nelsR);
    
end
if nen == 4
    for i=1:nelsC/2
        el=nelsC*(nelsR-1)+i;
        node1 = etpl(el,node(1));
        node2 = etpl(el,node(2));
        h=coord(node1,2)-coord(node2,2);
        dofx = [ed(el,5),ed(el,7)];
        dofxy = [ed(el,6),ed(el,8)];
        t=th;
        r=10./cos(t);
        f(dofx(2))=f(dofx(2))+(h/2)*pr*(1-(1.^2./r.^2).*((3.*cos(2.*t)./2) + cos(4.*t)) + (3.*(1.^4./r.^4).*cos(4.*t)./2));
        f(dofx(1))=f(dofx(1))+(h/2)*pr*(1-(1.^2./r.^2).*((3.*cos(2.*t)./2) + cos(4.*t)) + (3.*(1.^4./r.^4).*cos(4.*t)./2));
        f(dofxy(2))=f(dofxy(2))+(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        f(dofxy(1))=f(dofxy(1))+(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        th=th+(pi/2)/(nelsC); 
    end
   
    for i=1:nelsC/2
        el=nelsC*(nelsR-1)+(nelsC/2)+i;
        node1 = etpl(el,node(1));
        node2 = etpl(el,node(2));
        h=coord(node1,1)-coord(node2,1);
        dofy = [ed(el,5),ed(el,7)];
        dofxy = [ed(el,6),ed(el,8)];
        t=th;
        r=10./cos(t);
        t=t+pi/4;
        f(dofy(2))=f(dofy(2))+(h/2)*pr*(-(1.^2./r.^2).*((cos(2.*th)./2) - cos(4.*th)) - (3.*(1.^4./r.^4).*cos(4.*th)./2));
        f(dofy(1))=f(dofy(1))+(h/2)*pr*(-(1.^2./r.^2).*((cos(2.*th)./2) - cos(4.*th)) - (3.*(1.^4./r.^4).*cos(4.*th)./2));
        f(dofxy(2))=f(dofxy(2))+(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        f(dofxy(1))=f(dofxy(1))+(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        th=th+(pi/2)/(nelsC); 
    end
end
if nen == 8
    for i=1:nelsC/2
        el=nelsC*(nelsR-1)+i;
        node1 = etpl(el,node(1));
        node2 = etpl(el,node(2));
        h=coord(node1,2)-coord(node2,2);
        dofx = [ed(el,13),ed(el,11),ed(el,9)];
        dofxy= [ed(el,14),ed(el,12),ed(el,10)];
        w=[1,1];
        gp=[-sqrt(3)/3, sqrt(3)/3];
        for i2=1:2
            xsi=gp(i2);
            t=atan((tan(th))+xsi*(h/20));
            r=10./cos(t);
            n1=(-1/2) *xsi*(1-xsi);
            n2=(1-xsi)*(1+xsi);
            n3=(1/2)  *xsi*(1+xsi);
            f(dofx(3))=f(dofx(3))+n3*w(i2)*(h/2)*pr*(1-(1.^2./r.^2).*((3.*cos(2.*t)./2) + cos(4.*t)) + (3.*(1.^4./r.^4).*cos(4.*t)./2));
            f(dofx(2))=f(dofx(2))+n2*w(i2)*(h/2)*pr*(1-(1.^2./r.^2).*((3.*cos(2.*t)./2) + cos(4.*t)) + (3.*(1.^4./r.^4).*cos(4.*t)./2));
            f(dofx(1))=f(dofx(1))+n1*w(i2)*(h/2)*pr*(1-(1.^2./r.^2).*((3.*cos(2.*t)./2) + cos(4.*t)) + (3.*(1.^4./r.^4).*cos(4.*t)./2));
            f(dofxy(3))=f(dofxy(3))+n3*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
            f(dofxy(2))=f(dofxy(2))+n2*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
            f(dofxy(1))=f(dofxy(1))+n1*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        end
        th=th+(pi/2)/(nelsC);
    end
    
    
    
    th=(pi/2)/(nelsC);
    for i=1:nelsC/2
        el=nelsC*(nelsR-1)+(nelsC/2)+i;
        node1 = etpl(el,node(1));
        node2 = etpl(el,node(2));
        h=coord(node1,1)-coord(node2,1);
        dofy = [ed(el,13),ed(el,11),ed(el,9)];
        dofxy= [ed(el,14),ed(el,12),ed(el,10)];
        w=[1,1];
        gp=[-sqrt(3)/3, sqrt(3)/3];
        for i2=1:2
            xsi=gp(i2);
            t=atan((tan(th))+xsi*(h/20));
            r=10./cos(t);
            t=t+pi/4;
            n1=(-1/2) *xsi*(1-xsi);
            n2=(1-xsi)*(1+xsi);
            n3=(1/2)  *xsi*(1+xsi);
            f(dofy(3))=f(dofy(3))+n3*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*((cos(2.*th)./2) - cos(4.*th)) - (3.*(1.^4./r.^4).*cos(4.*th)./2));
            f(dofy(2))=f(dofy(2))+n2*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*((cos(2.*th)./2) - cos(4.*th)) - (3.*(1.^4./r.^4).*cos(4.*th)./2));
            f(dofy(1))=f(dofy(1))+n1*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*((cos(2.*th)./2) - cos(4.*th)) - (3.*(1.^4./r.^4).*cos(4.*th)./2));
            f(dofxy(3))=f(dofxy(3))+n3*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
            f(dofxy(2))=f(dofxy(2))+n2*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
            f(dofxy(1))=f(dofxy(1))+n1*w(i2)*(h/2)*pr*(-(1.^2./r.^2).*(sin(2.*t)./2 + sin(4.*t)) + (3.*(1.^4./r.^4).*sin(4.*t)./2));
        end
        th=th+(pi/2)/(nelsC);
    end
end
    
end

