function stressplotter(E,v,coord,etpl,nD,nen,d,nels)
[z_s]=stress_calc(E,v,coord,etpl,nD,nen,d,nels); 
trisurf(etpl,coord(:,1),coord(:,2),z_s,'EdgeColor','none');
view(0, 90);

