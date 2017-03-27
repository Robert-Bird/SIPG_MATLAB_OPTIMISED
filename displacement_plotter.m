function displacement_plotter(d,etpl,coord,nen,nelsC,nelsR)

x = zeros(nelsR*nelsC,nen); y=x;
dx = d(1:2:end);dy = d(2:2:end);  % Producing new positions of mesh dof
for i=1:(nen+1)
    if i > nen
        x(:,i) = coord(etpl(:,1),1)+dx(etpl(:,1));  y(:,i) = coord(etpl(:,1),2)+dy(etpl(:,1));
    else
        x(:,i) = coord(etpl(:,i),1)+dx(etpl(:,i));  y(:,i) = coord(etpl(:,i),2)+dy(etpl(:,i));
    end
end

figure(2); hold on
for nel=1:nelsC*nelsR
    plot(coord(etpl(nel,:),1),coord(etpl(nel,:),2),'b-'); % Plotting original element map
    plot(x(nel,:),y(nel,:),'-r');                         % Plotting resultant element map
end