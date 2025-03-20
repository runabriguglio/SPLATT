function [x_vtx, y_vtx, z_vtx, edg, fcs, sld]=get_geometry_solid(model)

clear vtx x_vtx y_vtx z_vtx sld edg fcs

vtx=model.geom('geom1').getVertexCoord';
x_vtx=vtx(:,1);
y_vtx=vtx(:,2);
z_vtx=vtx(:,3);

if model.geom('geom1').getNDomains>1
	tmp=model.geom('geom1').getAdj(3,0);
	for i=2:length(tmp)
		sld(i-1,1).n=i-1;
		sld(i-1,1).nodes=double(tmp{i});
		sld(i-1,:).coo=[x_vtx(double(tmp{i})) y_vtx(double(tmp{i})) z_vtx(double(tmp{i}))];
	end
	clear tmp
	tmp=model.geom('geom1').getAdj(3,2);
	for i=2:length(tmp)
		sld(i-1,1).face=double(tmp{i});
	end
	clear tmp
	tmp=model.geom('geom1').getAdj(3,1);
	for i=2:length(tmp)
		sld(i-1,1).edge=double(tmp{i});
	end
	clear tmp
elseif model.geom('geom1').getNDomains==1
	sld(1).n=1;
	sld(1).coo=[x_vtx y_vtx z_vtx];
	tmp=model.geom('geom1').getAdj(3,0);
	sld(1).nodes=unique([tmp{1};tmp{2}]);
	clear tmp
	sld(1).edge=unique(model.geom('geom1').getAdj(3,1),'rows')';
	sld(1).face=unique(model.geom('geom1').getAdj(3,2),'rows');
end

tmp=model.geom('geom1').getAdj(1,0);
if isempty(tmp)~=1
	for i=2:length(tmp)
	        edg(i-1,1).n=i-1;
	        edg(i-1,1).nodes=double(tmp{i});
	        edg(i-1,:).coo=[x_vtx(double(tmp{i})) y_vtx(double(tmp{i})) z_vtx(double(tmp{i}))];
	end
else
	edg=[];
end
clear tmp

tmp=model.geom('geom1').getAdj(2,1);
if isempty(tmp)~=1
	for i=2:length(tmp)
	        fcs(i-1,1).n=i-1;
	        fcs(i-1,1).edges=double(tmp{i});
	end
else
	fcs=[];
end
clear tmp

tmp=model.geom('geom1').getAdj(2,0);
if isempty(tmp)~=1
	for i=2:length(tmp)
	        fcs(i-1,:).coo=[x_vtx(double(tmp{i})) y_vtx(double(tmp{i})) z_vtx(double(tmp{i}))];
		fcs(i-1,1).nodes=double(tmp{i});
	end
else
	fcs=[];
end
clear tmp

end
