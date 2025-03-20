function hexagonal_shell() %out_path)
% Performs the FEA of a hexagonal shell from the parameters read in
% 'mech_parameters.txt' and the actuator coordinates in 'act_coords.txt'
%
% Therefore, the current folder must contain the SimInputFiles folder
% with the following files before running the simulation:
%
%       mech_parameters.txt: contains a vector of mechanical parameters of
%                            the shell in SI units ordered as follows:
%       [RoC, thk, side_len, E, rho, nu]
%   
%           - RoC:      the radius of curvature of the shell [m]
%           - thk:      the shell's thickness [m]
%           - side_len: the hexagon side length [m]
%           - E:        the material's Young modulus [Pa]
%           - rho:      the material's density [kg/m^3]
%           - nu:       the material's Poisson ratio [-]
%
%       act_coords.txt: contains the vector of actuator (x,y) coordinates
%                       of size [Nacts,2]
%
%
% The stiffness matrix and influence functions are then saved as .txt files
% in the SimOutputFiles folder
%
%

% Define paths
input_path = 'SimInputFiles/';
output_path = 'SimOutputFiles/';
comsol_filepath = 'c:/Program Files/COMSOL/COMSOL62/Multiphysics/'; % CHANGE THIS TO THE COMSOL PATH ON YOUR MACHINE!

% Kill any running mph process
system('taskkill /IM comsolmphserver.exe /F');
pause(2)

% Start server
system(sprintf('%s/bin/win64/comsolmphserver &',comsol_filepath));

% Add path
addpath(sprintf('%s/mli',comsol_filepath));
port = 2036;
mphstart(port);

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model'); 

% Data
mech_data = load(sprintf('%s/mech_parameters.txt',input_path));
RoC = mech_data(1); % [m] radius of curvature
thk = mech_data(2); % [m] shell's thickness
hex_side = mech_data(3); % [m] hexagon side
Ez = mech_data(4); % [Pa] Zerodur's Young modulus
rhoz = mech_data(5); % [kg/m^3] Zerodur's density
nuz = mech_data(6); % Zerodur's Poisson ratio

act_coords = load(sprintf('%s/act_coords.txt',input_path));
Nacts = length(act_coords);

displ = -1e-6; % [m]

% Set model parameters
model.param.set('RoC', sprintf('%1.1f[m]',RoC));
model.param.set('hex_side_len', sprintf('%1.1f[m]',hex_side));
model.param.set('E_z', sprintf('%1.1f[Pa]',Ez));
model.param.set('nu_z',sprintf('%1.2f',nuz));
model.param.set('rho_z', sprintf('%1.1f[kg/m^3]',rhoz));
model.param.set('thk',sprintf('%1.4f[m]',thk));
model.param.set('displ', sprintf('%1.2f[um]',displ*1e+6));
model.param.set('k_p', '100 [N/m]');

% Component geometry
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
model.component('comp1').geom('geom1').feature('sph1').set('rot', 45);
model.component('comp1').geom('geom1').feature('sph1').set('type', 'surface');
model.component('comp1').geom('geom1').feature('sph1').set('axis', [45 45]);
model.component('comp1').geom('geom1').feature('sph1').set('r', 'RoC');
model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl1').set('type', 'surface');
model.component('comp1').geom('geom1').feature('cyl1').set('r', '2*hex_side_len');
model.component('comp1').geom('geom1').feature('cyl1').set('h', '1+RoC');
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'cyl1' 'sph1'});
model.component('comp1').geom('geom1').create('del1', 'Delete');
model.component('comp1').geom('geom1').feature('del1').selection('input').set('uni1(1)', [1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17]);
model.component('comp1').geom('geom1').create('ps1', 'ParametricSurface');
model.component('comp1').geom('geom1').feature('ps1').set('parmin1', '-hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps1').set('parmax1', 'hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps1').set('parmin2', '-hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps1').set('parmax2', 'hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps1').set('coord', {'s1' 'hex_side_len*sin(pi/3)' 'RoC+s2'});
model.component('comp1').geom('geom1').create('ps3', 'ParametricSurface');
model.component('comp1').geom('geom1').feature('ps3').set('parmax1', 'hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps3').set('parmin2', '-hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps3').set('parmax2', 'hex_side_len*2');
model.component('comp1').geom('geom1').feature('ps3').set('coord', {'s1' 'hex_side_len*sin(pi/3) - (s1-hex_side_len/2)*tan(pi/3)' 'RoC+s2'});
model.component('comp1').geom('geom1').create('mir1', 'Mirror');
model.component('comp1').geom('geom1').feature('mir1').set('keep', true);
model.component('comp1').geom('geom1').feature('mir1').set('axis', [1 0 0]);
model.component('comp1').geom('geom1').feature('mir1').selection('input').set({'ps3'});
model.component('comp1').geom('geom1').create('mir2', 'Mirror');
model.component('comp1').geom('geom1').feature('mir2').set('keep', true);
model.component('comp1').geom('geom1').feature('mir2').set('axis', [0 1 0]);
model.component('comp1').geom('geom1').feature('mir2').selection('input').set({'mir1' 'ps1' 'ps3'});
model.component('comp1').geom('geom1').create('uni2', 'Union');
model.component('comp1').geom('geom1').feature('uni2').selection('input').set({'del1' 'mir1' 'mir2' 'ps1' 'ps3'});
model.component('comp1').geom('geom1').create('del2', 'Delete');
model.component('comp1').geom('geom1').feature('del2').selection('input').set('uni2(1)', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59]);
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickplane', 'zx');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').create('parf2', 'PartitionFaces');
model.component('comp1').geom('geom1').feature('parf2').set('partitionwith', 'workplane');
model.component('comp1').geom('geom1').feature('parf2').selection('face').set('del2(1)', 1);
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').set('planetype', 'normalvector');
model.component('comp1').geom('geom1').feature('wp2').set('normalvector', {'sin(pi/3)' 'cos(pi/3)' '0'});
model.component('comp1').geom('geom1').feature('wp2').set('unite', true);
model.component('comp1').geom('geom1').create('parf3', 'PartitionFaces');
model.component('comp1').geom('geom1').feature('parf3').set('partitionwith', 'workplane');
model.component('comp1').geom('geom1').feature('parf3').selection('face').set('parf2(1)', [1 2]);
model.component('comp1').geom('geom1').create('wp3', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp3').set('planetype', 'normalvector');
model.component('comp1').geom('geom1').feature('wp3').set('normalvector', {'-sin(pi/3)' 'cos(pi/3)' '0'});
model.component('comp1').geom('geom1').feature('wp3').set('unite', true);
model.component('comp1').geom('geom1').create('parf4', 'PartitionFaces');
model.component('comp1').geom('geom1').feature('parf4').set('partitionwith', 'workplane');
model.component('comp1').geom('geom1').feature('parf4').selection('face').set('parf3(1)', [1 3]);

% Add points to act coords locations
x = act_coords(:,1);
y = act_coords(:,2);
z = sqrt(RoC^2 - x.^2 - y.^2);
for i = 1:Nacts
    model.component('comp1').geom('geom1').create(sprintf('pt%d',i), 'Point');
    model.component('comp1').geom('geom1').feature(sprintf('pt%d',i)).set('p', {sprintf('%1.6f',x(i)) sprintf('%1.6f',y(i)) sprintf('%1.6f',z(i))});
end

model.component('comp1').geom('geom1').run;

% Create mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').run;

% Coordinate systems
model.component('comp1').coordSystem.create('sys2', 'Cylindrical');
model.component('comp1').coordSystem.create('sys3', 'Spherical');

% Shell and material physics
model.component('comp1').physics.create('shell', 'Shell', 'geom1');
model.component('comp1').physics('shell').feature('emm1').set('E_mat', 'userdef');
model.component('comp1').physics('shell').feature('emm1').set('E', 'E_z');
model.component('comp1').physics('shell').feature('emm1').set('nu_mat', 'userdef');
model.component('comp1').physics('shell').feature('emm1').set('nu', 'nu_z');
model.component('comp1').physics('shell').feature('emm1').set('rho_mat', 'userdef');
model.component('comp1').physics('shell').feature('emm1').set('rho', 'rho_z');
model.component('comp1').physics('shell').feature('to1').set('d', 'thk');

% Find actuator indices
[x_vtx, y_vtx, ~, edges, ~] = get_geometry_solid(model);

dist = @(x,y) sqrt(x.^2+y.^2);
act_ids = zeros(Nacts,1);

for i = 1:Nacts
    xi = act_coords(i,1);
    yi = act_coords(i,2);

    d_vtx = dist(x_vtx-xi,y_vtx-yi);
    [~,vtx_id] = min(d_vtx);
    act_ids(i) = vtx_id;
end

% find origin
[~,origin] = min(x_vtx.^2 + y_vtx.^2);

% Find edge boundary ids
clear boundary_edge_ids
k = 0;
for ik = 1:length(edges)
    xi = edges(ik).coo(:,1);
    yi = edges(ik).coo(:,2);

    if abs(xi) <= hex_side*0.5

        if abs(yi) >= (hex_side * sin(pi/3))*0.99
            k = k + 1;
            boundary_edge_ids(k) = ik;
%             xi,yi
        end

    else

        if abs(abs(yi) - (sin(pi/3) - (abs(xi)-0.5)*tan(pi/3))) <= 1e-13
            k = k + 1;
            boundary_edge_ids(k) = ik;
%             xi,yi
        end

    end

end


% No rotation on boundary
model.component('comp1').physics('shell').create('disp2', 'Displacement1', 1);
model.component('comp1').physics('shell').feature('disp2').selection.set(boundary_edge_ids);
model.component('comp1').physics('shell').feature('disp2').set('Direction', {'free'; 'prescribed'; 'free'});
model.component('comp1').physics('shell').feature('disp2').set('coordinateSystem', 'sys2');

% Fixed actuators in radial direction
model.component('comp1').physics('shell').create('disp5', 'Displacement0', 0);
model.component('comp1').physics('shell').feature('disp5').selection.set(act_ids);
model.component('comp1').physics('shell').feature('disp5').set('Direction', {'prescribed'; 'free'; 'free'});
model.component('comp1').physics('shell').feature('disp5').set('coordinateSystem', 'sys3');

% Moving actuator
model.component('comp1').physics('shell').create('disp6', 'Displacement0', 0);
model.component('comp1').physics('shell').feature('disp6').set('Direction', {'prescribed'; 'free'; 'free'});
model.component('comp1').physics('shell').feature('disp6').set('U0', {'displ'; '0'; '0'});
model.component('comp1').physics('shell').feature('disp6').set('coordinateSystem', 'sys3');

% Spring in x,y
model.component('comp1').physics('shell').create('spf1', 'SpringFoundation0', 0);
model.component('comp1').physics('shell').feature('spf1').selection.set(origin);
model.component('comp1').physics('shell').feature('spf1').set('SpringType', 'FDef');
model.component('comp1').physics('shell').feature('spf1').set('FDef', {'k_p*shell.uspring1_spf1'; 'k_p*shell.uspring2_spf1'; '0[N/m]*shell.uspring3_spf1'});

% Cycle on actuators
stiff_mat = zeros(Nacts,Nacts);
z_displacement = 'w';

if RoC == 0 || ~isfinite(RoC) || isnan(RoC)
    normal_force = 'shell.RFz';
else
    normal_force = 'shell.RFx*x/RoC + shell.RFy*y/RoC + shell.RFz*z/RoC';
end

for k = 1:Nacts

    fprintf('\nActuator %d ...\n',k);

    study_name = sprintf('std%d',k);
    solution_name = sprintf('dset%d',k);
    dataset_name = sprintf('dataset%d',k);

    model.component('comp1').physics('shell').feature('disp6').selection.set(act_ids(k));

    % Define study
    model.study.create(study_name);
    model.study(study_name).create('stat', 'Stationary');
    
    % Create solution
    model.sol.create(solution_name);
    model.sol(solution_name).study(study_name);
    model.sol(solution_name).attach(study_name);
    model.sol(solution_name).create('st1', 'StudyStep');
    model.sol(solution_name).create('v1', 'Variables');
    model.sol(solution_name).create('s1', 'Stationary');
    model.sol(solution_name).feature('s1').create('fc1', 'FullyCoupled');
    model.sol(solution_name).feature('s1').feature.remove('fcDef');
    
    % Define solution
    model.sol(solution_name).attach(study_name);
    model.sol(solution_name).feature('v1').feature('comp1_ar').set('scalemethod', 'manual');
    model.sol(solution_name).feature('v1').feature('comp1_ar').set('scaleval', 0.01);
    model.sol(solution_name).feature('s1').feature('aDef').set('cachepattern', false); %true);
    model.sol(solution_name).runAll;

    model.result.dataset.create(dataset_name, 'Solution');
    model.result.dataset(dataset_name).set('solution', solution_name);

    % Post-processing
    for i = 1:Nacts
        k_aux = mpheval(model, normal_force, 'edim', 0, 'selection', act_ids(i), 'dataset', dataset_name);
        stiff_mat(k,i) = k_aux.d1/displ;
    end

    % Iffs
    iffstruct = mpheval(model, z_displacement, 'edim', 2, 'dataset', dataset_name); 

    iffs_vec(k,1,:) = iffstruct.p(1,:);
    iffs_vec(k,2,:) = iffstruct.p(2,:);
    iffs_vec(k,3,:) = iffstruct.d1/displ;
end

% post-processing
N = 256;
Nacts = size(iffs_vec,1);
Npts = size(iffs_vec,3);
mask = zeros(N,N);
Zdata = zeros(Nacts,N,N);

for i = 1:Nacts
    x = reshape(iffs_vec(i,1,:),[Npts,1]);
    y = reshape(iffs_vec(i,2,:),[Npts,1]);
    z = reshape(iffs_vec(i,3,:),[Npts,1]);

    [X,Y] = meshgrid(linspace(min(x),max(x),N), linspace(min(y),max(y),N));
    Z = griddata(x,y,z,X,Y);

    mask = max(mask,isnan(Z));

    Zdata(i,:,:) = Z;
end

iffs = zeros(sum(1-mask,'all'),Nacts);
for i = 1:Nacts
    z = reshape(Zdata(i,:,:),[N,N]);
    iffs(:,i) = z(mask==0);
end

mesh(:,1) = X(mask==0);
mesh(:,2) = Y(mask==0);

writematrix(stiff_mat,sprintf('%s/stiffness_matrix.txt',output_path),'Delimiter',' ');
writematrix(iffs,sprintf('%s/iffs.txt',output_path),'Delimiter',' ');
writematrix(mesh,sprintf('%s/mesh.txt',output_path),'Delimiter',' ');

% Disconnect
ModelUtil.clear;
ModelUtil.disconnect;

% Kill the terminal
system('taskkill /IM cmd.exe /F');

end
