% forward kinematics of upper limb
% which has 7 dof 3 shoulder, 1 elbow, 1 forearm & 2 wrist
%

clear;

gl_rotm = [0, -1, 0; ...
    1, 0, 0; ...
    0, 0, 1];

markerNameFile = 'c:\throw2020\markername.csv';
markerName = importdata(markerNameFile);


n_rsh1 = 5;
n_rsh2 = 6;
n_rel1 = 8;
n_rel2 = 9;
n_rwr1 = 10;
n_rwr2 = 11;
n_rhd = 12;
n_lsh1 = 15;
n_lsh2 = 16;

n_clav = 19;
n_c7 = 20;
n_xip = 21;
n_xipb = 22;

staticFile = 'c:\throw2020\b\b_static1.trc';

motFile = 'c:\throw2020\b\b_con2_6_7_Trimmed-minibas_throw28.but';

if ~isfile(staticFile)
    error('There is no such a file.')
end

x_static = dlmread(staticFile, '\t', 6, 0);
x_mot = dlmread(motFile, '', 3, 0);

% determine number of frames
[nnn, nc] = size(x_mot);
[nnn2, nc2] = size(x_static);
nc2 = nc2 - 1;

% add virtual markers to marker data
% 29:n_rsh, 30:n_rel, 31:n_rwr, 32:n_lsh, 
x_static(:,nc2+1:nc2+3) = (x_static(:,2+(n_rsh1-1)*3+1:2+(n_rsh1-1)*3+3) + ...
    x_static(:,2+(n_rsh2-1)*3+1:2+(n_rsh2-1)*3+3))/2;
x_static(:,nc2+4:nc2+6) = (x_static(:,2+(n_rel1-1)*3+1:2+(n_rel1-1)*3+3) + ...
    x_static(:,2+(n_rel2-1)*3+1:2+(n_rel2-1)*3+3))/2;
x_static(:,nc2+7:nc2+9) = (x_static(:,2+(n_rwr1-1)*3+1:2+(n_rwr1-1)*3+3) + ...
    x_static(:,2+(n_rwr2-1)*3+1:2+(n_rwr2-1)*3+3))/2;
x_static(:,nc2+10:nc2+12) = (x_static(:,2+(n_lsh1-1)*3+1:2+(n_lsh1-1)*3+3) + ...
    x_static(:,2+(n_lsh2-1)*3+1:2+(n_lsh2-1)*3+3))/2;

x_mot(:,nc+1:nc+3) = (x_mot(:,1+(n_rsh1-1)*3+1:1+(n_rsh1-1)*3+3) + ...
    x_mot(:,1+(n_rsh2-1)*3+1:1+(n_rsh2-1)*3+3))/2;
x_mot(:,nc+4:nc+6) = (x_mot(:,1+(n_rel1-1)*3+1:1+(n_rel1-1)*3+3) + ...
    x_mot(:,1+(n_rel2-1)*3+1:1+(n_rel2-1)*3+3))/2;
x_mot(:,nc+7:nc+9) = (x_mot(:,1+(n_rwr1-1)*3+1:1+(n_rwr1-1)*3+3) + ...
    x_mot(:,1+(n_rwr2-1)*3+1:1+(n_rwr2-1)*3+3))/2;
x_mot(:,nc+10:nc+12) = (x_mot(:,1+(n_lsh1-1)*3+1:1+(n_lsh1-1)*3+3) + ...
    x_mot(:,1+(n_lsh2-1)*3+1:1+(n_lsh2-1)*3+3))/2;

% define body segment using structure
% 1: upper trunk 2: upper arm 3: forearm 4: hand
seg(1).nMarkers = 4;
seg(1).idMarker = [19 20 21 22];
seg(2).nMarkers = 5;
seg(2).idMarker = [5 6 7 8 9];
seg(3).nMarkers = 3;
seg(3).idMarker = [30 10 11];
seg(4).nMarkers = 3;
seg(4).idMarker = [10 11 12];

nSeg = size(seg, 2);
% loop for different frames (time instants)
nnn = 1;
for i = 1:nnn
    % dertermine the weighting with segmental residual errors
    % loop for different segments
    for iSeg = 1:nSeg
        for iMarker=1:seg(iSeg).nMarkers
            np = seg(iSeg).idMarker(iMarker);
            x(iMarker,:) = x_static(100, (2+np*3+1):(2+np*3+3));
            y(iMarker,:) = x_mot(i, 1+np*3+1:1+np*3+3);
        end
        [R, d] = RotationMatrix_SVD(x, y);

        % residual of segmental optimization
        for iMarker=1:seg(1).nMarkers
            err = y(iMarker, :)-(R*x(iMarker,:)')'-d';
            rsd(iMarker) = norm(err); 
        end
        wt(iSeg)=1/mean(rsd);
    end
end

[seg2, q_i] = upperlimb_model_parameter(x_static);

% determine transformation matrix for each joint rotation
% upper trunk (yxz)
axang = [0 1 0 q_i(1)]; Rut_y = axang2rotm(axang); 
axang = [1 0 0 q_i(2)]; Rut_x = axang2rotm(axang); 
axang = [0 0 1 q_i(3)]; Rut_z = axang2rotm(axang); 
Rut = gl_rotm*Rut_y*Rut_x*Rut_z;
tform_ut = rotm2tform(Rut);

% shoulder (xyz)
axang = [1 0 0 q_i(4)]; Rsh_x = axang2rotm(axang); 
axang = [0 1 0 q_i(5)]; Rsh_y = axang2rotm(axang); 
axang = [0 0 1 q_i(6)]; Rsh_z = axang2rotm(axang); 
Rsh = Rsh_x*Rsh_y*Rsh_z;
tform_sh = rotm2tform(Rsh);

% elbow (yxz)
axang = [0 1 0 q_i(7)]; Rel_y = axang2rotm(axang); 
axang = [1 0 0 q_i(8)]; Rel_x = axang2rotm(axang); 
axang = [0 0 1 q_i(9)]; Rel_z = axang2rotm(axang); 
Rel = Rel_y*Rel_x*Rel_z;
tform_el = rotm2tform(Rel);

% wrist (yxz)
axang = [0 1 0 q_i(10)]; Rwr_y = axang2rotm(axang); 
axang = [1 0 0 q_i(11)]; Rwr_x = axang2rotm(axang); 
axang = [0 0 1 q_i(12)]; Rwr_z = axang2rotm(axang); 
Rwr = Rwr_y*Rwr_x*Rwr_z;
tform_wr = rotm2tform(Rwr);

% homogeneous transformation
% upper trunk
p1 = seg2(1).origin;
tform_ut_tr = trvec2tform(p1/1000);

% shoulder
p2 = seg2(1).length_vector_local;
tform_sh_tr = trvec2tform(p2'/1000);

% elbow
p3 = seg2(2).length_vector_local;
tform_el_tr = trvec2tform(p3'/1000);

% wrist
p4 = seg2(3).length_vector_local;
tform_wr_tr = trvec2tform(p4'/1000);

% hand
p5 = seg2(4).length_vector_local;
tform_ha_tr = trvec2tform(p5'/1000);

% forward kinematics of the upper limb with the hand being end-point
% tform_ut_tr*tform_ut*tform_sh_tr*tform_sh*tform_el_tr*tform_el*tform_wr_tr*tform_wr*tform_ha_tr;

% loop for motion data frames
for iFrame = 1:nnn
    
    
end


% define body markers of static trial
xs_clav = x_static(:, 2+(n_clav-1)*3+1:2+(n_clav-1)*3+3);
xs_c7 = x_static(:, 2+(n_c7-1)*3+1:2+(n_c7-1)*3+3);
xs_xip = x_static(:, 2+(n_xip-1)*3+1:2+(n_xip-1)*3+3);
xs_xipb = x_static(:, 2+(n_xipb-1)*3+1:2+(n_xipb-1)*3+3);

xs_ut = (xs_xip + xs_xipb)/2;

xs_rsh1 = x_static(:, 2+(n_rsh1-1)*3+1:2+(n_rsh1-1)*3+3);
xs_rsh2 = x_static(:, 2+(n_rsh2-1)*3+1:2+(n_rsh2-1)*3+3);
xs_rel1 = x_static(:, 2+(n_rel1-1)*3+1:2+(n_rel1-1)*3+3);
xs_rel2 = x_static(:, 2+(n_rel2-1)*3+1:2+(n_rel2-1)*3+3);
xs_rwr1 = x_static(:, 2+(n_rwr1-1)*3+1:2+(n_rwr1-1)*3+3);
xs_rwr2 = x_static(:, 2+(n_rwr2-1)*3+1:2+(n_rwr2-1)*3+3);

xs_rhd = x_static(:, 2+(n_rhd-1)*3+1:2+(n_rhd-1)*3+3);

xs_rsh = (xs_rsh1 + xs_rsh2)/2;
xs_rel = (xs_rel1 + xs_rel2)/2;
xs_rwr = (xs_rwr1 + xs_rwr2)/2;

