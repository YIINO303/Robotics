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

% define the origin of upper trunk in the motion file
x_clav = x_mot(:, 1+(n_clav-1)*3+1:1+(n_clav-1)*3+3);
x_7 = x_mot(:, 1+(n_c7-1)*3+1:1+(n_c7-1)*3+3);
x_xip = x_mot(:, 1+(n_xip-1)*3+1:1+(n_xip-1)*3+3);
x_xipb = x_mot(:, 1+(n_xipb-1)*3+1:1+(n_xipb-1)*3+3);
x_xipm = (x_xip+x_xipb)/2;
x_rhd = x_mot(:, 1+(n_rhd-1)*3+1:1+(n_rhd-1)*3+3);

x_rsh = x_mot(:, 1+(29-1)*3+1:1+(29-1)*3+3);
x_rel = x_mot(:, 1+(30-1)*3+1:1+(30-1)*3+3);
x_rwr = x_mot(:, 1+(31-1)*3+1:1+(31-1)*3+3);

p_ua = x_rsh(1,:) - x_rel(1,:);
p_ua = p_ua/norm(p_ua);
p_fa = x_rel(1,:) - x_rwr(1,:);
p_fa = p_fa/norm(p_fa);

% determine segment length
for i = 1:nnn
    length_ua(i) = norm(x_rsh(i,:)-x_rel(i,:))/1000;
    length_fa(i) = norm(x_rel(i,:)-x_rwr(i,:))/1000;
    length_ha(i) = norm(x_rwr(i,:)-x_rhd(i,:))/1000;
end

% define body segment using structure
% 1: upper trunk 2: upper arm 3: forearm 4: hand
seg(1).nMarkers = 4;
seg(1).idMarker = [19 20 21 22];
seg(2).nMarkers = 5;
seg(2).idMarker = [5 6 7 8 9];
seg(3).nMarkers = 4;
seg(3).idMarker = [30 31 10 11];
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
seg2(1).tform_trnsl = tform_sh_tr;

% elbow
p3 = seg2(2).length_vector_local;
tform_el_tr = trvec2tform(p3'/1000);
seg2(2).tform_trnsl = tform_el_tr;

% wrist
p4 = seg2(3).length_vector_local;
tform_wr_tr = trvec2tform(p4'/1000);
seg2(3).tform_trnsl = tform_wr_tr;

% hand
p5 = seg2(4).length_vector_local;
tform_ha_tr = trvec2tform(p5'/1000);
seg2(4).tform_trnsl = tform_ha_tr;

% forward kinematics of the upper limb with the hand being end-point
% tform_ut_tr*tform_ut*tform_sh_tr*tform_sh*tform_el_tr*tform_el*tform_wr_tr*tform_wr*tform_ha_tr;

% determine joint angles for motion file
q = upperlimb_JointAngle(x_mot);

% loop for different frames & 
% estimate each marker coordinates using relationship in static trial
for iFrame = 1:nnn
    
    % determine transformation matrix for each joint rotation
    % upper trunk (yxz)
    axang = [0 1 0 q(iFrame, 1)]; Rut_y = axang2rotm(axang); 
    axang = [1 0 0 q(iFrame, 2)]; Rut_x = axang2rotm(axang); 
    axang = [0 0 1 q(iFrame, 3)]; Rut_z = axang2rotm(axang); 
    Rut = gl_rotm*Rut_y*Rut_x*Rut_z;
    tform_ut = rotm2tform(Rut);
    seg_mot(1).tform_rot = tform_ut;

    % shoulder (xyz)
    axang = [1 0 0 q(iFrame, 4)]; Rsh_x = axang2rotm(axang); 
    axang = [0 1 0 q(iFrame, 5)]; Rsh_y = axang2rotm(axang); 
    axang = [0 0 1 q(iFrame, 6)]; Rsh_z = axang2rotm(axang); 
    Rsh = Rsh_x*Rsh_y*Rsh_z;
    tform_sh = rotm2tform(Rsh);
    seg_mot(2).tform_rot = tform_sh;

    % elbow (yxz)
    axang = [0 1 0 q(iFrame, 7)]; Rel_y = axang2rotm(axang); 
    axang = [1 0 0 q(iFrame, 8)]; Rel_x = axang2rotm(axang); 
    axang = [0 0 1 q(iFrame, 9)]; Rel_z = axang2rotm(axang); 
    Rel = Rel_y*Rel_x*Rel_z;
    tform_el = rotm2tform(Rel);
    seg_mot(3).tform_rot = tform_el;

    % wrist (yxz)
    axang = [0 1 0 q(iFrame, 10)]; Rwr_y = axang2rotm(axang); 
    axang = [1 0 0 q(iFrame, 11)]; Rwr_x = axang2rotm(axang); 
    axang = [0 0 1 q(iFrame, 12)]; Rwr_z = axang2rotm(axang); 
    Rwr = Rwr_y*Rwr_x*Rwr_z;
    tform_wr = rotm2tform(Rwr);
    seg_mot(4).tform_rot = tform_wr;
    
    % initialize tform martix
    tform_ut_tr = trvec2tform(x_xipm(iFrame,:)/1000);

    % set tform_sh_tr from motion file martix
    p_2_1 = x_rsh(iFrame, :) - x_xipm(iFrame,:);
    tform_sh_tr = trvec2tform((Rut'*(p_2_1/1000)')');

    % loop for different segments
    for iSeg = 1:nSeg
        % loop for different markers
        if iSeg == 1
            seg_mot(iSeg).tform = tform_ut_tr*tform_ut;
        elseif iSeg == 2
            seg_mot(iSeg).tform = tform_ut_tr*tform_ut* ...
                tform_sh_tr*tform_sh;
%             seg_mot(iSeg).tform = seg_mot(iSeg-1).tform* ...
%                 seg2(iSeg-1).tform_trnsl*seg_mot(iSeg).tform_rot;
        elseif iSeg == 3
            seg_mot(iSeg).tform = tform_ut_tr*tform_ut* ...
                tform_sh_tr*tform_sh* ...
                tform_el_tr*tform_el;
        elseif iSeg == 4
            seg_mot(iSeg).tform = tform_ut_tr*tform_ut* ...
                tform_sh_tr*tform_sh* ...
                tform_el_tr*tform_el* ...
                tform_wr_tr*tform_wr;
        else
        end
        
        for iMarker = 1:seg2(iSeg).nMarkers
            % local coordinates of iMarker attached to iSeg
            p_local = seg2(iSeg).marker(iMarker).p_local;
            % unit conversion from mm to m
            p_local = p_local/1000;
            tform_marker = trvec2tform(p_local');
            tform_marker_global = seg_mot(iSeg).tform*tform_marker;
            seg_mot(iSeg).marker(iMarker).x = tform_marker_global(1:3,4);
            
            n_iMarker = seg2(iSeg).idMarker(iMarker);
            x_iMarker = x_mot(iFrame,1+(n_iMarker-1)*3+1:1+(n_iMarker-1)*3+3)/1000;
            seg_mot(iSeg).marker(iMarker).x_measured = x_iMarker;
        end
        
        
        
    end


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

