function {q] = upperlimb_JointAngle(seg, x_mot)
%UPPERLIMB_JOINTANGLE 
%   input      
%       x_mot           marker coordinates for a motion trial
%   output
%       q[nFrmae, 12]   joint angle vector

n_clav = 19; n_c7 = 20; n_xip = 21; n_xipb = 22;
n_rsh1 = 5; n_rsh2 = 6; n_rel1 = 8; n_rel2 = 9;
n_rwr1 = 10; n_rwr2 = 11; n_rhd = 12;
n_lsh1 = 15; n_lsh2 = 16;
n_rsh = 29; n_rel = 30; n_rwr = 31; n_lsh = 32;

x_clav = x_mot(:, 1+(n_clav-1)*3+1:1+(n_clav-1)*3+3);
x_c7 = x_mot(:, 1+(n_c7-1)*3+1:1+(n_c7-1)*3+3);
x_xip = x_mot(:, 1+(n_xip-1)*3+1:1+(n_xip-1)*3+3);
x_xipb = x_mot(:, 1+(n_xipb-1)*3+1:1+(n_xipb-1)*3+3);

x_rsh = x_mot(:, 1+(n_rsh-1)*3+1:1+(n_rsh-1)*3+3);
x_rel = x_mot(:, 1+(n_rel-1)*3+1:1+(n_rel-1)*3+3);
x_rwr = x_mot(:, 1+(n_rwr-1)*3+1:1+(n_rwr-1)*3+3);

x_rsh1 = x_mot(:, 1+(n_rsh1-1)*3+1:1+(n_rsh1-1)*3+3);
x_rsh2 = x_mot(:, 1+(n_rsh2-1)*3+1:1+(n_rsh2-1)*3+3);
x_rel1 = x_mot(:, 1+(n_rel1-1)*3+1:1+(n_rel1-1)*3+3);
x_rel2 = x_mot(:, 1+(n_rel2-1)*3+1:1+(n_rel2-1)*3+3);
x_rwr1 = x_mot(:, 1+(n_rwr1-1)*3+1:1+(n_rwr1-1)*3+3);
x_rwr2 = x_mot(:, 1+(n_rwr2-1)*3+1:1+(n_rwr2-1)*3+3);
x_rhd = x_mot(:, 1+(n_rhd-1)*3+1:+(n_rhd-1)*3+3);

nFrames = size(x_mot, 1);

Rut = zeros(nFrames, 3, 3);
Rua = zeros(nFrames, 3, 3);
Rfa = zeros(nFrames, 3, 3);
Rha = zeros(nFrames, 3, 3);
Rotm_gl = zeros(nFrames, 3, 3);
for iFrame = 1:nFrames
    Rotm_gl(iFrame, :, :) = [0, -1, 0; ...
    1, 0, 0; ...
    0, 0, 1];
end 

% rotation matrix of upper trunk 
x_ut = ((x_xip-x_xipb) + (x_clav-x_c7))/2;
z_ut = (x_clav+x_c7)/2 - (x_xip+x_xipb)/2; z_ut = z_ut./vecnorm(z_ut, 2, 2);
y_ut = cross(z_ut, x_ut); y_ut = y_ut./vecnorm(y_ut, 2, 2);
x_ut = cross(y_ut, z_ut);
for iFrame = 1:nFrames
    Rut(iFrame,:,:) = [x_ut(iFrame,:)', y_ut(iFrame,:)', z_ut(iFrame,:)'];
end

% rotation matrix of upper arm
z_ua = x_rsh - x_rel; z_ua = z_ua./vecnorm(z_ua, 2, 2);
y_ua = x_rel2 - x_rel1;
x_ua = cross(y_ua, z_ua); x_ua = x_ua./vecnorm(x_ua, 2, 2);
y_ua = cross(z_ua, x_ua);
for iFrame = 1:nFrames
    Rua(iFrame,:,:) = [x_ua(iFrame,:)', y_ua(iFrame,:)', z_ua(iFrame,:)'];
end

% rotation matrix of forearm
z_fa = x_rel - x_rwr; z_fa = z_fa./vecnorm(z_fa, 2, 2);
y_fa = x_rwr2 - x_rwr1;
x_fa = cross(y_fa, z_fa); x_fa = x_fa./vecnorm(x_fa, 2, 2);
y_fa = cross(z_fa, x_fa);
for iFrame = 1:nFrames
    Rfa(iFrame,:,:) = [x_fa(iFrame,:)', y_fa(iFrame,:)', z_fa(iFrame,:)'];
end

% rotation matrix of hand
z_ha = x_rwr - x_rhd; z_ha = z_ha./vecnorm(z_ha, 2, 2);
y_ha = x_rwr2 - x_rwr1;
x_ha = cross(y_ha, z_ha); x_ha = x_ha./vecnorm(x_ha, 2, 2);
y_ha = cross(z_ha, x_ha);
for iFrame = 1:nFrames
    Rha(iFrame,:,:) = [x_ha(iFrame,:)', y_ha(iFrame,:)', z_ha(iFrame,:)'];
end

[q_ut] = eua2(Rotm_gl, Rut);
[q_sh] = eua2(Rut, Rua, 'xyz');
[q_el] = eua2(Rua, Rfa);
[q_wr] = eua2(Rfa, Rha);


end

