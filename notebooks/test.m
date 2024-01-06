%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) NICKEL MINING: SUPPLY CHAIN OF METALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except L Y Y_region A Ext d index_matrix_P Y_global year

% nickel mining
z = 0;
for n = 1:189;
for i = 27;
q = i + (163 * (n-1));
z= z+1;
index_nickel(z)= q;
end
end

d_BD = zeros(1,30807);  % d_BD = np.zeros(30807)
d_BD(1,index_nickel) = d(4,index_nickel);  % d_BD[index_nickel - 1] = d[3, index_nickel - 1]

% DERIVE INDEX FOR TARGET ETC
index_t_r = [1:189];

metals = [24:31,69:82]; 
index_t_s = [metals];

% derive the index of target-sector-regions (t) and non-target sector-regions (o)
index_all_s = [1:163];
index_all_r = [1:189];

index_o_s = setdiff(index_all_s, index_t_s);
index_o_r = setdiff(index_all_r, index_t_r);

z = 0;
for n = index_t_r;
for i = index_t_s;
q = i + (163 * (n-1));
z= z+1;
index_t(z)= q;
end
end
n_t = length(index_t);

index_all = [1:30807];
index_o = setdiff(index_all, index_t);
n_o = length(index_o);

n_t_s = length(index_t_s);
n_t_r = length(index_t_r);
n_o_s = length(index_o_s);

% Derive the Leontief inverse for the remaining economy 
I = eye(size(A));
L_oo_dash = zeros(30807,30807);
L_oo_dash(index_o,index_o) = inv(I(index_o,index_o)-A(index_o,index_o));


% Link Preg Treg
E_P_T = zeros(30807,30807);
E_P_T(:,index_t) = diag(d_BD) * L(:,index_t) * diag((Y_global(index_t) + A(index_t,index_o) * L_oo_dash(index_o,index_o) * Y_global(index_o)));

for r1 = 1:189
    for r2 = 1:189
        E_Preg_Treg(r1,r2) = sum(sum(E_P_T(index_matrix_P(r1,:),index_matrix_P(r2,:))));
    end
end

% Link Treg FDreg
E_T_FDreg = zeros(30807,189);
E_T_FDreg(index_t,:) = diag(d_BD * L(:,index_t)) * (Y_region(index_t,:) + A(index_t,index_o) * L_oo_dash(index_o,index_o) * Y_region(index_o,:));

for r1 = 1:189
        E_Treg_FDreg(r1,:) = sum(E_T_FDreg(index_matrix_P(r1,:),:));
end

% Link Preg FDreg 
E_P_FDreg = diag(d_BD) * L(:,index_t) * (Y_region(index_t,:) + A(index_t,index_o) * L_oo_dash(index_o,index_o) * Y_region(index_o,:));

for r1 = 1:189
        E_Preg_FDreg(r1,:) = sum(E_P_FDreg(index_matrix_P(r1,:),:));
end

% Link Tsec Psec 
for s1 = 1:163
    for s2 = 1:163
        E_Psec_Tsec(s1,s2) = sum(sum(E_P_T(index_matrix_P(:,s1),index_matrix_P(:,s2))));
    end
end

E_Tsec_Psec = E_Psec_Tsec';

% Link Psec Preg
E_P = sum(E_P_T,2);  % sum along the columns (i.e. 'axis=1' in numpy) --> column vector
E_Psec_Preg = vec2mat(E_P,163)'; % reshape to 189 (rows) x 163 (cols) and then transpose

% link FDreg FSsec 
E_FS_FDreg_t = zeros(30807,189);
E_FS_FDreg_o = zeros(30807,189);
E_FS_FDreg_t(index_t,:) = diag(d_BD * L(:,index_t)) * Y_region(index_t,:);
E_FS_FDreg_o(index_o,:) = diag(d_BD * L(:,index_t) * A(index_t,index_o) * L_oo_dash(index_o,index_o)) * Y_region(index_o,:);

E_FS_FDreg = E_FS_FDreg_t + E_FS_FDreg_o;

for s = 1:163
        E_FSsec_FDreg(s,:) = sum(E_FS_FDreg(index_matrix_P(:,s),:));  % index_matrix_P = sector_labels
end

E_FDreg_FSsec = E_FSsec_FDreg';

% Link Tsec Treg
E_T = sum(E_P_T)';  % sum along the rows (i.e. 'axis=0' in numpy) --> row vector which then transposed
E_Tsec_Treg = vec2mat(E_T,163)'; % reshape to 189 (rows) x 163 (cols) and then transpose

% store data 
datapath = ['BD_mining_study_results/Nickel_mining/'];
save([datapath 'E_Psec_Preg'],'E_Psec_Preg');
save([datapath 'E_Preg_Treg'],'E_Preg_Treg');
save([datapath 'E_Treg_FDreg'],'E_Treg_FDreg');
save([datapath 'E_Preg_FDreg'],'E_Preg_FDreg');
save([datapath 'E_Tsec_Psec'],'E_Tsec_Psec');
save([datapath 'E_FDreg_FSsec'],'E_FDreg_FSsec');
save([datapath 'E_Tsec_Treg'],'E_Tsec_Treg');

Compiled_reg(:,1) = sum(E_Preg_FDreg,2) ./ sum(sum(E_Preg_FDreg));
Compiled_reg(:,2) = sum(E_Treg_FDreg,2) ./ sum(sum(E_Treg_FDreg));
Compiled_reg(:,3) = sum(E_Treg_FDreg)' ./ sum(sum(E_Treg_FDreg));

New_Caledonia = 66;
Philippines = 70;
Cuba = 94;
Japan= 30;
Australia = 38;
China = [31,41];
USA = 29;
EU = 1:27;
Indonesia = 43;
Otherreg = setdiff(1:189,[New_Caledonia,Indonesia,Philippines,Cuba,Japan,Australia,China,USA,EU]);

% matrix of 189 rows x 10 columns, where columns is 1 for the selected region
Agg_Preg = zeros(189,10);
Agg_Preg(New_Caledonia,1) = 1;  % [1, 0, 0, 0 ... 0] (column vector)
Agg_Preg(Philippines,2) = 1;    % [0, 1, 0, 0 ... 0] (column vector)
Agg_Preg(Indonesia,3) = 1;      % [0, 0, 1, 0 ... 0] (column vector)
Agg_Preg(Cuba,4) = 1;
Agg_Preg(Japan,5) = 1;
Agg_Preg(Australia,6) = 1;
Agg_Preg(China,7) = 1;
Agg_Preg(USA,8) = 1;
Agg_Preg(EU,9) = 1;
Agg_Preg(Otherreg,10) = 1;


E_Preg_Treg_agg_shares = Agg_Preg' * E_Preg_Treg * Agg_Preg ./ sum(sum(E_Preg_Treg));
E_Treg_FDreg_agg_shares = Agg_Preg' * E_Treg_FDreg * Agg_Preg ./ sum(sum(E_Treg_FDreg));
E_Preg_FDreg_agg_shares = Agg_Preg' * E_Preg_FDreg * Agg_Preg ./ sum(sum(E_Preg_FDreg));

Compiled_sec(:,2) = sum(E_FDreg_FSsec)' ./ sum(sum(E_FDreg_FSsec));
Compiled_sec(:,1) = sum(E_Tsec_Treg,2) ./ sum(sum(E_Tsec_Treg));

% aggregate FSsec

construction = [112:113];
electronics = [84:87];
machinery = 83;
auto = [88:89];
service = [114:143];
food = [1:14,19,35:45];
manuf = [46:49,92];
otherFSsec = setdiff(1:163,[construction,electronics,machinery,auto,service,food,manuf]);

sum(sum(E_FDreg_FSsec(:,[service]))) ./ sum(sum(E_FDreg_FSsec))
sum(sum(E_FDreg_FSsec(:,[electronics]))) ./ sum(sum(E_FDreg_FSsec))
sum(sum(E_FDreg_FSsec(:,[machinery]))) ./ sum(sum(E_FDreg_FSsec))
sum(sum(E_FDreg_FSsec(:,[auto]))) ./ sum(sum(E_FDreg_FSsec))
sum(sum(E_FDreg_FSsec(:,[construction]))) ./ sum(sum(E_FDreg_FSsec))
sum(sum(E_FDreg_FSsec(:,[food,manuf,otherFSsec]))) ./ sum(sum(E_FDreg_FSsec))

Agg_FSsec = zeros(163,6);
Agg_FSsec([service],1) = 1;
Agg_FSsec(construction,2) = 1;
Agg_FSsec([electronics],3) = 1;
Agg_FSsec([machinery],5) = 1;
Agg_FSsec([auto],4) = 1;
Agg_FSsec([food,manuf,otherFSsec],6) = 1;

E_FDreg_FSsec_agg_shares = Agg_Preg' * E_FDreg_FSsec * Agg_FSsec ./ sum(sum(E_FDreg_FSsec));

% convert 2d array into a column vector
E_Preg_Treg_agg_shares_1D = E_Preg_Treg_agg_shares(:);
E_Treg_FDreg_agg_shares_1D = E_Treg_FDreg_agg_shares(:);
E_FDreg_FSsec_agg_shares_1D = E_FDreg_FSsec_agg_shares(:);
