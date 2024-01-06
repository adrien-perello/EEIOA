
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) SUPPLY CHAIN OF METALS: GLOBAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except L Y Y_region A Ext d index_matrix_P Y_global year

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

d_BD = d(4,:);  % ! DIFF COMPARED TO NICKEL  -->  d_BD = d[3, :] --> extra row

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


% link Preg vs FSsec
E_P_FS_t = zeros(30807,30807);   % ! DIFF COMPARED TO NICKEL zeros(30807,189)
E_P_FS_o = zeros(30807,30807);   % ! DIFF COMPARED TO NICKEL zeros(30807,189)
E_P_FS_t(:,index_t) = diag(d_BD) * L(:,index_t) * diag(Y_global(index_t));
E_P_FS_o(:,index_o) = diag(d_BD) * L(:,index_t) * A(index_t,index_o) * L_oo_dash(index_o,index_o) * diag(Y_global(index_o));

E_P_FS = E_P_FS_t + E_P_FS_o;

for s1 = 1:163
    for s2 = 1:163
        E_Psec_FSsec(s1,s2) = sum(sum(E_P_FS(index_matrix_P(:,s1),index_matrix_P(:,s2))));
    end
end

% Link Tsec Treg 
E_T = sum(E_P_T)';
E_Tsec_Treg = vec2mat(E_T,163)';

% store data 
datapath = ['BD_mining_study_results/Metals_production/'];
save([datapath 'E_Psec_Preg'],'E_Psec_Preg');
save([datapath 'E_Preg_Treg'],'E_Preg_Treg');
save([datapath 'E_Treg_FDreg'],'E_Treg_FDreg');
save([datapath 'E_Preg_FDreg'],'E_Preg_FDreg');
save([datapath 'E_Tsec_Psec'],'E_Tsec_Psec');
save([datapath 'E_FDreg_FSsec'],'E_FDreg_FSsec');
save([datapath 'E_Tsec_Treg'],'E_Tsec_Treg');
save([datapath 'E_Psec_FSsec'],'E_Psec_FSsec');  % ! save new section


%%
nickel = [27];
coal = [20];
alu = [28,74,75];
cu = [26,78,79];
iron = [25,70,71];
PM = [29,72,73];
metals = [24:31,69:82]; 
leadzinctin = [30,76,77];
othermet = setdiff(metals,[alu,cu,iron,nickel,PM,leadzinctin]);
index_o_s = setdiff(1:163, metals);
othermining = setdiff(index_o_s,[coal]);

Agg_Psec = zeros(163,7);
Agg_Psec(nickel,1) = 1;
Agg_Psec(PM,2) = 1;
Agg_Psec(alu,3) = 1;
Agg_Psec(iron,4) = 1;
Agg_Psec(cu,5) = 1;
Agg_Psec([othermet,leadzinctin],6) = 1;
Agg_Psec([coal,othermining],7) = 1;

Agg_Tsec = zeros(163,6);
Agg_Tsec(nickel,1) = 1;g
Agg_Tsec(PM,2) = 1;
Agg_Tsec(alu,3) = 1;
Agg_Tsec(iron,4) = 1;
Agg_Tsec(cu,5) = 1;
Agg_Tsec([othermet,leadzinctin],6) = 1;

%G20 = [1:38,40,41,43,44,82,186];
%nonG20 = setdiff(1:189,G20);
newcaledonia = 66;
australia = 38;
indonesia = 43;
brazil = 34;
philippines = 70;
peru = 107;
chile = 91;
china = [31,41];
japan = 30;
usa = 29;
eu = [1:27];

otherreg = setdiff(1:189,[newcaledonia,brazil,philippines,peru,chile,china,japan,indonesia,australia,usa,eu]);

Agg_Preg = zeros(189,7);
Agg_Preg(newcaledonia,1) = 1;
Agg_Preg(australia,2) = 1;
Agg_Preg(indonesia,3) = 1;
Agg_Preg(brazil,4) = 1;
%Agg_Preg(philippines,5) = 1;
%Agg_Preg(peru,6) = 1;
%Agg_Preg(chile,7) = 1;
Agg_Preg(china,5) = 1;
Agg_Preg(japan,6) = 1;
%Agg_Preg(usa,7) = 1;
%Agg_Preg(eu,8) = 1;
Agg_Preg([otherreg,philippines,peru,chile,eu,usa],7) = 1;


% aggregate FSsec
construction = [112:113];
electronics = [84:87];
machinery = 83;
auto = [88:89];
service = [114:143];
food = [1:14,19,35:45];
manuf = [46:49,93];
otherFSsec = setdiff(1:163,[construction,electronics,machinery,auto,service,food,manuf]);

Agg_FSsec = zeros(163,6);
Agg_FSsec(construction,1) = 1;
Agg_FSsec([service],2) = 1;
Agg_FSsec([machinery],5) = 1;
Agg_FSsec([electronics],3) = 1;
Agg_FSsec([auto],4) = 1;
Agg_FSsec([food,manuf,otherFSsec],6) = 1;

% Sankey no 1
TOT = sum(sum(E_Preg_Treg));

%%
%Sankey_Tsec_Psec = Agg_Psec' * E_Tsec_Psec' * Agg_Tsec;
%Sankey_Psec_Preg = Agg_Psec' * E_Psec_Preg * Agg_Preg;
Sankey_Preg_Treg = Agg_Preg' * E_Preg_Treg * Agg_Preg;
Sankey_Treg_FDreg = Agg_Preg' * E_Treg_FDreg * Agg_Preg;
Sankey_FDreg_FSsec = Agg_Preg' * E_FDreg_FSsec * Agg_FSsec;

%shares_Sankey_Tsec_Psec = Sankey_Tsec_Psec ./ TOT;
%shares_Sankey_Psec_Preg = Sankey_Psec_Preg ./ TOT;
shares_Sankey_Preg_Treg = Sankey_Preg_Treg ./ TOT;
shares_Sankey_Treg_FDreg = Sankey_Treg_FDreg ./ TOT;
shares_Sankey_FDreg_FSsec = Sankey_FDreg_FSsec ./ TOT;

Perspectives = [sum(shares_Sankey_Preg_Treg,2);sum(shares_Sankey_Treg_FDreg,2);sum(shares_Sankey_Treg_FDreg)';sum(shares_Sankey_FDreg_FSsec)'];
Linkages = [shares_Sankey_Preg_Treg(:);shares_Sankey_Treg_FDreg(:);shares_Sankey_FDreg_FSsec(:)];


%% Sankey no 2
Sankey_Preg_Psec = Agg_Preg' * E_Psec_Preg' * Agg_Psec;
Sankey_Psec_Tsec = Agg_Psec' * E_Tsec_Psec' * Agg_Tsec;
Sankey_Tsec_Treg = Agg_Tsec' * E_Tsec_Treg * Agg_Preg;

shares_Sankey_FSsec_Psec = Sankey_Preg_Psec ./ TOT;
shares_Sankey_Psec_Tsec = Sankey_Psec_Tsec ./ TOT;
shares_Sankey_Tsec_Treg = Sankey_Tsec_Treg ./ TOT;

Perspectives = [sum(shares_Sankey_FSsec_Psec,2);sum(shares_Sankey_Psec_Tsec,2);sum(shares_Sankey_Psec_Tsec)';sum(shares_Sankey_Tsec_Treg)'];
Linkages = [shares_Sankey_FSsec_Psec(:);shares_Sankey_Psec_Tsec(:);shares_Sankey_Tsec_Treg(:)];


%% Sankey no 3
Sankey_FSsec_Psec = Agg_FSsec' * E_Psec_FSsec' * Agg_Psec;
%Sankey_Preg_Psec = Agg_Preg' * E_Psec_Preg' * Agg_Psec;
Sankey_Psec_Tsec = Agg_Psec' * E_Tsec_Psec' * Agg_Tsec;
Sankey_Tsec_Treg = Agg_Tsec' * E_Tsec_Treg * Agg_Preg;

shares_Sankey_FSsec_Psec = Sankey_FSsec_Psec ./ TOT;
shares_Sankey_Psec_Tsec = Sankey_Psec_Tsec ./ TOT;
shares_Sankey_Tsec_Treg = Sankey_Tsec_Treg ./ TOT;

Perspectives = [sum(shares_Sankey_FSsec_Psec,2);sum(shares_Sankey_Psec_Tsec,2);sum(shares_Sankey_Psec_Tsec)';sum(shares_Sankey_Tsec_Treg)'];
Linkages = [shares_Sankey_FSsec_Psec(:);shares_Sankey_Psec_Tsec(:);shares_Sankey_Tsec_Treg(:)];

