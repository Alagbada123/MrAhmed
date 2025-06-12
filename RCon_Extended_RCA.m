%% RCon_Model_RCA.m
% This code is for modelling the homogenized stiffness and strength
% prediction of concrete with the composition of brick aggregates.
% VERSION: Incorporates Extended Powers' Model for both paste properties
% and final concrete volume fractions.

clearvars; close all; clc;
addpath('functions');
addpath('preCAL');
addpath('experimental_data');

%% 1). LOAD INPUT

% precision
steps_angleHYD=11;
steps_angleCON=31;
steps_rHZ=11;
choose_Eonly=0; % 1...only elasticity prediction, 0... also strength

% initialization of output
Estruc=struct('E',NaN,'nu',NaN);
scaleEstruct=struct('bcp',Estruc,'Icp1',Estruc,'Icp2',Estruc,'Icp3',Estruc,'conc',Estruc);
critstruct=struct('SIult',NaN,'I',NaN,'zeniHYD',NaN,'aziCON',NaN,'zeni_con',NaN,'r',NaN);
RCsol=struct('failure',critstruct,'siDP',NaN,'crit',NaN,'vol',NaN,'elasticity',scaleEstruct,'Rlist',NaN,'discrete',NaN);
critstruct=struct('maxSI',NaN,'zeniHYD',NaN,'aziCON',NaN,'zeniCON',NaN,'r',NaN);
RCsol.crit=repmat({critstruct}, 1,6);
siDPstruct=zeros(steps_angleHYD,steps_angleCON);
siDPstruct2=zeros(steps_angleHYD,steps_angleCON,steps_angleCON,steps_rHZ);
siDP_cell=repmat({siDPstruct}, 1,5);
siDP_cell{2}=siDPstruct2;siDP_cell{3}=siDPstruct2;siDP_cell{4}=siDPstruct2;siDP_cell{5}=siDPstruct2;
RCsol.siDP=siDP_cell;
RCsol.discrete=struct('angleHYD',NaN,'angleCON',NaN,'r',NaN);


% --- MODIFIED: required parameters to extract needed properties ---
req_wc = 0.55;
req_xi = 0.05:0.05:1;
req_Fpor = [1.0, 1.4, 1.3, 1.2]; % Fpor values for ITZ layers

%% --- NEW: Define which Extended Model case to run ---
% These values MUST exist in the pre-calculated lists
req_wa0a = 0.0;
req_alpha = 0.0;
fair_conc = 0.02; % Assumed volume fraction of entrapped air in concrete (e.g., 2%)
% To run the "Normal" Powers' model, set req_wa0a=0.0 and req_alpha=0.0
% ---

fileoutput='sensFpor_extended_v2'; % Give a new output name
RC_sens=cell(length(req_wc),length(req_xi),length(req_Fpor));
strength_values = NaN(1, length(req_xi));
compressive_strength = NaN(length(req_xi), 5);

% --- MODIFIED: load precalculated BCP properties (from "full" / no Fpor run) ---
disp('Loading BCP properties from pre-calculated extended model file...');
load('precalc_extended_model_full.mat');
output_cell_bcp = output_cell;
clear output_cell;
[~, ~, length_wa0a_bcp, length_alpha_bcp] = size(output_cell_bcp);
wc_calclist_bcp = cellfun(@(s) s.wc, output_cell_bcp(:,1,1,1));
xi_calclist_bcp = cellfun(@(s) s.xi, output_cell_bcp(1,:,1,1));
wa0a_calclist_bcp = cellfun(@(s) s.wa0_a, output_cell_bcp(1,1,:,1));
alpha_calclist_bcp = cellfun(@(s) s.alpha, output_cell_bcp(1,1,1,:));


% --- MODIFIED: load precalculated ITZ properties (from run with Fpor) ---
disp('Loading ITZ properties from pre-calculated extended model file...');
load('precalc_extended_model_with_Fpor_FULL.mat');
output_cell_itz = output_cell;
clear output_cell;
[~, ~, ~, length_wa0a_itz, length_alpha_itz] = size(output_cell_itz);
wc_calclist_itz = cellfun(@(s) s.wc, output_cell_itz(:,1,1,1,1));
xi_calclist_itz = cellfun(@(s) s.xi, output_cell_itz(1,:,1,1,1));
Fpor_calclist_itz = cellfun(@(s) s.Fpor, output_cell_itz(1,1,:,1,1));
wa0a_calclist_itz = cellfun(@(s) s.wa0_a, output_cell_itz(1,1,1,:,1));
alpha_calclist_itz = cellfun(@(s) s.alpha, output_cell_itz(1,1,1,1,:));


for wcit=1:length(req_wc)
    wc_ncp=req_wc(wcit);
    for xiit=1:length(req_xi)
        xi_ncp_max=min(1,wc_ncp/0.42);
        xi_ncp=min(req_xi(xiit),xi_ncp_max);
        
        disp(['STARTING WITH TEST wc=',num2str(wc_ncp),'  xi=',num2str(xi_ncp), ...
            '  wa0/a=', num2str(req_wa0a), '  alpha=', num2str(req_alpha)]);
        
        %% 2) HOMOGENIZATION OF CEMENT PASTE FROM PRECALCULATION
        
        % --- BCP properties ---
        disp('... extracting BCP properties');
        find_tol = 1e-8;
        indwc_bcp = find(abs(wc_calclist_bcp-wc_ncp) < find_tol, 1);
        indxi_bcp = find(abs(xi_calclist_bcp-xi_ncp) < find_tol, 1);
        indwa0a_bcp = find(abs(wa0a_calclist_bcp-req_wa0a) < find_tol, 1);
        indalpha_bcp = find(abs(alpha_calclist_bcp-req_alpha) < find_tol, 1);
        if isempty(indwc_bcp)||isempty(indxi_bcp)||isempty(indwa0a_bcp)||isempty(indalpha_bcp), error('Could not find BCP parameters.'); end
        bcp_data = output_cell_bcp{indwc_bcp, indxi_bcp, indwa0a_bcp, indalpha_bcp};
        Ebcp=bcp_data.calc_cp.E; nubcp=bcp_data.calc_cp.nu;
        diffQ_vol_bcp=bcp_data.calc_cp.diffQvol; diffQ_dev_bcp=bcp_data.calc_cp.diffQdev;
        fhf_bcp=bcp_data.vol_cp.hf; fhyddiff_bcp=bcp_data.numerics(2);
        [kbcp,mubcp]=fun_kmu_from_Enu(Ebcp,nubcp); Cbcp=fun_CfromEnu(Ebcp,nubcp);
        RCsol.elasticity.bcp.E=Ebcp; RCsol.elasticity.bcp.nu=nubcp;
        
        % --- ITZ properties ---
        disp('... extracting ITZ properties');
        indwc_itz = find(abs(wc_calclist_itz-wc_ncp) < find_tol, 1);
        indxi_itz = find(abs(xi_calclist_itz-xi_ncp) < find_tol, 1);
        indwa0a_itz = find(abs(wa0a_calclist_itz-req_wa0a) < find_tol, 1);
        indalpha_itz = find(abs(alpha_calclist_itz-req_alpha) < find_tol, 1);
        if isempty(indwc_itz)||isempty(indxi_itz)||isempty(indwa0a_itz)||isempty(indalpha_itz), error('Could not find ITZ parameters.'); end
        
        indFpor=find(abs(Fpor_calclist_itz-req_Fpor(2))<find_tol,1); if isempty(indFpor), error('Fpor not found for ITZ1'); end
        itz1_data = output_cell_itz{indwc_itz, indxi_itz, indFpor, indwa0a_itz, indalpha_itz};
        EIcp1=itz1_data.calc_cp.E; nuIcp1=itz1_data.calc_cp.nu; diffQ_vol_Icp1=itz1_data.calc_cp.diffQvol; diffQ_dev_Icp1=itz1_data.calc_cp.diffQdev; fhf_Icp1=itz1_data.vol_cp.hf; fhyddiff_Icp1=itz1_data.numerics(2);
        [kIcp1,muIcp1]=fun_kmu_from_Enu(EIcp1,nuIcp1); CIcp1=fun_CfromEnu(EIcp1,nuIcp1); RCsol.elasticity.Icp1.E=EIcp1; RCsol.elasticity.Icp1.nu=nuIcp1;

        indFpor=find(abs(Fpor_calclist_itz-req_Fpor(3))<find_tol,1); if isempty(indFpor), error('Fpor not found for ITZ2'); end
        itz2_data = output_cell_itz{indwc_itz, indxi_itz, indFpor, indwa0a_itz, indalpha_itz};
        EIcp2=itz2_data.calc_cp.E; nuIcp2=itz2_data.calc_cp.nu; diffQ_vol_Icp2=itz2_data.calc_cp.diffQvol; diffQ_dev_Icp2=itz2_data.calc_cp.diffQdev; fhf_Icp2=itz2_data.vol_cp.hf; fhyddiff_Icp2=itz2_data.numerics(2);
        [kIcp2,muIcp2]=fun_kmu_from_Enu(EIcp2,nuIcp2); CIcp2=fun_CfromEnu(EIcp2,nuIcp2); RCsol.elasticity.Icp2.E=EIcp2; RCsol.elasticity.Icp2.nu=nuIcp2;
 
        indFpor=find(abs(Fpor_calclist_itz-req_Fpor(4))<find_tol,1); if isempty(indFpor), error('Fpor not found for ITZ3'); end
        itz3_data = output_cell_itz{indwc_itz, indxi_itz, indFpor, indwa0a_itz, indalpha_itz};
        EIcp3=itz3_data.calc_cp.E; nuIcp3=itz3_data.calc_cp.nu; diffQ_vol_Icp3=itz3_data.calc_cp.diffQvol; diffQ_dev_Icp3=itz3_data.calc_cp.diffQdev; fhf_Icp3=itz3_data.vol_cp.hf; fhyddiff_Icp3=itz3_data.numerics(2);
        [kIcp3,muIcp3]=fun_kmu_from_Enu(EIcp3,nuIcp3); CIcp3=fun_CfromEnu(EIcp3,nuIcp3); RCsol.elasticity.Icp3.E=EIcp3; RCsol.elasticity.Icp3.nu=nuIcp3;

        %% --- MODIFIED: Section 3 using Extended Model (Eqs. 11-13) ---
        disp('... calculating concrete volume fractions using extended model');
        
        % DENSITIES AND MASS RATIOS
        rhoH2O=1000; rhocem=3150; rho_sand=2240; rho_ragg=2242;
        sc=3.522; % sand-to-cement mass ratio
        raggc=1.288; % recycled aggregate-to-cement mass ratio
        
        % 1. Calculate the initial effective w/c ratio for the entire concrete mix
        ac_total = sc + raggc; % Total aggregate to cement mass ratio
        wc0_eff_concrete = wc_ncp - req_wa0a * ac_total;
        if wc0_eff_concrete <= 0, error('Initial effective w/c ratio for the concrete mix is non-positive.'); end
        
        % 2. Calculate the volumes of each component (normalized by cement mass 'c')
        V_cp_norm = (1/rhocem) + (wc0_eff_concrete / rhoH2O);
        V_sand_norm = sc / rho_sand;
        V_ragg_norm = raggc / rho_ragg;
        
        % 3. Calculate the "no-air" volume fractions, as per Eq. 11
        V_total_no_air = V_cp_norm + V_sand_norm + V_ragg_norm;
        fcp_conc_no_air = V_cp_norm / V_total_no_air;
        fsand_conc_no_air = V_sand_norm / V_total_no_air;
        fragg_conc_no_air = V_ragg_norm / V_total_no_air;
        
        % 4. Adjust for entrapped air to get the final volume fractions, as per Eq. 13
        fcp_conc = fcp_conc_no_air * (1 - fair_conc);
        fsand_conc = fsand_conc_no_air * (1 - fair_conc);
        fragg_conc = fragg_conc_no_air * (1 - fair_conc);
        
        % 5. The rest of the logic for ITZ layers remains the same
        fragg_ITZ = calculate_ITZ_volume_fraction(fragg_conc, 4, 16, 0.03, 0);
        fITZ_layer = fragg_ITZ/3;
        fbcp_conc = fcp_conc - fragg_ITZ;
        
        %% 4) HOMOGENIZATION OF RECYCLED CONCRETE
        Eclin=139.9; nuclin=0.3; [kclin,muclin]=fun_kmu_from_Enu(Eclin,nuclin); Cclin=fun_CfromEnu(Eclin,nuclin);
        Ehyd=29.15786664; nuhyd=0.24; [khyd,muhyd]=fun_kmu_from_Enu(Ehyd,nuhyd); Chyd=fun_CfromEnu(Ehyd,nuhyd);
        ENAC=47.66 ; nuNAC=0.159; [kNAC,muNAC]=fun_kmu_from_Enu(ENAC,nuNAC); CNAC=fun_Cfromkmu(kNAC,muNAC);
        kNA=35.35; muNA=29.91; CNA = fun_Cfromkmu(kNA,muNA);
        I = eye(6); J = (1/3) * [ones(3,3), zeros(3,3); zeros(3,6)]; K = I - J;
        R_ragg = (1/2)*6^(1/3)*fragg_conc^(1/3)/pi^(1/3);
        R_ITZ1 = (1/2)*6^(1/3)*(fragg_conc+fITZ_layer)^(1/3)/pi^(1/3);
        R_ITZ2 = (1/2)*6^(1/3)*(fragg_conc+fITZ_layer*2)^(1/3)/pi^(1/3);
        R_ITZ3 = (1/2)*6^(1/3)*(fragg_conc+fITZ_layer*3)^(1/3)/pi^(1/3);
        R_bcp = (1/2)*6^(1/3)*(fragg_conc+fragg_ITZ+fbcp_conc)^(1/3)/pi^(1/3);
        RCsol.Rlist=[R_ragg,R_ITZ1,R_ITZ2,R_ITZ3,R_bcp];
        in_m=[kNAC,muNAC,R_ragg; kIcp1,muIcp1,R_ITZ1; kIcp2,muIcp2,R_ITZ2; kIcp3,muIcp3,R_ITZ3; kbcp,mubcp,R_bcp];
        sol_int=fun_HZ_int(in_m);   
        Ainf_ragg=sol_int(1,5)*J+sol_int(1,7)*K; Ainf_ITZ1=sol_int(2,5)*J+sol_int(2,7)*K;
        Ainf_ITZ2=sol_int(3,5)*J+sol_int(3,7)*K; Ainf_ITZ3=sol_int(4,5)*J+sol_int(4,7)*K; Ainf_bcp=I;
        P_sph=fun_P_sphere_iso(Cbcp); Ainf_NA=inv(I+P_sph*(CNA-Cbcp));
        EEinfty_con=inv(fsand_conc*Ainf_NA + fragg_conc*Ainf_ragg + fITZ_layer*Ainf_ITZ1 + fITZ_layer * Ainf_ITZ2 + fITZ_layer*Ainf_ITZ3 + fbcp_conc*Ainf_bcp);
        A_NA=Ainf_NA*EEinfty_con; A_ragg=Ainf_ragg*EEinfty_con; A_ITZ1=Ainf_ITZ1*EEinfty_con;
        A_ITZ2=Ainf_ITZ2*EEinfty_con; A_ITZ3=Ainf_ITZ3*EEinfty_con; A_bcp=Ainf_bcp*EEinfty_con;
        Ccon=fsand_conc*CNA*A_NA + fragg_conc*CNAC*A_ragg + fITZ_layer*CIcp1*A_ITZ1 + fITZ_layer*CIcp2*A_ITZ2 + fITZ_layer*CIcp3*A_ITZ3 + fbcp_conc*Cbcp*A_bcp;
        kcon = Ccon(1,1)-4/3*(Ccon(6,6)/2); mucon = (Ccon(6,6)/2);
        [Econ,nucon]=fun_Enu_from_kmu(kcon,mucon);
        RCsol.elasticity.conc.E=Econ; RCsol.elasticity.conc.nu=nucon;
        RC_sens{wcit, xiit, 1} = RCsol.elasticity.conc.E;

        %% 5) STRESS CONCENTRATION AND FAILURE ANALYSIS
        phi_hyd_degree=12; c_hyd=0.050; phi_hyd=phi_hyd_degree*pi/180;
        fc_MC=2*c_hyd*cos(phi_hyd)/(1-sin(phi_hyd));
        alpha_DP=sqrt(3)*fc_MC*tan(phi_hyd)/(3*c_hyd+fc_MC*tan(phi_hyd));
        k_DP=c_hyd*alpha_DP/tan(phi_hyd); fc_DP=3*k_DP/(sqrt(3)-alpha_DP);

        if choose_Eonly==0
            SI_macro=[0,0,-1,0,0,0]';
            angle_HYD=0:pi/(2*(steps_angleHYD-1)):pi/2;
            angle_CON=0:pi/(2*(steps_angleCON-1)):pi/2;
            RCsol.discrete.angleHYD=angle_HYD;
            RCsol.discrete.angleCON=angle_CON;
            B_sand=CNA*A_NA*inv(Ccon);
            BITZsand_sph=fun_BITZagg(kNA,kbcp,muNA,mubcp);
    
            for zeniHYDit=1:steps_angleHYD
                zeniL=angle_HYD(zeniHYDit); aziL=0; Q4=fun_Q4_bp(aziL,zeniL); Q4t=transpose(Q4);
                SI_macro_rot=Q4t*SI_macro; E_macro_rot=inv(Ccon)*SI_macro_rot;
                
                % 5a) failure around non-covered aggregates
                for zeniCONit=1:length(angle_CON)
                    psiL=angle_CON(zeniCONit); Q4=fun_Q4_bp(0,psiL); Q4t=transpose(Q4);
                    BITZsand=Q4t*BITZsand_sph*Q4;
                    eps_ITZ=inv(Cbcp)*BITZsand*B_sand*SI_macro_rot;
                    si_hyd_vol2=khyd*sqrt(3/(fhyddiff_bcp*fhf_bcp)*eps_ITZ'*(diffQ_vol_bcp*eps_ITZ));
                    si_hyd_dev2=muhyd*sqrt(2/(fhyddiff_bcp*fhf_bcp)*eps_ITZ'*(diffQ_dev_bcp*eps_ITZ));
                    siDP=si_hyd_dev2/sqrt(2)-alpha_DP*si_hyd_vol2/sqrt(3);
                    RCsol.siDP{1}(zeniHYDit,zeniCONit)=siDP;
                end
        
                % 5b) failure around covered aggregates
                E_HZ_rot=EEinfty_con*E_macro_rot;
                [eps_cart_cell,~,r_HZ,angle_HZ] = fun_HZ_field(in_m,E_HZ_rot,sol_int,steps_rHZ,steps_angleCON);
                RCsol.discrete.r = [r_HZ(1:steps_rHZ); r_HZ((steps_rHZ+1):(2*steps_rHZ)); r_HZ((2*steps_rHZ+1):(3*steps_rHZ)); r_HZ((3*steps_rHZ+1):(4*steps_rHZ)); r_HZ((4*steps_rHZ+1):(5*steps_rHZ))];
                
                for aziCONit=1:steps_angleCON
                    for zenCONit=1:steps_angleCON
                        for rit_con=1:steps_rHZ
                            eps_ITZ=eps_cart_cell{aziCONit,zenCONit}(:,steps_rHZ+rit_con);
                            si_hyd_vol2=khyd*sqrt(3/(fhyddiff_Icp1*fhf_Icp1)*eps_ITZ'*(diffQ_vol_Icp1*eps_ITZ)); si_hyd_dev2=muhyd*sqrt(2/(fhyddiff_Icp1*fhf_Icp1)*eps_ITZ'*(diffQ_dev_Icp1*eps_ITZ)); siDP=si_hyd_dev2/sqrt(2)-alpha_DP*si_hyd_vol2/sqrt(3); RCsol.siDP{2}(zeniHYDit,aziCONit,zenCONit,rit_con)=siDP;
                            
                            eps_ITZ=eps_cart_cell{aziCONit,zenCONit}(:,2*steps_rHZ+rit_con);
                            si_hyd_vol2=khyd*sqrt(3/(fhyddiff_Icp2*fhf_Icp2)*eps_ITZ'*(diffQ_vol_Icp2*eps_ITZ)); si_hyd_dev2=muhyd*sqrt(2/(fhyddiff_Icp2*fhf_Icp2)*eps_ITZ'*(diffQ_dev_Icp2*eps_ITZ)); siDP=si_hyd_dev2/sqrt(2)-alpha_DP*si_hyd_vol2/sqrt(3); RCsol.siDP{3}(zeniHYDit,aziCONit,zenCONit,rit_con)=siDP;
        
                            eps_ITZ=eps_cart_cell{aziCONit,zenCONit}(:,3*steps_rHZ+rit_con);
                            si_hyd_vol2=khyd*sqrt(3/(fhyddiff_Icp3*fhf_Icp3)*eps_ITZ'*(diffQ_vol_Icp3*eps_ITZ)); si_hyd_dev2=muhyd*sqrt(2/(fhyddiff_Icp3*fhf_Icp3)*eps_ITZ'*(diffQ_dev_Icp3*eps_ITZ)); siDP=si_hyd_dev2/sqrt(2)-alpha_DP*si_hyd_vol2/sqrt(3); RCsol.siDP{4}(zeniHYDit,aziCONit,zenCONit,rit_con)=siDP;
        
                            eps_ITZ=eps_cart_cell{aziCONit,zenCONit}(:,4*steps_rHZ+rit_con);
                            si_hyd_vol2=khyd*sqrt(3/(fhyddiff_bcp*fhf_bcp)*eps_ITZ'*(diffQ_vol_bcp*eps_ITZ)); si_hyd_dev2=muhyd*sqrt(2/(fhyddiff_bcp*fhf_bcp)*eps_ITZ'*(diffQ_dev_bcp*eps_ITZ)); siDP=si_hyd_dev2/sqrt(2)-alpha_DP*si_hyd_vol2/sqrt(3); RCsol.siDP{5}(zeniHYDit,aziCONit,zenCONit,rit_con)=siDP;
                        end
                    end
                end
            end
            
            % 6) get maxima of stresses and corresponding strength
            for i=1:5
                if i<=1, tmp=RCsol.siDP{i}(:,:); [RCsol.crit{i}.maxSI,k]=max(tmp(:)); [zeniHYDi,zeniCONi]=ind2sub(size(tmp),k); RCsol.crit{i}.zeniHYD=angle_HYD(zeniHYDi)*180/pi; RCsol.crit{i}.zeniCON=angle_CON(zeniCONi)*180/pi;
                else, tmp=RCsol.siDP{i}(:,:,:,:); [RCsol.crit{i}.maxSI,k]=max(tmp(:)); [zeniHYDi,aziCONi,zeniCONi,ri]=ind2sub(size(tmp),k); RCsol.crit{i}.zeniHYD=angle_HYD(zeniHYDi)*180/pi; RCsol.crit{i}.aziCON=angle_CON(aziCONi)*180/pi; RCsol.crit{i}.zeniCON=angle_CON(zeniCONi)*180/pi; RCsol.crit{i}.r=r_HZ(ri); end
            end
            
            tmp=[RCsol.crit{1}.maxSI,RCsol.crit{2}.maxSI,RCsol.crit{3}.maxSI,RCsol.crit{4}.maxSI,RCsol.crit{5}.maxSI];
            compressive_strength(xiit, :) = (k_DP ./ tmp) * 1000;
            % Corrected line
            [maxSI,k]=max(tmp); 
            RCsol.failure.SIult=k_DP/maxSI*1000;
            strength_values(xiit) = RCsol.failure.SIult;
            RCsol.failure.I=k; RCsol.failure.zeniHYD=RCsol.crit{k}.zeniHYD;
            RCsol.failure.aziCON=RCsol.crit{k}.aziCON; RCsol.failure.zeniCON=RCsol.crit{k}.zeniCON; RCsol.failure.r=RCsol.crit{k}.r;
        end
        
    end
end

%% 7) OUTPUT
if choose_Eonly==1
    filename = ['RC_Extended_Essd141312_',fileoutput,'.mat'];
    save(filename,'RC_sens');
else
    filename = ['RC_Extended_strengthssd141312_',fileoutput,'.mat'];
    save(filename,'strength_values','compressive_strength');
end

disp(['results for all compositions written to ',filename])

%% DISPLAY
RCsol.failure
RCsol.elasticity.conc