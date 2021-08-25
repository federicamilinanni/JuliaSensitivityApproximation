#############
# PKA model #
#############

function f_ODE(vf_, x_, p_, t)

    RiiP       = x_[1];
    RiiP_cAMP  = x_[2];
    RiiP_C     = x_[3];
    RiiP_C_cAMP = x_[4];
    C          = x_[5];
    Rii_cAMP   = x_[6];
    Rii_C_cAMP = x_[7];
    RiiP_CaN   = x_[8];
    RiiP_cAMP_CaN = x_[9];
    AKAR4_C    = x_[10];
    AKAR4p     = x_[11];
    kf_Rii_C__RiiP_C = p_[1];
    kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[2];
    kb_RiiP_CxcAMP__RiiP_C_cAMP = p_[3];
    kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[4];
    kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[5];
    kb_RiiPXcAMP__RiiP_cAMP = p_[6];
    kf_RiiPXcAMP__RiiP_cAMP = p_[7];
    kf_RiiPxC__RiiP_C = p_[8];
    kb_RiiPxC__RiiP_C = p_[9];
    kf_cAMPxRii__Rii_cAMP = p_[10];
    kb_cAMPxRii__Rii_cAMP = p_[11];
    kf_Rii_CxcAMP__Rii_C_cAMP = p_[12];
    kb_Rii_CxcAMP__Rii_C_cAMP = p_[13];
    kf_RiixC__Rii_C = p_[14];
    kf_Rii_cAMPxC__Rii_C_cAMP = p_[15];
    kb_Rii_cAMPxC__Rii_C_cAMP = p_[16];
    kf_Rii_C_cAMP__RiiP_C_cAMP = p_[17];
    kb_RiixC__Rii_C = p_[18];
    kf_C_AKAR4 = p_[19];
    kb_C_AKAR4 = p_[20];
    kcat_AKARp = p_[21];
    k18off     = p_[22];
    k20off     = p_[23];
    k21off     = p_[24];
    k23off     = p_[25];
    k18on      = p_[26];
    k20on      = p_[27];
    k21on      = p_[28];
    k23on      = p_[29];
    AKAPon     = p_[30];
    AKAR4_ConservedConst = p_[31];
    CaN_ConservedConst = p_[32];
    Rii_C_ConservedConst = p_[33];
    cAMP_ConservedConst = p_[34];
    Rii_ConservedConst = p_[35];
    AKAR4 = -AKAR4_C-AKAR4p+AKAR4_ConservedConst;
    CaN = -RiiP_cAMP_CaN+CaN_ConservedConst-RiiP_CaN;
    Rii_C = -Rii_C_cAMP-AKAR4_C+Rii_C_ConservedConst-RiiP_C-C-RiiP_C_cAMP;
    cAMP = -Rii_C_cAMP-Rii_cAMP-RiiP_cAMP_CaN-RiiP_cAMP+cAMP_ConservedConst-RiiP_C_cAMP;
    Rii = AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP;
    kf_RiiP_cAMP_CaN__CaNXRii_cAMP = AKAPon*k18on-k18off*(-1+AKAPon);
    kb_RiiPxCaN__RiiP_CaN = k20on*AKAPon-k20off*(-1+AKAPon);
    kf_RiiP_CaN__RiixCaN = AKAPon*k21on-k21off*(-1+AKAPon);
    kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = -k23off*(-1+AKAPon)+AKAPon*k23on;
    Km1 = 1;
    Km2 = 100;
    kf_RiiPxCaN__RiiP_CaN = (kb_RiiPxCaN__RiiP_CaN+kf_RiiP_cAMP_CaN__CaNXRii_cAMP)*Km1^(-1);
    kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = Km2^(-1)*(kf_RiiP_CaN__RiixCaN+kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN);
    reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
    reaction_14 = -RiiP_C*kb_RiiPxC__RiiP_C+kf_RiiPxC__RiiP_C*C*RiiP;
    reaction_12 = -kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP+RiiP_C*cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP-RiiP_cAMP*kb_RiiPXcAMP__RiiP_cAMP;
    reaction_23 = -kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP+kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C;
    reaction_78 = Rii*kf_cAMPxRii__Rii_cAMP*cAMP-kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
    reaction_56 = cAMP*Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP-Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP;
    reaction_76 = Rii_cAMP*C*kf_Rii_cAMPxC__Rii_C_cAMP-Rii_C_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP;
    reaction_62 = Rii_C_cAMP*kf_Rii_C_cAMP__RiiP_C_cAMP;
    reaction_58 = -Rii_C*kb_RiixC__Rii_C+kf_RiixC__Rii_C*Rii*C;
    reaction_44_ = CaN*RiiP*kf_RiiPxCaN__RiiP_CaN-kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
    reaction_33_ = CaN*RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN-RiiP_cAMP_CaN*kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN;
    reaction_4_8 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
    reaction_3_7 = RiiP_cAMP_CaN*kf_RiiP_cAMP_CaN__CaNXRii_cAMP;
    reaction_1 = -AKAR4_C*kb_C_AKAR4+AKAR4*kf_C_AKAR4*C;
    reaction_2 = AKAR4_C*kcat_AKARp;

    vf_[1] = -reaction_44_-reaction_43-reaction_14;
    vf_[2] = -reaction_23+reaction_43-reaction_33_;
    vf_[3] = reaction_51-reaction_12+reaction_14;
    vf_[4] = reaction_23+reaction_12+reaction_62;
    vf_[5] = -reaction_23-reaction_76-reaction_58+reaction_2-reaction_1-reaction_14;
    vf_[6] = -reaction_76+reaction_3_7+reaction_78;
    vf_[7] = reaction_76+reaction_56-reaction_62;
    vf_[8] = reaction_44_-reaction_4_8;
    vf_[9] = -reaction_3_7+reaction_33_;
    vf_[10] = -reaction_2+reaction_1;
    vf_[11] = reaction_2;


return vf_
end

##########################################

function f_ODE(x_, p_, t)
vf_ = similar(x_);
    RiiP       = x_[1];
    RiiP_cAMP  = x_[2];
    RiiP_C     = x_[3];
    RiiP_C_cAMP = x_[4];
    C          = x_[5];
    Rii_cAMP   = x_[6];
    Rii_C_cAMP = x_[7];
    RiiP_CaN   = x_[8];
    RiiP_cAMP_CaN = x_[9];
    AKAR4_C    = x_[10];
    AKAR4p     = x_[11];
    kf_Rii_C__RiiP_C = p_[1];
    kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[2];
    kb_RiiP_CxcAMP__RiiP_C_cAMP = p_[3];
    kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[4];
    kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[5];
    kb_RiiPXcAMP__RiiP_cAMP = p_[6];
    kf_RiiPXcAMP__RiiP_cAMP = p_[7];
    kf_RiiPxC__RiiP_C = p_[8];
    kb_RiiPxC__RiiP_C = p_[9];
    kf_cAMPxRii__Rii_cAMP = p_[10];
    kb_cAMPxRii__Rii_cAMP = p_[11];
    kf_Rii_CxcAMP__Rii_C_cAMP = p_[12];
    kb_Rii_CxcAMP__Rii_C_cAMP = p_[13];
    kf_RiixC__Rii_C = p_[14];
    kf_Rii_cAMPxC__Rii_C_cAMP = p_[15];
    kb_Rii_cAMPxC__Rii_C_cAMP = p_[16];
    kf_Rii_C_cAMP__RiiP_C_cAMP = p_[17];
    kb_RiixC__Rii_C = p_[18];
    kf_C_AKAR4 = p_[19];
    kb_C_AKAR4 = p_[20];
    kcat_AKARp = p_[21];
    k18off     = p_[22];
    k20off     = p_[23];
    k21off     = p_[24];
    k23off     = p_[25];
    k18on      = p_[26];
    k20on      = p_[27];
    k21on      = p_[28];
    k23on      = p_[29];
    AKAPon     = p_[30];
    AKAR4_ConservedConst = p_[31];
    CaN_ConservedConst = p_[32];
    Rii_C_ConservedConst = p_[33];
    cAMP_ConservedConst = p_[34];
    Rii_ConservedConst = p_[35];
    AKAR4 = -AKAR4_C-AKAR4p+AKAR4_ConservedConst;
    CaN = -RiiP_cAMP_CaN+CaN_ConservedConst-RiiP_CaN;
    Rii_C = -Rii_C_cAMP-AKAR4_C+Rii_C_ConservedConst-RiiP_C-C-RiiP_C_cAMP;
    cAMP = -Rii_C_cAMP-Rii_cAMP-RiiP_cAMP_CaN-RiiP_cAMP+cAMP_ConservedConst-RiiP_C_cAMP;
    Rii = AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP;
    kf_RiiP_cAMP_CaN__CaNXRii_cAMP = AKAPon*k18on-k18off*(-1+AKAPon);
    kb_RiiPxCaN__RiiP_CaN = k20on*AKAPon-k20off*(-1+AKAPon);
    kf_RiiP_CaN__RiixCaN = AKAPon*k21on-k21off*(-1+AKAPon);
    kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = -k23off*(-1+AKAPon)+AKAPon*k23on;
    Km1 = 1;
    Km2 = 100;
    kf_RiiPxCaN__RiiP_CaN = (kb_RiiPxCaN__RiiP_CaN+kf_RiiP_cAMP_CaN__CaNXRii_cAMP)*Km1^(-1);
    kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = Km2^(-1)*(kf_RiiP_CaN__RiixCaN+kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN);
    reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
    reaction_14 = -RiiP_C*kb_RiiPxC__RiiP_C+kf_RiiPxC__RiiP_C*C*RiiP;
    reaction_12 = -kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP+RiiP_C*cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP-RiiP_cAMP*kb_RiiPXcAMP__RiiP_cAMP;
    reaction_23 = -kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP+kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C;
    reaction_78 = Rii*kf_cAMPxRii__Rii_cAMP*cAMP-kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
    reaction_56 = cAMP*Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP-Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP;
    reaction_76 = Rii_cAMP*C*kf_Rii_cAMPxC__Rii_C_cAMP-Rii_C_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP;
    reaction_62 = Rii_C_cAMP*kf_Rii_C_cAMP__RiiP_C_cAMP;
    reaction_58 = -Rii_C*kb_RiixC__Rii_C+kf_RiixC__Rii_C*Rii*C;
    reaction_44_ = CaN*RiiP*kf_RiiPxCaN__RiiP_CaN-kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
    reaction_33_ = CaN*RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN-RiiP_cAMP_CaN*kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN;
    reaction_4_8 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
    reaction_3_7 = RiiP_cAMP_CaN*kf_RiiP_cAMP_CaN__CaNXRii_cAMP;
    reaction_1 = -AKAR4_C*kb_C_AKAR4+AKAR4*kf_C_AKAR4*C;
    reaction_2 = AKAR4_C*kcat_AKARp;
    vf_ = zeros(11,1);
    vf_[1] = -reaction_44_-reaction_43-reaction_14;
    vf_[2] = -reaction_23+reaction_43-reaction_33_;
    vf_[3] = reaction_51-reaction_12+reaction_14;
    vf_[4] = reaction_23+reaction_12+reaction_62;
    vf_[5] = -reaction_23-reaction_76-reaction_58+reaction_2-reaction_1-reaction_14;
    vf_[6] = -reaction_76+reaction_3_7+reaction_78;
    vf_[7] = reaction_76+reaction_56-reaction_62;
    vf_[8] = reaction_44_-reaction_4_8;
    vf_[9] = -reaction_3_7+reaction_33_;
    vf_[10] = -reaction_2+reaction_1;
    vf_[11] = reaction_2;


return vf_
end

##########################################

function Jacobian_x(x_, p_, t)
    RiiP       = x_[1];
    RiiP_cAMP  = x_[2];
    RiiP_C     = x_[3];
    RiiP_C_cAMP = x_[4];
    C          = x_[5];
    Rii_cAMP   = x_[6];
    Rii_C_cAMP = x_[7];
    RiiP_CaN   = x_[8];
    RiiP_cAMP_CaN = x_[9];
    AKAR4_C    = x_[10];
    AKAR4p     = x_[11];
    kf_Rii_C__RiiP_C = p_[1];
    kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[2];
    kb_RiiP_CxcAMP__RiiP_C_cAMP = p_[3];
    kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[4];
    kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[5];
    kb_RiiPXcAMP__RiiP_cAMP = p_[6];
    kf_RiiPXcAMP__RiiP_cAMP = p_[7];
    kf_RiiPxC__RiiP_C = p_[8];
    kb_RiiPxC__RiiP_C = p_[9];
    kf_cAMPxRii__Rii_cAMP = p_[10];
    kb_cAMPxRii__Rii_cAMP = p_[11];
    kf_Rii_CxcAMP__Rii_C_cAMP = p_[12];
    kb_Rii_CxcAMP__Rii_C_cAMP = p_[13];
    kf_RiixC__Rii_C = p_[14];
    kf_Rii_cAMPxC__Rii_C_cAMP = p_[15];
    kb_Rii_cAMPxC__Rii_C_cAMP = p_[16];
    kf_Rii_C_cAMP__RiiP_C_cAMP = p_[17];
    kb_RiixC__Rii_C = p_[18];
    kf_C_AKAR4 = p_[19];
    kb_C_AKAR4 = p_[20];
    kcat_AKARp = p_[21];
    k18off     = p_[22];
    k20off     = p_[23];
    k21off     = p_[24];
    k23off     = p_[25];
    k18on      = p_[26];
    k20on      = p_[27];
    k21on      = p_[28];
    k23on      = p_[29];
    AKAPon     = p_[30];
    AKAR4_ConservedConst = p_[31];
    CaN_ConservedConst = p_[32];
    Rii_C_ConservedConst = p_[33];
    cAMP_ConservedConst = p_[34];
    Rii_ConservedConst = p_[35];
    jac_ = zeros(11,11);
    jac_[1,1] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP-kf_RiiPxC__RiiP_C*C+(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon));
    jac_[1,2] = kf_RiiPXcAMP__RiiP_cAMP*RiiP+kb_RiiPXcAMP__RiiP_cAMP;
    jac_[1,3] = kb_RiiPxC__RiiP_C;
    jac_[1,4] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[1,5] = -kf_RiiPxC__RiiP_C*RiiP;
    jac_[1,6] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[1,7] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[1,8] = k20on*AKAPon+(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP-k20off*(-1+AKAPon);
    jac_[1,9] = kf_RiiPXcAMP__RiiP_cAMP*RiiP+(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP;
    jac_[2,1] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP;
    jac_[2,2] = 1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon))-kf_RiiPXcAMP__RiiP_cAMP*RiiP-kb_RiiPXcAMP__RiiP_cAMP-kf_RiiP_cAMPxC__RiiP_C_cAMP*C;
    jac_[2,4] = kb_RiiP_cAMPxC__RiiP_C_cAMP-kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[2,5] = -kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP;
    jac_[2,6] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[2,7] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jac_[2,8] = 1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jac_[2,9] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP-k23off*(-1+AKAPon)+AKAPon*k23on+1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jac_[3,1] = kf_RiiPxC__RiiP_C*C;
    jac_[3,2] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[3,3] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP-kf_Rii_C__RiiP_C-kb_RiiPxC__RiiP_C;
    jac_[3,4] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-kf_Rii_C__RiiP_C+kb_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[3,5] = kf_RiiPxC__RiiP_C*RiiP-kf_Rii_C__RiiP_C;
    jac_[3,6] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[3,7] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-kf_Rii_C__RiiP_C;
    jac_[3,9] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[3,10] = -kf_Rii_C__RiiP_C;
    jac_[4,2] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP+kf_RiiP_cAMPxC__RiiP_C_cAMP*C;
    jac_[4,3] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[4,4] = -kb_RiiP_cAMPxC__RiiP_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-kb_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[4,5] = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP;
    jac_[4,6] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[4,7] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP+kf_Rii_C_cAMP__RiiP_C_cAMP;
    jac_[4,9] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jac_[5,1] = kf_RiixC__Rii_C*C-kf_RiiPxC__RiiP_C*C;
    jac_[5,2] = kf_RiixC__Rii_C*C-kf_RiiP_cAMPxC__RiiP_C_cAMP*C;
    jac_[5,3] = kb_RiiPxC__RiiP_C-kb_RiixC__Rii_C;
    jac_[5,4] = kb_RiiP_cAMPxC__RiiP_C_cAMP-kb_RiixC__Rii_C;
    jac_[5,5] = -kf_RiiPxC__RiiP_C*RiiP-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-kf_RiixC__Rii_C*C-kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP-(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_RiixC__Rii_C+kf_C_AKAR4*(AKAR4_C+AKAR4p-AKAR4_ConservedConst)-kb_RiixC__Rii_C;
    jac_[5,6] = kf_RiixC__Rii_C*C-C*kf_Rii_cAMPxC__Rii_C_cAMP;
    jac_[5,7] = kb_Rii_cAMPxC__Rii_C_cAMP-kb_RiixC__Rii_C;
    jac_[5,8] = kf_RiixC__Rii_C*C;
    jac_[5,9] = kf_RiixC__Rii_C*C;
    jac_[5,10] = -kf_RiixC__Rii_C*C+kb_C_AKAR4+kf_C_AKAR4*C-kb_RiixC__Rii_C+kcat_AKARp;
    jac_[5,11] = kf_C_AKAR4*C;
    jac_[6,1] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP;
    jac_[6,2] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP-(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP;
    jac_[6,4] = -(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP;
    jac_[6,5] = -Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP;
    jac_[6,6] = -kb_cAMPxRii__Rii_cAMP+(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP-(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP-C*kf_Rii_cAMPxC__Rii_C_cAMP;
    jac_[6,7] = -(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP+kb_Rii_cAMPxC__Rii_C_cAMP;
    jac_[6,8] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP;
    jac_[6,9] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP+AKAPon*k18on-(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP-k18off*(-1+AKAPon);
    jac_[6,10] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP;
    jac_[7,2] = (Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[7,3] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[7,4] = (Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP+(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[7,5] = Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP+(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[7,6] = (Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP+C*kf_Rii_cAMPxC__Rii_C_cAMP;
    jac_[7,7] = -kf_Rii_C_cAMP__RiiP_C_cAMP+(Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP+(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP-kb_Rii_cAMPxC__Rii_C_cAMP;
    jac_[7,9] = (Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[7,10] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jac_[8,1] = -(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon));
    jac_[8,8] = -AKAPon*k21on-k20on*AKAPon-(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP+k21off*(-1+AKAPon)+k20off*(-1+AKAPon);
    jac_[8,9] = -(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP;
    jac_[9,2] = -1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jac_[9,8] = -1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jac_[9,9] = -AKAPon*k18on+k23off*(-1+AKAPon)-AKAPon*k23on-1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon))+k18off*(-1+AKAPon);
    jac_[10,5] = -kf_C_AKAR4*(AKAR4_C+AKAR4p-AKAR4_ConservedConst);
    jac_[10,10] = -kb_C_AKAR4-kf_C_AKAR4*C-kcat_AKARp;
    jac_[10,11] = -kf_C_AKAR4*C;
    jac_[11,10] = kcat_AKARp;


return jac_
end

##########################################

function Jacobian_p(x_, p_, t)
    RiiP       = x_[1];
    RiiP_cAMP  = x_[2];
    RiiP_C     = x_[3];
    RiiP_C_cAMP = x_[4];
    C          = x_[5];
    Rii_cAMP   = x_[6];
    Rii_C_cAMP = x_[7];
    RiiP_CaN   = x_[8];
    RiiP_cAMP_CaN = x_[9];
    AKAR4_C    = x_[10];
    AKAR4p     = x_[11];
    kf_Rii_C__RiiP_C = p_[1];
    kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[2];
    kb_RiiP_CxcAMP__RiiP_C_cAMP = p_[3];
    kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[4];
    kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[5];
    kb_RiiPXcAMP__RiiP_cAMP = p_[6];
    kf_RiiPXcAMP__RiiP_cAMP = p_[7];
    kf_RiiPxC__RiiP_C = p_[8];
    kb_RiiPxC__RiiP_C = p_[9];
    kf_cAMPxRii__Rii_cAMP = p_[10];
    kb_cAMPxRii__Rii_cAMP = p_[11];
    kf_Rii_CxcAMP__Rii_C_cAMP = p_[12];
    kb_Rii_CxcAMP__Rii_C_cAMP = p_[13];
    kf_RiixC__Rii_C = p_[14];
    kf_Rii_cAMPxC__Rii_C_cAMP = p_[15];
    kb_Rii_cAMPxC__Rii_C_cAMP = p_[16];
    kf_Rii_C_cAMP__RiiP_C_cAMP = p_[17];
    kb_RiixC__Rii_C = p_[18];
    kf_C_AKAR4 = p_[19];
    kb_C_AKAR4 = p_[20];
    kcat_AKARp = p_[21];
    k18off     = p_[22];
    k20off     = p_[23];
    k21off     = p_[24];
    k23off     = p_[25];
    k18on      = p_[26];
    k20on      = p_[27];
    k21on      = p_[28];
    k23on      = p_[29];
    AKAPon     = p_[30];
    AKAR4_ConservedConst = p_[31];
    CaN_ConservedConst = p_[32];
    Rii_C_ConservedConst = p_[33];
    cAMP_ConservedConst = p_[34];
    Rii_ConservedConst = p_[35];
    jacp_ = zeros(11,35);
    jacp_[1,6] = RiiP_cAMP;
    jacp_[1,7] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*RiiP;
    jacp_[1,8] = -C*RiiP;
    jacp_[1,9] = RiiP_C;
    jacp_[1,22] = -(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(-1+AKAPon)*RiiP;
    jacp_[1,23] = -RiiP_CaN*(-1+AKAPon)-(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(-1+AKAPon)*RiiP;
    jacp_[1,26] = AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP;
    jacp_[1,27] = AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP+AKAPon*RiiP_CaN;
    jacp_[1,30] = (k20on-k20off+k18on-k18off)*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP+(k20on-k20off)*RiiP_CaN;
    jacp_[1,32] = -(AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP;
    jacp_[1,34] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jacp_[2,4] = -RiiP_cAMP*C;
    jacp_[2,5] = RiiP_C_cAMP;
    jacp_[2,6] = -RiiP_cAMP;
    jacp_[2,7] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*RiiP;
    jacp_[2,24] = -1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP*(-1+AKAPon);
    jacp_[2,25] = -RiiP_cAMP_CaN*(-1+AKAPon)-1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP*(-1+AKAPon);
    jacp_[2,28] = 1/100*AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP;
    jacp_[2,29] = AKAPon*RiiP_cAMP_CaN+1/100*AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP;
    jacp_[2,30] = RiiP_cAMP_CaN*(k23on-k23off)-1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(k21off-k23on+k23off-k21on)*RiiP_cAMP;
    jacp_[2,32] = -1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jacp_[2,34] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
    jacp_[3,1] = -Rii_C_cAMP-AKAR4_C+Rii_C_ConservedConst-RiiP_C-C-RiiP_C_cAMP;
    jacp_[3,2] = RiiP_C*(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP);
    jacp_[3,3] = RiiP_C_cAMP;
    jacp_[3,8] = C*RiiP;
    jacp_[3,9] = -RiiP_C;
    jacp_[3,33] = kf_Rii_C__RiiP_C;
    jacp_[3,34] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jacp_[4,2] = -RiiP_C*(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP);
    jacp_[4,3] = -RiiP_C_cAMP;
    jacp_[4,4] = RiiP_cAMP*C;
    jacp_[4,5] = -RiiP_C_cAMP;
    jacp_[4,17] = Rii_C_cAMP;
    jacp_[4,34] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP;
    jacp_[5,4] = -RiiP_cAMP*C;
    jacp_[5,5] = RiiP_C_cAMP;
    jacp_[5,8] = -C*RiiP;
    jacp_[5,9] = RiiP_C;
    jacp_[5,14] = -(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*C;
    jacp_[5,15] = -Rii_cAMP*C;
    jacp_[5,16] = Rii_C_cAMP;
    jacp_[5,18] = -Rii_C_cAMP-AKAR4_C+Rii_C_ConservedConst-RiiP_C-C-RiiP_C_cAMP;
    jacp_[5,19] = C*(AKAR4_C+AKAR4p-AKAR4_ConservedConst);
    jacp_[5,20] = AKAR4_C;
    jacp_[5,21] = AKAR4_C;
    jacp_[5,31] = -kf_C_AKAR4*C;
    jacp_[5,33] = kb_RiixC__Rii_C;
    jacp_[5,35] = -kf_RiixC__Rii_C*C;
    jacp_[6,10] = -(AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP);
    jacp_[6,11] = -Rii_cAMP;
    jacp_[6,15] = -Rii_cAMP*C;
    jacp_[6,16] = Rii_C_cAMP;
    jacp_[6,22] = -RiiP_cAMP_CaN*(-1+AKAPon);
    jacp_[6,26] = AKAPon*RiiP_cAMP_CaN;
    jacp_[6,30] = RiiP_cAMP_CaN*(k18on-k18off);
    jacp_[6,34] = (AKAR4_C-Rii_cAMP-RiiP_cAMP_CaN+Rii_ConservedConst-RiiP_cAMP+C-RiiP_CaN-RiiP)*kf_cAMPxRii__Rii_cAMP;
    jacp_[6,35] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP;
    jacp_[7,12] = (Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*(Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP);
    jacp_[7,13] = -Rii_C_cAMP;
    jacp_[7,15] = Rii_cAMP*C;
    jacp_[7,16] = -Rii_C_cAMP;
    jacp_[7,17] = -Rii_C_cAMP;
    jacp_[7,33] = -(Rii_C_cAMP+Rii_cAMP+RiiP_cAMP_CaN+RiiP_cAMP-cAMP_ConservedConst+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jacp_[7,34] = -(Rii_C_cAMP+AKAR4_C-Rii_C_ConservedConst+RiiP_C+C+RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP;
    jacp_[8,22] = (RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(-1+AKAPon)*RiiP;
    jacp_[8,23] = RiiP_CaN*(-1+AKAPon)+(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(-1+AKAPon)*RiiP;
    jacp_[8,24] = RiiP_CaN*(-1+AKAPon);
    jacp_[8,26] = -AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP;
    jacp_[8,27] = -AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP-AKAPon*RiiP_CaN;
    jacp_[8,28] = -AKAPon*RiiP_CaN;
    jacp_[8,30] = -(k20on-k20off+k18on-k18off)*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP-(k20on-k20off)*RiiP_CaN+RiiP_CaN*(k21off-k21on);
    jacp_[8,32] = (AKAPon*k18on+k20on*AKAPon-k20off*(-1+AKAPon)-k18off*(-1+AKAPon))*RiiP;
    jacp_[9,22] = RiiP_cAMP_CaN*(-1+AKAPon);
    jacp_[9,24] = 1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP*(-1+AKAPon);
    jacp_[9,25] = RiiP_cAMP_CaN*(-1+AKAPon)+1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP*(-1+AKAPon);
    jacp_[9,26] = -AKAPon*RiiP_cAMP_CaN;
    jacp_[9,28] = -1/100*AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP;
    jacp_[9,29] = -AKAPon*RiiP_cAMP_CaN-1/100*AKAPon*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*RiiP_cAMP;
    jacp_[9,30] = -RiiP_cAMP_CaN*(k23on-k23off)-RiiP_cAMP_CaN*(k18on-k18off)+1/100*(RiiP_cAMP_CaN-CaN_ConservedConst+RiiP_CaN)*(k21off-k23on+k23off-k21on)*RiiP_cAMP;
    jacp_[9,32] = 1/100*RiiP_cAMP*(AKAPon*k21on-k23off*(-1+AKAPon)+AKAPon*k23on-k21off*(-1+AKAPon));
    jacp_[10,19] = -C*(AKAR4_C+AKAR4p-AKAR4_ConservedConst);
    jacp_[10,20] = -AKAR4_C;
    jacp_[10,21] = -AKAR4_C;
    jacp_[10,31] = kf_C_AKAR4*C;
    jacp_[11,21] = AKAR4_C;


return jacp_
end

p = [3.300E+01,
4.960E-01,
2.226E+00,
5.450E-02,
1.560E-01,
1.600E-03,
1.500E-02,
3.800E-02,
2.600E-03,
1.500E-01,
1.600E-03,
4.960E-01,
2.217E+00,
2.100E+00,
2.984E-01,
1.800E-02,
3.300E+01,
3.000E-04,
1.800E-02,
1.000E+00,
1.028E+01,
2.600E+00,
2.000E+01,
2.600E+00,
2.000E+01,
3.300E-01,
2.000E+00,
3.300E-01,
2.000E+00,
0,
2.0,
0,
0.6,
1.0,
6.3]

x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
