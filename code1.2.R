results_df <- read.xlsx("最佳分布.xlsx")
#"exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1","FP2","RCS","RP-hazard","RP-odds","RP-normal","GAM","mix-cure"
os_cg <- c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1","RP-hazard","RP-odds","RP-normal")#
os_tg <- c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
pfs_cg <- c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1","RP-hazard","RP-odds","RP-normal")
pfs_tg <- c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma","FP1","FP2","RCS","RP-hazard","RP-odds","RP-normal","GAM","mix-cure")

surv_PFS_TG <- read.xlsx("surv_PFS_TG.xlsx")
surv_OS_TG <- read.xlsx("surv_OS_TG.xlsx")
surv_PFS_CG <- read.xlsx("surv_PFS_CG.xlsx")
surv_OS_CG <- read.xlsx("surv_OS_CG.xlsx")
# 初始化一个空列表来保存TG结果
TG_results <- list()

for(i in os_tg) {#
  for(j in pfs_tg) {#
    print(paste0(i, "_", j))
    Difference <- surv_OS_TG[, i] - surv_PFS_TG[, j]#
    if(any(Difference < 0)){
      print(paste0("PFS高于OS：",i, "_", j))
    }
    if(all(Difference >=0) == TRUE) {
      results_df[2, 2] <- i#4
      results_df[1, 2] <- j#3
      
      PFS <- surv_PFS_TG
      OS <- surv_OS_TG
      #head(PFS)
      #head(OS)
      #PFS$llogis <- pmin(PFS$llogis, OS$exp)
      #OS$exp <- pmax(OS$exp, PFS$llogis)
      
      ages <- 60  # 患者年龄
      md_pfs <- results_df$distribution[results_df$name == "dfGOF_PFS_TG.xlsx"]  # 选择模型分布名称
      md_os <- results_df$distribution[results_df$name == "dfGOF_OS_TG.xlsx"] 
      
      # 计算TG
      TG <- pXTX(path1, target_country, years, ages, t_cycle, cycle, PFS, OS)
      
      # 将TG结果保存到列表中
      TG_results[[paste0(i, "_", j)]] <- TG
    }
  }
}

CG_results <- list()
for(i in os_cg) {#
  for(j in pfs_cg) {#
    print(paste0(i, "_", j))
    Difference <- surv_OS_CG[, i] - surv_PFS_CG[, j]#
    if(any(Difference < 0)){
      print(paste0("PFS高于OS：",i, "_", j))
    }
    if(all(Difference >=0) == TRUE) {
      results_df[4, 2] <- i#4
      results_df[3, 2] <- j#3
      
      PFS <- surv_PFS_CG
      OS <- surv_OS_CG
      
      ages <- 60  # 患者年龄
      md_pfs <- results_df$distribution[results_df$name == "dfGOF_PFS_CG.xlsx"]  # 选择模型分布名称
      md_os <- results_df$distribution[results_df$name == "dfGOF_OS_CG.xlsx"] 
      
      # 计算TG
      CG <- pXTX(path1, target_country, years, ages, t_cycle, cycle, PFS, OS)
      
      # 将TG结果保存到列表中
      CG_results[[paste0(i, "_", j)]] <- CG
    }
  }
}

results1 <- list()
for(i in 1:34){
  for(j in 1:72){
    print(paste0("TG_",names(TG_results)[i], "_CG_", names(CG_results)[j]))
    TG <- TG_results[[i]]
    CG <- CG_results[[j]]
    
    TG$pPTP[1] <- 1
    CG$pPTP[1] <- 1
    
    mat_TG <- define_transition(
      state_names = c("pfs", "pd", "d"),
      TG[model_time,"pFTF"], TG[model_time,"pFTP"], TG[model_time,"pFTD"],
      0, TG[model_time,"pPTP"],TG[model_time,"pPTD"],
      0,0,1)
    mat_CG <- define_transition(
      state_names = c("pfs", "pd", "d"),
      CG[model_time,"pFTF"], CG[model_time,"pFTP"], CG[model_time,"pFTD"],
      0, CG[model_time,"pPTP"],CG[model_time,"pPTD"],
      0,0,1)
    
    state_pfs <- define_state(
      cost_treat = discount(
        dispatch_strategy(
          TG = ifelse(model_time <= 6, C_imu + C_oxa + C_cap + C_bev, C_cap + C_bev),
          CG = ifelse(model_time <= 6, C_oxa + C_cap + C_bev, C_cap + C_bev)
        ),
        r = dr / (365.24 / 21)
      ),
      # 计算治疗成本，并将其贴现到现值
      cost_other = discount(
        ifelse(
          model_time <= 18 * 2,
          ifelse((model_time - 1) %% 9 == 0, C_test, 0),
          ifelse((model_time - 1) %% 18 == 0, C_test, 0)
        ) +
          ifelse(
            model_time <= 18 * 2,
            ifelse((model_time - 1) %% 9 == 0, C_imag, 0),
            ifelse((model_time - 1) %% 18 == 0, C_imag, 0)
          ),
        r = dr / (365.24 / 21)
      ),
      cost_ae = discount(
        ifelse(
          model_time == 1,
          dispatch_strategy(
            TG = C_ane * R_ane2 + C_leuk * R_leuk2 + C_throm * R_throm2 + C_hfsyn * R_hfsyn2 + C_dia * R_dia2 + C_ormu * R_ormu2,
            CG = C_ane * R_ane1 + C_leuk * R_leuk1 + C_throm * R_throm1 + C_hfsyn * R_hfsyn1 + C_dia * R_dia1 + C_ormu * R_ormu1
          ),
          0
        ),
        r = dr / (365.24 / 21)
      ),
      cost_eol = 0,
      cost_total = cost_treat + cost_other + cost_ae + cost_eol,
      # 效用值计算
      utility = discount(
        ifelse(
          model_time == 1,
          dispatch_strategy(
            TG = (U_pfs2 - (U_ane * R_ane2 + U_leuk * R_leuk2 + U_throm * R_throm2 + U_hfsyn * R_hfsyn2 + U_diar * R_dia2 + U_ormu * R_ormu2)) / (365.24 / 21),
            CG = (U_pfs1 - (U_ane * R_ane1 + U_leuk * R_leuk1 + U_throm * R_throm1 + U_hfsyn * R_hfsyn1 + U_diar * R_dia1 + U_ormu * R_ormu1)) / (365.24 / 21)
          ),
          ifelse(
            TRUE,  # 只要 model_time != 1，直接返回以下两个值
            U_pfs2 / (365.24 / 21),
            U_pfs1 / (365.24 / 21)
          )
        ),
        r = dr / (365.24 / 21)
      )
    )
    # 定义疾病进展（PD）状态
    state_pd <- define_state(
      cost_treat = discount(
        dispatch_strategy(
          TG = C_bst * R_bst2 +
            ifelse((model_time - 1) %% 2 == 0, C_foriB * R_foriB2 + C_foriC * R_foriC2 + C_IriB * R_IriCaB2 + C_XeB * R_XeB2, 0) +
            ifelse((model_time) %% 4 == 0, 0, C_Fru * R_Fru2) +
            C_fori * R_fori2 + C_IriCaB * R_IriCaB2 + C_CaB * R_CaB2,
          CG = C_bst * R_bst1 +
            ifelse((model_time - 1) %% 2 == 0, C_foriB * R_foriB1 + C_foriC * R_foriC1 + C_IriB * R_IriCaB1 + C_XeB * R_XeB1, 0) +
            ifelse((model_time) %% 4 == 0, 0, C_Fru * R_Fru1) +
            C_fori * R_fori1 + C_IriCaB * R_IriCaB1 + C_CaB * R_CaB1
        ),
        r = dr / (365.24 / 21)
      ),
      cost_other = discount(
        ifelse(
          model_time <= 18 * 2,
          ifelse((model_time - 1) %% 9 == 0, C_test, 0),
          ifelse((model_time - 1) %% 18 == 0, C_test, 0)
        ) +
          ifelse(
            model_time <= 18 * 2,
            ifelse((model_time - 1) %% 9 == 0, C_imag, 0),
            ifelse((model_time - 1) %% 18 == 0, C_imag, 0)
          ),
        r = dr / (365.24 / 21)
      ),
      cost_ae = 0,
      cost_eol = 0,
      cost_total = cost_treat + cost_other + cost_ae + cost_eol,
      utility = discount(U_pd / (365.24 / 21), r = dr / (365.24 / 21))
    )
    # 定义死亡（D）状态
    state_d <- define_state(
      cost_treat = 0,
      cost_other = 0,
      cost_ae = 0,
      cost_eol = discount(ifelse(state_time == 1, C_eol, 0), r = dr / (365.24 / 21)),
      cost_total = cost_treat + cost_other + cost_ae + cost_eol,
      utility = 0
    )
    # 定义策略
    strat_TG <- define_strategy(
      transition = mat_TG,
      pfs = state_pfs,
      pd = state_pd,
      d = state_d
    )
    strat_CG <- define_strategy(
      transition = mat_CG,
      pfs = state_pfs,
      pd = state_pd,
      d = state_d
    )
    res_mod2<- run_model(
      TG = strat_TG, 
      CG = strat_CG,
      parameters =  define_parameters (
        ############
        C_imu = par1$base[par1$par == "C_imu"],
        C_oxa = par1$base[par1$par == "C_oxa"],
        C_cap = par1$base[par1$par == "C_cap"],
        C_bev = par1$base[par1$par == "C_bev"],
        C_foriB = par1$base[par1$par == "C_foriB"],
        C_foriC = par1$base[par1$par == "C_foriC"],
        C_fori = par1$base[par1$par == "C_fori"],
        C_IriCaB = par1$base[par1$par == "C_IriCaB"],
        C_IriB = par1$base[par1$par == "C_IriB"],
        C_XeB = par1$base[par1$par == "C_XeB"],
        C_CaB = par1$base[par1$par == "C_CaB"],
        C_Fru = par1$base[par1$par == "C_Fru"],
        C_bst = par1$base[par1$par == "C_bst"],
        C_test = par1$base[par1$par == "C_test"],
        C_imag = par1$base[par1$par == "C_imag"],
        C_eol = par1$base[par1$par == "C_eol"],
        C_ane = par1$base[par1$par == "C_ane"],
        C_leuk = par1$base[par1$par == "C_leuk"],
        C_throm = par1$base[par1$par == "C_throm"],
        C_hfsyn = par1$base[par1$par == "C_hfsyn"],
        C_dia = par1$base[par1$par == "C_dia"],
        C_ormu = par1$base[par1$par == "C_ormu"],
        U_ane = par1$base[par1$par == "U_ane"],
        U_leuk = par1$base[par1$par == "U_leuk"],
        U_throm = par1$base[par1$par == "U_throm"],
        U_hfsyn = par1$base[par1$par == "U_hfsyn"],
        U_diar = par1$base[par1$par == "U_diar"],
        U_ormu = par1$base[par1$par == "U_ormu"],
        U_pfs2 = par1$base[par1$par == "U_pfs2"],
        U_pfs1 = par1$base[par1$par == "U_pfs1"],
        U_pd = par1$base[par1$par == "U_pd"],
        dr = par1$base[par1$par == "dr"],
        R_leuk2 = par1$base[par1$par == "R_leuk2"],
        R_leuk1 = par1$base[par1$par == "R_leuk1"],
        R_throm2 = par1$base[par1$par == "R_throm2"],
        R_throm1 = par1$base[par1$par == "R_throm1"],
        R_ane2 = par1$base[par1$par == "R_ane2"],
        R_ane1 = par1$base[par1$par == "R_ane1"],
        R_hfsyn2 = par1$base[par1$par == "R_hfsyn2"],
        R_hfsyn1 = par1$base[par1$par == "R_hfsyn1"],
        R_dia2 = par1$base[par1$par == "R_dia2"],
        R_dia1 = par1$base[par1$par == "R_dia1"],
        R_ormu2 = par1$base[par1$par == "R_ormu2"],
        R_ormu1 = par1$base[par1$par == "R_ormu1"],
        R_foriB2 = par1$base[par1$par == "R_foriB2"],
        R_foriC2 = par1$base[par1$par == "R_foriC2"],
        R_fori2 = par1$base[par1$par == "R_fori2"],
        R_IriCaB2 = par1$base[par1$par == "R_IriCaB2"],
        R_IriB2 = par1$base[par1$par == "R_IriB2"],
        R_XeB2 = par1$base[par1$par == "R_XeB2"],
        R_CaB2 = par1$base[par1$par == "R_CaB2"],
        R_Fru2 = par1$base[par1$par == "R_Fru2"],
        R_bst2 = par1$base[par1$par == "R_bst2"],
        R_foriB1 = par1$base[par1$par == "R_foriB1"],
        R_foriC1 = par1$base[par1$par == "R_foriC1"],
        R_fori1 = par1$base[par1$par == "R_fori1"],
        R_IriCaB1 = par1$base[par1$par == "R_IriCaB1"],
        R_IriB1 = par1$base[par1$par == "R_IriB1"],
        R_XeB1 = par1$base[par1$par == "R_XeB1"],
        R_CaB1 = par1$base[par1$par == "R_CaB1"],
        R_Fru1 = par1$base[par1$par == "R_Fru1"],
        R_bst1 = par1$base[par1$par == "R_bst1"],
        ###############
      ),
      cycles = 347,
      cost = cost_total, effect = utility,
      state_time_limit = c(pd=340,d=1),
      method = "life-table",
      init = c(1,0,0))
    a2 <- summary(res_mod2)
    a2[["res_values"]]$ICER <- (a2[["res_values"]][["cost_total"]][1]-a2[["res_values"]][["cost_total"]][2])/
      (a2[["res_values"]][["utility"]][1]-a2[["res_values"]][["utility"]][2])
    
    results1[[paste0("OS_",names(TG_results)[i], "_PFS_", names(CG_results)[j])]] <- a2[["res_values"]]
  }
}




