# 加载所需的库
library(writexl)
library(carcass)
library(survminer)
library(ggplot2)
library(dplyr)
library(heemod)
library(openxlsx)
library(cowplot)
library(gridExtra)
library(extrafont)
library(readxl)
library(showtext)
library(sysfonts)
library(scales)  
library(xml2)

# 基础分析 ##########################################################
setwd("C:/PHD/CRC/外推20年基线分析")  # 设置工作目录

# 读取治疗组和对照组的概率数据
TG <- read.xlsx("TG_prob.xlsx")
CG <- read.xlsx("CG_prob.xlsx")

# Pptp 表示 PD 不发生状态转移的概率
TG$pPTP[1] <- 1
CG$pPTP[1] <- 1

# 读取参数表
par <- read.xlsx("参数1.xlsx")

# 定义治疗组的状态转移矩阵
mat_TG <- define_transition(
  state_names = c("pfs", "pd", "d"),
  TG[model_time, "pFTF"], TG[model_time, "pFTP"], TG[model_time, "pFTD"],
  0, TG[model_time, "pPTP"], TG[model_time, "pPTD"],
  0, 0, 1
)

# 定义对照组的状态转移矩阵
mat_CG <- define_transition(
  state_names = c("pfs", "pd", "d"),
  CG[model_time, "pFTF"], CG[model_time, "pFTP"], CG[model_time, "pFTD"],
  0, CG[model_time, "pPTP"], CG[model_time, "pPTD"],
  0, 0, 1
)

# 定义模型中的状态和成本计算逻辑
state_pfs <- define_state(
  cost_treat = discount(
    dispatch_strategy(
      TG = ifelse(model_time <= 6, C_imu + C_oxa + C_cap + C_bev, C_cap + C_bev),
      CG = ifelse(model_time <= 6, C_oxa + C_cap + C_bev, C_cap + C_bev)
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
  utility = discount(
    ifelse(
      model_time == 1,
      dispatch_strategy(
        TG = (U_pfs2 - (U_ane * R_ane2 + U_leuk * R_leuk2 + U_throm * R_throm2 + U_hfsyn * R_hfsyn2 + U_diar * R_dia2 + U_ormu * R_ormu2)) / (365.24 / 21),
        CG = (U_pfs1 - (U_ane * R_ane1 + U_leuk * R_leuk1 + U_throm * R_throm1 + U_hfsyn * R_hfsyn1 + U_diar * R_dia1 + U_ormu * R_ormu1)) / (365.24 / 21)
      ),
      ifelse(
        TRUE,
        U_pfs2 / (365.24 / 21),
        U_pfs1 / (365.24 / 21)
      )
    ),
    r = dr / (365.24 / 21)
  )
)

state_pd <- define_state(
  cost_treat = discount(
    dispatch_strategy(
      TG = C_bst * R_bst2 +
        ifelse((state_time - 1) %% 2 == 0, C_foriB * R_foriB2 + C_foriC * R_foriC2 + C_IriB * R_IriCaB2 + C_XeB * R_XeB2, 0) +
        ifelse((state_time) %% 4 == 0, 0, C_Fru * R_Fru2) +
        C_fori * R_fori2 + C_IriCaB * R_IriCaB2 + C_CaB * R_CaB2,
      CG = C_bst * R_bst1 +
        ifelse((state_time - 1) %% 2 == 0, C_foriB * R_foriB1 + C_foriC * R_foriC1 + C_IriB * R_IriCaB1 + C_XeB * R_XeB1, 0) +
        ifelse((state_time) %% 4 == 0, 0, C_Fru * R_Fru1) +
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

# 运行模型
res_mod2 <- run_model(
  TG = strat_TG, 
  CG = strat_CG,
  parameters = define_parameters(
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
  ),
  cycles = 347,
  cost = cost_total,
  effect = utility,
  state_time_limit = c(pd = 340, d = 1),
  method = "life-table",
  init = c(1, 0, 0)
)

# 输出基础分析结果
summary(res_mod2)

par_sentences1 <- apply(par, 1, function(row) {
  paste0(row["par"], "=par1$base[par1$par ==",row["par"],"]")
})
par_sentences1_df <- data.frame(Formula = unlist(par_sentences1))
write.xlsx(par_sentences1_df,file = "par_sentences1_df1.xlsx")
a2 <- summary(res_mod2)
a2[["res_values"]]$ICER[1]=(a2[["res_values"]]$cost_total[1]-a2[["res_values"]]$cost_total[2])/
  (a2[["res_values"]]$utility[1]-a2[["res_values"]]$utility[2])
write.xlsx(a2[["res_values"]], file="基础分析.xlsx")
a2 <- summary(res_mod2)
human1 <- get_counts(res_mod2)
cost_df2 <- get_values(res_mod2)
write.xlsx(human1, file="human.xlsx")
write.xlsx(cost_df2, file="cost_df.xlsx")


# 单因素敏感性分析 #################################################
setwd("C:/PHD/CRC/外推20年基线分析")
par <- read.xlsx("参数1.xlsx",sheet = "par")

#17:53
# 定义一个函数来修正置信区间的上下限
fix_confidence_interval <- function(low, high) {
  low <- max(0, min(1, low))
  high <- max(0, min(1, high))
  return(list(low = low, high = high))
}

for (i in 1:nrow(par)) {
  # 如果par列首字母不是C
  if (par$distribution[i] == "beta") {
    if (par$base[i]!= 0) {
      # 修正置信区间的上下限
      conf_interval <- fix_confidence_interval(par$low[i], par$hight[i])
      par$a[i] <- shapeparameter(m = par$base[i], lwr = conf_interval$low, upr = conf_interval$high)$a
      par$b[i] <- shapeparameter(m = par$base[i], lwr = conf_interval$low, upr = conf_interval$high)$b
    }
  } else {
    par$sd[i] <- (par$hight[i] - par$low[i]) / (1.96 * 2)
  }
}

par$dis <- apply(par, 1, function(row) {
  paste0("D_",row["par"])
})

par_sentences1 <- apply(par, 1, function(row) {
  paste0(row["par"], "~", row["distribution"],"(","mean=", row["base"], ", ","sd=", row["sd"], "),")
})

par_sentences2 <- apply(par, 1, function(row) {
  paste0(row["par"], "~beta(","shape1=", row["a"], ", ","shape2=", row["b"], "),")
})
par_sentences1_df <- data.frame(Formula = unlist(par_sentences1))
par_sentences2_df <- data.frame(Formula = unlist(par_sentences2))

write.xlsx(par_sentences1_df,file = "par_sentences1_df.xlsx")
write.xlsx(par_sentences2_df,file = "par_sentences2_df.xlsx")

write.xlsx(par,file = "par3.xlsx")


########概率敏感性分析########################################
rsp <- define_psa(
  ###########
  C_imu~gamma(mean=6850.00000, sd=873.724490),
  C_oxa~gamma(mean= 190.23830, sd= 24.265089),
  C_cap~gamma(mean=  44.35000, sd=  5.656888),
  C_bev~gamma(mean=1293.72980, sd=165.016556),
  C_foriB~gamma(mean= 982.31543, sd=125.295336),
  C_foriC~gamma(mean=1765.11643, sd=225.142401),
  C_fori~gamma(mean= 119.82890, sd= 15.284298),
  C_IriCaB~gamma(mean= 908.08018, sd=115.826554),
  C_IriB~gamma(mean= 882.11096, sd=112.514154),
  C_XeB~gamma(mean=1078.97177, sd=137.623950),
  C_CaB~gamma(mean= 888.45575, sd=113.323438),
  C_Fru~gamma(mean=1034.63892, sd=131.969250),
  C_bst~gamma(mean= 293.77264, sd= 37.471000),
  C_test~gamma(mean=1531.83800, sd=195.387500),
  C_imag~gamma(mean=  42.30562, sd=  5.396125),
  C_eol~gamma(mean=4469.47659, sd=570.086300),
  C_ane~gamma(mean=1089.50520, sd=138.967500),
  C_leuk~gamma(mean=2309.39156, sd=294.565250),
  C_throm~gamma(mean=3809.85329, sd=485.950675),
  C_hfsyn~gamma(mean=2042.33588, sd=260.502025),
  C_dia~gamma(mean=1005.97647, sd=128.313325),
  C_ormu~gamma(mean=1481.46776, sd=188.962725),
  
  U_ane~beta(shape1=40.581667, shape2= 436.849706),
  U_leuk~beta(shape1=40.354444, shape2= 408.028272),
  U_throm~beta(shape1=42.763000, shape2=1112.993757),
  U_hfsyn~beta(shape1=43.899111, shape2=3614.360148),
  U_diar~beta(shape1=40.354444, shape2= 408.028272),
  U_ormu~beta(shape1=41.717778, shape2= 653.578519),
  U_pfs2~beta(shape1= 9.721559, shape2=   1.715569),
  U_pfs1~beta(shape1= 9.778512, shape2=   2.444628),
  U_pd~beta(shape1=11.270000, shape2=   4.168356),
  dr~beta(shape1= 5.887500, shape2= 111.862500),
  R_leuk2~beta(shape1=60.100000, shape2= 941.566667),
  R_leuk1~beta(shape1=60.815000, shape2=1180.307449),
  R_throm2~beta(shape1=62.700000, shape2=3072.300000),
  R_throm1~beta(shape1=61.465000, shape2=1514.560641),
  R_ane2~beta(shape1=62.050000, shape2=2006.283333),
  R_ane1~beta(shape1=62.700000, shape2=3072.300000),
  R_hfsyn2~beta(shape1=61.400000, shape2=1473.600000),
  R_hfsyn1~beta(shape1=61.465000, shape2=1514.560641),
  R_dia2~beta(shape1=62.700000, shape2=3072.300000),
  R_dia1~beta(shape1=60.815000, shape2=1180.307449),
  R_ormu2~beta(shape1=63.350000, shape2=6271.650000),
  R_ormu1~beta(shape1=62.700000, shape2=3072.300000),
  R_foriB2~beta(shape1=50.350000, shape2= 189.411905),
  R_foriC2~beta(shape1=62.050000, shape2=2006.283333),
  R_fori2~beta(shape1=63.350000, shape2=6271.650000),
  R_IriCaB2~beta(shape1=63.350000, shape2=6271.650000),
  R_IriB2~beta(shape1=62.700000, shape2=3072.300000),
  R_XeB2~beta(shape1=63.350000, shape2=6271.650000),
  R_CaB2~beta(shape1=62.050000, shape2=2006.283333),
  R_Fru2~beta(shape1=62.700000, shape2=3072.300000),
  R_bst2~beta(shape1=56.850000, shape2= 459.968182),
  R_foriB1~beta(shape1=45.540000, shape2= 114.812113),
  R_foriC1~beta(shape1=62.700000, shape2=3072.300000),
  R_fori1~beta(shape1=62.115000, shape2=2079.781552),
  R_IriCaB1~beta(shape1=62.115000, shape2=2079.781552),
  R_IriB1~beta(shape1=62.700000, shape2=3072.300000),
  R_XeB1~beta(shape1=62.700000, shape2=3072.300000),
  R_CaB1~beta(shape1=62.700000, shape2=3072.300000),
  R_Fru1~beta(shape1=62.700000, shape2=3072.300000),
  R_bst1~beta(shape1=60.425000, shape2=1038.211364)
  
  
)

# 设置总的 PSA 运行次数
num_psa <- 1000 

# 初始化一个空的数据框用于存储最终结果
final_results <- data.frame(Incremental_effect = numeric(0), 
                            Incremental_cost = numeric(0), 
                            ICER = numeric(0))
# 循环执行 PSA
for (i in 1:num_psa) {
  # 逐个执行 PSA
  ndt1 <- run_psa(res_mod2, psa = rsp, N = 1)
  
  # 提取所需的结果
  a <- ndt1$psa[1:(nrow(ndt1$psa) / 2), ]
  b <- ndt1$psa[(nrow(ndt1$psa) / 2 + 1):nrow(ndt1$psa), ]
  
  tab2 <- data.frame(
    Incremental_effect = (a$.effect - b$.effect),
    Incremental_cost = (a$.cost - b$.cost)
  )
  tab2$ICER <- tab2$Incremental_cost / tab2$Incremental_effect
  
  final_results <- rbind(final_results, tab2)
  
  rm(ndt1, tab2, a, b)
  gc()
}
mean(final_results$ICER)
write_xlsx(final_results, path = "散点图数据.xlsx")


######单因素敏感性分析########################
par <- read.xlsx("参数1.xlsx",sheet = "par")
a1 <- apply(par, 1, function(row) {
  paste0(row["par"],  ", ", row["low"], ",",row["hight"],",")
})
a1_df <- data.frame(Formula = unlist(a1))
write.xlsx(a1_df,file = "a1_df.xlsx")
dsa <- define_dsa(
  
  C_imu, 5137.50000,8562.50000,
  C_oxa,  142.67872, 237.79787,
  C_cap,   33.26250,  55.43750,
  C_bev,  970.29735,1617.16225,
  C_foriB,  736.73658,1227.89429,
  C_foriC, 1323.83732,2206.39553,
  C_fori,   89.87167, 149.78612,
  C_IriCaB,  681.06013,1135.10022,
  C_IriB,  661.58322,1102.63870,
  C_XeB,  809.22883,1348.71471,
  C_CaB,  666.34181,1110.56969,
  C_Fru,  775.97919,1293.29865,
  C_bst,  220.32948, 367.21580,
  C_test, 1148.87850,1914.79750,
  C_imag,   31.72921,  52.88202,
  C_eol, 3352.10744,5586.84574,
  C_ane,  817.12890,1361.88150,
  C_leuk, 1732.04367,2886.73945,
  C_throm, 2857.38997,4762.31661,
  C_hfsyn, 1531.75191,2552.91984,
  C_dia,  754.48235,1257.47059,
  C_ormu, 1111.10082,1851.83470,
  U_ane,    0.05950,   0.11050,
  U_leuk,    0.06300,   0.11700,
  U_throm,    0.02590,   0.04810,
  U_hfsyn,    0.00840,   0.01560,
  U_diar,    0.06300,   0.11700,
  U_ormu,    0.04200,   0.07800,
  U_pfs2,    0.59500,   1.10500,
  U_pfs1,    0.56000,   1.04000,
  U_pd,    0.51100,   0.94900,
  dr,    0.00000,   0.08000,
  R_leuk2,    0.04500,   0.07500,
  R_leuk1,    0.03675,   0.06125,
  R_throm2,    0.01500,   0.02500,
  R_throm1,    0.02925,   0.04875,
  R_ane2,    0.02250,   0.03750,
  R_ane1,    0.01500,   0.02500,
  R_hfsyn2,    0.03000,   0.05000,
  R_hfsyn1,    0.02925,   0.04875,
  R_dia2,    0.01500,   0.02500,
  R_dia1,    0.03675,   0.06125,
  R_ormu2,    0.00750,   0.01250,
  R_ormu1,    0.01500,   0.02500,
  R_foriB2,    0.15750,   0.26250,
  R_foriC2,    0.02250,   0.03750,
  R_fori2,    0.00750,   0.01250,
  R_IriCaB2,    0.00750,   0.01250,
  R_IriB2,    0.01500,   0.02500,
  R_XeB2,    0.00750,   0.01250,
  R_CaB2,    0.02250,   0.03750,
  R_Fru2,    0.01500,   0.02500,
  R_bst2,    0.08250,   0.13750,
  R_foriB1,    0.21300,   0.35500,
  R_foriC1,    0.01500,   0.02500,
  R_fori1,    0.02175,   0.03625,
  R_IriCaB1,    0.02175,   0.03625,
  R_IriB1,    0.01500,   0.02500,
  R_XeB1,    0.01500,   0.02500,
  R_CaB1,    0.01500,   0.02500,
  R_Fru1,    0.01500,   0.02500,
  R_bst1,    0.04125,   0.06875,
  
)
ndsa1 <- run_dsa(
  res_mod2,
  dsa
)
ev <- (res_mod2$run_model$.cost[1]-res_mod2$run_model$.cost[2])/(res_mod2$run_model$.effect[1]-res_mod2$run_model$.effect[2])
df2 <- subset(ndsa1$dsa, .strategy_names == "TG")
df0 <- subset(ndsa1$dsa, .strategy_names == "CG")
df20 <- data.frame(par=df2$.par_names,color=df2$.par_value, ICER=(df2$.cost-df0$.cost)/(df2$.effect-df0$.effect))
df20 <- df20 %>%
  group_by(par) %>%
  mutate(
    color = as.numeric(color),  # 将color转换为数值类型
    label = ifelse(color == min(color), "low", "high")
  )
df20 <- df20 %>%
  mutate(ICER_adjusted = ICER - ev)
df20 <- df20 %>%
  mutate(color = ifelse(label == "high", "#F8766D", "#00BFC4"))
df_sorted20 <- df20 %>%
  group_by(par) %>%
  mutate(diff = abs(ICER - lag(ICER))) %>%
  filter(!is.na(diff)) %>%
  arrange(diff) %>%
  ungroup()
df_sorted20 <- df20 %>%
  inner_join(df_sorted20 %>% select(par, diff), by = "par") %>%
  arrange(diff, par)

write.xlsx(df_sorted20,file="龙卷风数据.xlsx")
