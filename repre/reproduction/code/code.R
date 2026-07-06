
# tidyverse  : 数据读取/清洗/汇总（dplyr, readr, purrr, tidyr, ggplot2）
# lme4       : 拟合(广义)线性混合效应模型 —— 原文核心方法
# lmerTest   : 为 lmer 提供 Satterthwaite 自由度 p 值（与原文口径一致）
# broom.mixed: tidy() 将模型结果整理成数据框，便于提取/导出
# gridExtra  : 拼接 ggplot 图（Figure 2A/2B、3A/3B）
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(gridExtra)
})


ROOT       <- "/Users/songshengcamellia/Desktop/R/pre"
STROOP_DIR <- file.path(ROOT, "stroop")           # 每名被试 1 个 Stroop CSV
SWITCH_DIR <- file.path(ROOT, "task_switching")   # 每名被试 2 个 CSV（高/低奖赏）
Q_PATH     <- file.path(ROOT, "questionnaire_data.csv")  # 问卷：nfc/gad/bis-11/bis/bas
OUT        <- file.path(ROOT, "reproduction_outputs")    # 结果输出目录
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)  # 不存在则创建

# 安静读取 CSV：抑制列类型提示信息，自动修复重名列
quiet_read_csv <- function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE,
                                   name_repair = "unique_quiet"))
}


## --- 2.1 Stroop 任务 ---------------------------------------------------------
# 文件名形如 "001_stroop_2016_Nov_11_0958.csv"，开头数字即被试编号。
read_stroop <- function(f) {
  sid <- as.integer(sub("_.*", "", basename(f)))   # 取首个下划线前的数字 = 被试编号
  d   <- quiet_read_csv(f)
  d %>%
    # 行筛选：trials.thisN 非 NA = 正式试次；practice_1.thisN 为 NA = 排除练习阶段
    filter(!is.na(trials.thisN), is.na(practice_1.thisN)) %>%
    transmute(
      subject     = sid,
      trial       = as.integer(trials.thisN),  # 试次序号（离散）
      congruent   = as.integer(congruent),     # 类别：1=一致, 0=不一致
      incongruent = 1L - congruent,            # 类别：1=不一致（论文关注的预测项）
      correct     = as.integer(resp.corr),     # 类别：1=正确, 0=错误(含漏答)
      rt          = as.numeric(resp.rt),        # 连续：反应时(秒)，NA=漏答
      rt_ms       = rt * 1000                   # 连续：反应时(毫秒)，用于报告/建模
    )
}
stroop <- map_dfr(list.files(STROOP_DIR, "\\.csv$", full.names = TRUE), read_stroop)

## --- 2.2 task-switching 任务 -------------------------------------------------
# 文件名 "001_taskSwitch_..."(低奖赏) 或 "001_taskSwitchHR_..."(高奖赏)。
# 文件内 condition 列："HR"=高奖赏(5c)，其它=低奖赏(1c)。
read_switch <- function(f) {
  base    <- basename(f)
  sid     <- as.integer(sub("_.*", "", base))              # 被试编号
  # 提取时间戳(如 2016_Nov_11_1024)，用于恢复被试内 block 先后顺序
  file_ts <- sub("^[0-9]+_(taskSwitchHR|taskSwitch)_(.*)\\.csv$", "\\2", base)
  d       <- quiet_read_csv(f)
  cond    <- na.omit(d$condition)[1]                        # 同文件内 condition 恒定
  d %>%
    # trials.thisN 非 NA = 正式试次；trials_2.thisN 为 NA = 排除练习块
    filter(!is.na(trials.thisN), is.na(trials_2.thisN)) %>%
    transmute(
      subject = sid,
      file_ts = file_ts,
      reward  = if_else(cond == "HR", 1L, 0L),  # 类别：1=高奖赏, 0=低奖赏
      trial   = as.integer(trials.thisN),
      task    = as.integer(spades_or_hearts),    # 类别：当前试次任务类型(黑桃/红心)
      correct = as.integer(resp.corr),           # 1=正确, 0=错误(含漏答)
      rt      = as.numeric(resp.rt),
      rt_ms   = rt * 1000
    )
}
ts_raw <- map_dfr(list.files(SWITCH_DIR, "\\.csv$", full.names = TRUE), read_switch) %>%
  arrange(subject, file_ts, trial)

## --- 2.3 恢复 block 顺序，并判定 switch / repeat ------------------------------
# 每名被试两个 block。按时间戳先后给 block 编号(1=先做,2=后做)，
# 以便把顺序/练习效应作为协变量(block)纳入模型。
block_lookup <- ts_raw %>%
  distinct(subject, file_ts) %>%
  group_by(subject) %>%
  mutate(block = rank(file_ts, ties.method = "first")) %>%
  ungroup()
ts_raw <- ts_raw %>% left_join(block_lookup, by = c("subject", "file_ts"))

# 试次类型：与前一试次任务比较。相同=repeat，不同=switch；
# 每个 block 首个试次无前序任务，记为 NA(不进入分析)。
ts_raw <- ts_raw %>%
  group_by(subject, block) %>%
  arrange(trial, .by_group = TRUE) %>%
  mutate(prev_task  = lag(task),
         trial_type = if_else(is.na(prev_task), NA_character_,
                              if_else(task != prev_task, "switch", "repeat"))) %>%
  ungroup()

## --- 2.4 问卷数据 ------------------------------------------------------------
# 列：subject, nfc(认知需求), gad, bis-11, bis, bas(行为激活)；均为连续变量。
quest <- quiet_read_csv(Q_PATH) %>%
  rename(bis11 = `bis-11`) %>%
  mutate(subject = as.integer(subject))


## 3. 被试纳入与排除（★ 缺 NFC 按需剔除，对齐原文）
have_stroop <- unique(stroop$subject)
have_switch <- unique(ts_raw$subject)
have_nfc    <- quest %>% filter(!is.na(nfc)) %>% pull(subject)

# (a) 基础人群：两任务都有数据（★ 不再要求有 NFC）
candidate <- sort(intersect(have_stroop, have_switch))

# 先把数据限制在候选被试范围内
stroop <- filter(stroop, subject %in% candidate)
ts_raw <- filter(ts_raw, subject %in% candidate)

# (b) 任一任务正确率 < .80 即剔除 → 主样本
stroop_acc <- stroop %>% group_by(subject) %>% summarise(acc = mean(correct, na.rm = TRUE))
ts_acc     <- ts_raw %>% group_by(subject) %>% summarise(acc = mean(correct, na.rm = TRUE))
excl_low_acc <- sort(union(stroop_acc$subject[stroop_acc$acc < 0.80],
                           ts_acc$subject[ts_acc$acc     < 0.80]))

stroop <- filter(stroop, !subject %in% excl_low_acc)
ts_raw <- filter(ts_raw, !subject %in% excl_low_acc)

# (c) 主样本 / NFC 分析样本
main_subjects <- sort(unique(ts_raw$subject))
nfc_subjects  <- sort(intersect(main_subjects, have_nfc))
no_nfc        <- setdiff(main_subjects, nfc_subjects)

cat("=== 样本审计 ===\n")
cat("Stroop 唯一被试数 :", length(have_stroop), "\n")
cat("task-switching 唯一被试数 :", length(have_switch), "\n")
cat("两任务都有数据的候选被试 :", length(candidate), "\n")
cat("因正确率<80%剔除 :", excl_low_acc, "（", length(excl_low_acc), "人）\n")
cat("★ 主样本 N（所有非 NFC 分析）:", length(main_subjects), "人\n")
cat("  其中缺 NFC、仅在 NFC 分析中剔除 :", no_nfc, "\n")
cat("★ NFC 分析样本 N :", length(nfc_subjects), "人\n")
cat(strrep("=", 64), "\n\n")


## ===========================================================================
## 4. 描述性统计（复现论文 Table 1）
stroop_correct <- stroop %>% filter(correct == 1, !is.na(rt_ms))

# Stroop：按一致性分条件
table1_stroop_rt <- stroop_correct %>%
  group_by(subject, congruent) %>% summarise(med = median(rt_ms), .groups = "drop") %>%
  group_by(congruent) %>% summarise(RT_mean = mean(med), RT_sd = sd(med), .groups = "drop")
table1_stroop_acc <- stroop %>%
  group_by(subject, congruent) %>% summarise(a = mean(correct), .groups = "drop") %>%
  group_by(congruent) %>% summarise(ACC_mean = mean(a), ACC_sd = sd(a), .groups = "drop")

# task-switching：按奖赏 × 试次类型分条件
ts_corr <- ts_raw %>% filter(correct == 1, !is.na(rt_ms), !is.na(trial_type))
table1_ts_rt <- ts_corr %>%
  group_by(subject, reward, trial_type) %>% summarise(med = median(rt_ms), .groups = "drop") %>%
  group_by(reward, trial_type) %>% summarise(RT_mean = mean(med), RT_sd = sd(med), .groups = "drop")
table1_ts_acc <- ts_raw %>% filter(!is.na(trial_type)) %>%
  group_by(subject, reward, trial_type) %>% summarise(a = mean(correct), .groups = "drop") %>%
  group_by(reward, trial_type) %>% summarise(ACC_mean = mean(a), ACC_sd = sd(a), .groups = "drop")

cat("=== Table 1：Stroop 描述统计 ===\n")
print(left_join(table1_stroop_rt, table1_stroop_acc, by = "congruent"))
cat("\n=== Table 1：task-switching 描述统计 ===\n")
print(left_join(table1_ts_rt, table1_ts_acc, by = c("reward", "trial_type")))

# 导出描述统计
desc_out <- bind_rows(
  left_join(table1_stroop_rt, table1_stroop_acc, by = "congruent") %>%
    mutate(block_label = if_else(congruent == 1, "Stroop-Congruent", "Stroop-Incongruent")) %>%
    select(condition = block_label, RT_mean, RT_sd, ACC_mean, ACC_sd),
  left_join(table1_ts_rt, table1_ts_acc, by = c("reward", "trial_type")) %>%
    mutate(condition = paste0("TS-", if_else(reward == 1, "High", "Low"), "-",
                              str_to_title(trial_type))) %>%
    select(condition, RT_mean, RT_sd, ACC_mean, ACC_sd)
)
write_csv(desc_out, file.path(OUT, "descriptive_stats.csv"))


## ===========================================================================
## 5. Stroop 干扰效应（论文 3.1 节）
## ---------------------------------------------------------------------------
## RT：log(RT_ms) ~ 不一致 + (1|被试)
## 正确率：logistic 混合模型 correct ~ 不一致 + (1|被试)
## ===========================================================================
m_stroop_rt  <- lmer(log(rt_ms) ~ incongruent + (1 | subject), data = stroop_correct)
m_stroop_acc <- glmer(correct ~ incongruent + (1 | subject), data = stroop,
                      family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("\n=== Stroop 干扰效应 ===\n")
cat("RT(log)  β =", round(fixef(m_stroop_rt)["incongruent"], 4), "\n")
cat("正确率(logit) β =", round(fixef(m_stroop_acc)["incongruent"], 3), "\n")


## ===========================================================================
## 6. 构造 EF 指标（Stroop RT cost）
ef_cost <- stroop_correct %>%
  group_by(subject) %>%
  mutate(zrt = scale(rt_ms)[, 1]) %>%
  summarise(stroop_cost = mean(zrt[congruent == 0]) - mean(zrt[congruent == 1]),
            .groups = "drop")


## ===========================================================================
## 7. switch cost 主效应（论文 3.1 节）
## ===========================================================================
ts_model_data <- ts_corr %>%
  mutate(switch  = as.integer(trial_type == "switch"),  # 重复=0, 切换=1
         log_rt  = log(rt_ms),
         block_c = block - mean(block))                  # block 中心化(线性协变量)
m_switch_rt <- lmer(log_rt ~ switch + (1 | subject), data = ts_model_data)

ts_acc_data <- ts_raw %>% filter(!is.na(trial_type)) %>%
  mutate(switch = as.integer(trial_type == "switch"))
m_switch_acc <- glmer(correct ~ switch + (1 | subject), data = ts_acc_data,
                      family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("\n=== switch cost 主效应 ===\n")
cat("RT(log)  β =", round(fixef(m_switch_rt)["switch"], 4), "\n")
cat("正确率(logit) β =", round(fixef(m_switch_acc)["switch"], 3), "\n")


## ===========================================================================
## 8. 组水平奖赏效应（论文 3.2 节）
## ---------------------------------------------------------------------------
## 不含个体差异变量时，奖赏对 switch cost 的影响(switch×reward)。
## ===========================================================================
m_group_reward <- lmer(log_rt ~ switch * reward + block_c + (1 | subject),
                       data = ts_model_data %>% mutate(reward = ts_corr$reward))
cat("\n=== 组水平 switch×reward 交互 ===\n")
print(round(summary(m_group_reward)$coefficients["switch:reward", ], 4))


## ===========================================================================
## 9. 相关分析（论文 3.1 节）
sub_rt <- ts_corr %>%
  group_by(subject, trial_type) %>% summarise(med = median(rt_ms), .groups = "drop") %>%
  pivot_wider(names_from = trial_type, values_from = med) %>%
  mutate(switch_rt_cost = switch - `repeat`)             # RT 切换代价 = switch − repeat
sub_acc <- ts_raw %>% filter(!is.na(trial_type)) %>%
  group_by(subject, trial_type) %>% summarise(a = mean(correct), .groups = "drop") %>%
  pivot_wider(names_from = trial_type, values_from = a) %>%
  mutate(switch_acc_cost = `repeat` - switch)            # 正确率切换代价 = repeat − switch

subj_tab <- ef_cost %>%
  left_join(sub_rt  %>% select(subject, switch_rt_cost),  by = "subject") %>%
  left_join(sub_acc %>% select(subject, switch_acc_cost), by = "subject") %>%
  left_join(quest,  by = "subject")

cat("\n=== 相关 ===\n")
cat(sprintf("r(EF, switch RT cost)  = %.4f, p = %.3f\n",
            cor(subj_tab$stroop_cost, subj_tab$switch_rt_cost, use = "complete.obs"),
            cor.test(subj_tab$stroop_cost, subj_tab$switch_rt_cost)$p.value))
cat(sprintf("r(EF, switch ACC cost) = %.5f, p = %.3f\n",
            cor(subj_tab$stroop_cost, subj_tab$switch_acc_cost, use = "complete.obs"),
            cor.test(subj_tab$stroop_cost, subj_tab$switch_acc_cost)$p.value))
cat(sprintf("r(EF, NFC)             = %.3f, p = %.3f\n",
            cor(subj_tab$stroop_cost, subj_tab$nfc, use = "complete.obs"),
            with(na.omit(subj_tab[, c("stroop_cost", "nfc")]),
                 cor.test(stroop_cost, nfc)$p.value)))


## ===========================================================================
## 10. 试次级建模数据框（用于 Table 2 / Table 3）
## ---------------------------------------------------------------------------
## 因变量：正确试次的 log(RT_ms)；预测变量：switch、reward、block(线性)、
## 以及标准化的个体差异变量(EF_z 或 NFC_z)；随机效应：随机截距 (1|subject)。
## ===========================================================================
ana <- ts_corr %>%
  mutate(switch  = as.integer(trial_type == "switch"),
         log_rt  = log(rt_ms),
         block_c = block - mean(block)) %>%
  left_join(ef_cost, by = "subject") %>%
  left_join(quest %>% select(subject, nfc, bas), by = "subject") %>%
  # EF 在主样本内标准化；NFC / BAS 用 na.rm 标准化，使缺 NFC 的被试得到 NA，
  # 从而在 Table 3 / NFC 相关 / BAS 控制模型中被自动排除（→ NFC 分析样本 51 人）。
  mutate(EF_z  = as.numeric(scale(stroop_cost)),
         NFC_z = (nfc - mean(nfc, na.rm = TRUE)) / sd(nfc, na.rm = TRUE),
         BAS_z = (bas - mean(bas, na.rm = TRUE)) / sd(bas, na.rm = TRUE))

## --- Table 2：EF 模型 ---
m_tab2 <- lmer(log_rt ~ switch * reward * EF_z + block_c + (1 | subject), data = ana)
cat("\n=== Table 2：log RT ~ switch × reward × EF + block ===\n")
print(round(summary(m_tab2)$coefficients, 4))

## --- Table 3：NFC 模型 ---
m_tab3 <- lmer(log_rt ~ switch * reward * NFC_z + block_c + (1 | subject),
               data = filter(ana, !is.na(NFC_z)))
cat("\n=== Table 3：log RT ~ switch × reward × NFC + block ===\n")
print(round(summary(m_tab3)$coefficients, 4))


## ===========================================================================
## 11. 正确率三重交互模型（论文 3.3 节）
## ===========================================================================
acc_ana <- ts_raw %>% filter(!is.na(trial_type)) %>%
  mutate(switch  = as.integer(trial_type == "switch"),
         block_c = block - mean(block)) %>%
  left_join(ef_cost, by = "subject") %>%
  mutate(EF_z = as.numeric(scale(stroop_cost)))
m_acc3way <- glmer(correct ~ switch * reward * EF_z + block_c + (1 | subject),
                   family = binomial, data = acc_ana,
                   control = glmerControl(optimizer = "bobyqa"))
cat("\n=== 正确率三重交互（switch×reward×EF）===\n")
print(round(summary(m_acc3way)$coefficients["switch:reward:EF_z", ], 4))


## ===========================================================================
## 12. BAS 控制模型（论文 3.5 节）
## ---------------------------------------------------------------------------
## 加入奖赏敏感性(BAS)及其与 switch/reward 的交互后，核心三重交互是否仍稳健。
## ===========================================================================
m_tab2_bas <- lmer(log_rt ~ switch * reward * EF_z + block_c +
                     switch * reward * BAS_z + (1 | subject), data = ana)
m_tab3_bas <- lmer(log_rt ~ switch * reward * NFC_z + block_c +
                     switch * reward * BAS_z + (1 | subject),
                   data = filter(ana, !is.na(NFC_z), !is.na(BAS_z)))
cat("\n=== 控制 BAS 后的核心三重交互 ===\n")
cat("Stroop(EF) 模型 switch×reward×EF :\n")
print(round(coef(summary(m_tab2_bas))["switch:reward:EF_z", ], 4))
cat("NFC 模型 switch×reward×NFC :\n")
print(round(coef(summary(m_tab3_bas))["switch:reward:NFC_z", ], 4))


## ===========================================================================
## 13. Figure 2 / Figure 3 复现（中位数分组柱状图 + 散点回归）
## ===========================================================================
# 每名被试在高/低奖赏 block 下的 switch cost(ms)
sub_block_cost <- ts_corr %>%
  group_by(subject, reward, trial_type) %>% summarise(med = median(rt_ms), .groups = "drop") %>%
  pivot_wider(names_from = trial_type, values_from = med) %>%
  mutate(sw_cost_rt = switch - `repeat`) %>%
  left_join(subj_tab %>% select(subject, stroop_cost, nfc), by = "subject")

# 配色（沿用本项目粉色方案，不照搬原文灰/白）
REW_FILL <- c("high reward block" = "#ED3CB0", "low reward block" = "#F7B9DD")

## --- Figure 2A：EF 中位数分组 + 全体被试（与原文一样含 all subjects 组）---
# 低 EF = 高 Stroop cost（>= 中位数）；高 EF = 低 Stroop cost。
ef_split <- sub_block_cost %>%
  mutate(grp = if_else(stroop_cost >= median(stroop_cost),
                       "low EF\n(high Stroop cost)", "high EF\n(low Stroop cost)"))
ef_all <- sub_block_cost %>% mutate(grp = "all subjects")   # 全体被试组
fig2_bar <- bind_rows(ef_split, ef_all) %>%
  mutate(grp = factor(grp, levels = c("low EF\n(high Stroop cost)",
                                      "high EF\n(low Stroop cost)", "all subjects")),
         reward = factor(reward, levels = c(1, 0),
                         labels = c("high reward block", "low reward block"))) %>%
  group_by(grp, reward) %>%
  summarise(m = mean(sw_cost_rt), se = sd(sw_cost_rt) / sqrt(n()), .groups = "drop")
p2a <- ggplot(fig2_bar, aes(grp, m, fill = reward)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = REW_FILL) +
  labs(x = NULL, y = "task switch cost (ms)", fill = NULL, title = "Fig. 2A") +
  theme_bw() + theme(legend.position = "top")
reward_effect_ef <- sub_block_cost %>%
  pivot_wider(id_cols = c(subject, stroop_cost), names_from = reward, values_from = sw_cost_rt) %>%
  mutate(reward_effect = `1` - `0`)            # 奖赏诱发的 switch cost 变化
p2b <- ggplot(reward_effect_ef, aes(stroop_cost, reward_effect)) +
  geom_point(color = "#ED3CB0") + geom_smooth(method = "lm", se = TRUE, color = "#2E5BFF") +
  labs(x = "EF capacity (Stroop RT cost, z)",
       y = "Reward-induced change in switch cost (ms)", title = "Fig. 2B") + theme_bw()

## --- Figure 3A：NFC 中位数分组 + 全体被试 ---
nfc_split <- sub_block_cost %>% filter(!is.na(nfc)) %>%
  mutate(grp = if_else(nfc <= median(nfc), "low\nneed for cognition", "high\nneed for cognition"))
nfc_all <- sub_block_cost %>% mutate(grp = "all subjects")
fig3_bar <- bind_rows(nfc_split, nfc_all) %>%
  mutate(grp = factor(grp, levels = c("low\nneed for cognition",
                                      "high\nneed for cognition", "all subjects")),
         reward = factor(reward, levels = c(1, 0),
                         labels = c("high reward block", "low reward block"))) %>%
  group_by(grp, reward) %>%
  summarise(m = mean(sw_cost_rt), se = sd(sw_cost_rt) / sqrt(n()), .groups = "drop")
p3a <- ggplot(fig3_bar, aes(grp, m, fill = reward)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = REW_FILL) +
  labs(x = NULL, y = "task switch cost (ms)", fill = NULL, title = "Fig. 3A") +
  theme_bw() + theme(legend.position = "top")
reward_effect_nfc <- sub_block_cost %>% filter(!is.na(nfc)) %>%
  pivot_wider(id_cols = c(subject, nfc), names_from = reward, values_from = sw_cost_rt) %>%
  mutate(reward_effect = `1` - `0`)
p3b <- ggplot(reward_effect_nfc, aes(nfc, reward_effect)) +
  geom_point(color = "#ED3CB0") + geom_smooth(method = "lm", se = TRUE, color = "#2E5BFF") +
  labs(x = "Need for Cognition score",
       y = "Reward-induced change in switch cost (ms)", title = "Fig. 3B") + theme_bw()


## ===========================================================================
## 14. 保存图形、被试层数据
ggsave(file.path(OUT, "figure2.png"),
       grid.arrange(p2a, p2b, ncol = 2), width = 10, height = 4)
ggsave(file.path(OUT, "figure3.png"),
       grid.arrange(p3a, p3b, ncol = 2), width = 10, height = 4)
write_csv(subj_tab, file.path(OUT, "subject_level_data.csv"))
