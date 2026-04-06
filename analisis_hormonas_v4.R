# Title: R scripts for steroid data processing and analysis
# Author: Irene Acosta
# Date: April 2026
# Description: Data cleaning, processing, and statistical analysis of steroid metabolites
# Related manuscript: Childbirth Experience and Maternal Steroid Profiles in the Early Postpartum: An Integrative Analysis of Associations with Maternal-Infant Bonding, Early Infant Neurodevelopment and Perinatal Depression



# Slow sequential processing was used to improve stability and avoid crashes


# Load libreries ----
library(xlsx) 
library(dplyr) 
library(stringr)
library(ggplot2)
library(kableExtra)
library(janitor)
library(car)
library(MASS)
library(sandwich)
library(tidyr)
library(robustbase)
library(ggpubr)


setwd("X") # Set working directory

if (!dir.exists("results")) {
  dir.create("results")
}

# Read bd
env1 <- new.env()
env2 <- new.env()

load("all_oct25.RData", envir = env1)
load("Nuevos.RData", envir = env2)

df1 <- env1$dbs
df2 <- env2$dbs

# Verifications ----
  identical(colnames(df1), colnames(df2))
  dif_col_db1<- setdiff(colnames(df1), colnames(df2))
  dif_col_db2<- setdiff(colnames(df2), colnames(df1))

#  Fix df1 and unify .x and .y
# 1.-Remove  empty columns
df1_s_na <- df1 %>% 
  dplyr::select(where(~ !all(is.na(.))))

# 2.-Detect .x columns and their base names
cols_x <- grep("\\.x$", names(df1_s_na), value = TRUE)
cols_base <- sub("\\.x$", "", cols_x)

# 3.-Combine .x and .y columns
for (col_base in cols_base) {
  
  col_x <- paste0(col_base, ".x")
  col_y <- paste0(col_base, ".y")
  
  if (col_base %in% names(df1_s_na)) {
    message("The column '", col_base, "' already exists. It will not be overwritten.")
  } else {
    values_y <- if (col_y %in% names(df1_s_na)) df1_s_na[[col_y]] else NULL
    
    # Convert to character so coalesce does not fail
    vec_x <- as.character(df1_s_na[[col_x]])
    vec_y <- if (!is.null(values_y)) as.character(values_y) else NULL

    if (!is.null(vec_y)) {
      df1_s_na[[col_base]] <- dplyr::coalesce(vec_x, vec_y)
    } else {
      df1_s_na[[col_base]] <- vec_x
    }
  }
}

# 4.-Remove .x and .y columns
df1_s_na <- df1_s_na %>% dplyr::select(-matches("\\.x$|\\.y$"))
# View(df1_s_na)

# Fix df2 and unify .x and .y
# 1.-Remove  empty columns
df2_s_na <- df2 %>% 
  dplyr::select(where(~ !all(is.na(.))))

# 2.-Detect .x columns and their base names
cols_x <- grep("\\.x$", names(df2_s_na), value = TRUE)
cols_base <- sub("\\.x$", "", cols_x)

# 3.-Combine .x and .y columns
for (col_base in cols_base) {
  
  col_x <- paste0(col_base, ".x")
  col_y <- paste0(col_base, ".y")
  
  if (col_base %in% names(df2_s_na)) {
    message("The column '", col_base, "' already exists. It will not be overwritten.")
  } else {
    values_y <- if (col_y %in% names(df2_s_na)) df2_s_na[[col_y]] else NULL
    
    # Convert to character so coalesce does not fail
    vec_x <- as.character(df2_s_na[[col_x]])
    vec_y <- if (!is.null(values_y)) as.character(values_y) else NULL
    
    if (!is.null(vec_y)) {
      df2_s_na[[col_base]] <- dplyr::coalesce(vec_x, vec_y)
    } else {
      df2_s_na[[col_base]] <- vec_x
    }
  }
}

# 4.-Remove .x and .y columns
df2_s_na <- df2_s_na %>% dplyr::select(-matches("\\.x$|\\.y$"))
# View(df2_s_na)


# Merge the two datasets; remove duplicate columns created by the join (.x, .y) and keep only one

df_union <- dplyr::full_join(df1_s_na, df2_s_na, by = "code_id", suffix = c(".x", ".y"))

df_union_s_na <- df_union %>% 
  dplyr::select(where(~ !all(is.na(.))))

cols_x <- grep("\\.x$", names(df_union_s_na), value = TRUE)
cols_base <- sub("\\.x$", "", cols_x)

for (col_base in cols_base) {
  
  col_x <- paste0(col_base, ".x")
  col_y <- paste0(col_base, ".y")

  if (col_base %in% names(df_union_s_na)) {
    message("The column '", col_base, "' already exists. It will not be overwritten.")
  } else {
    values_y <- if (col_y %in% names(df_union_s_na)) df_union_s_na[[col_y]] else NULL
    
    vec_x <- as.character(df_union_s_na[[col_x]])
    vec_y <- if (!is.null(values_y)) as.character(values_y) else NULL
    
    if (!is.null(vec_y)) {
      df_union_s_na[[col_base]] <- dplyr::coalesce(vec_x, vec_y)
    } else {
      df_union_s_na[[col_base]] <- vec_x
    }
  }
}

df_union_s_na <- df_union_s_na %>% dplyr::select(-matches("\\.x$|\\.y$"))
# View(df_union_s_na)


# Remove duplicates code_id
df_union_s_na <- df_union_s_na %>%
  dplyr::distinct(code_id, .keep_all = TRUE)


db_horm <- read.xlsx2(file = "DEPPOST-BTB_SAL_BBDD-1.xlsx",sheetName = 1)
colnames(db_horm) 
# View(db_horm)

# Rename Excel fields, drop NA values and remove rows where code_id is empty
db_horm_normalitzada <- db_horm %>% 
  rename(code_id = Sample_ID) %>% 
  filter(!(is.na(code_id) | code_id == "")) %>% 
  arrange(code_id)
# View(db_horm_normalitzada)


# Unify code_id values in the clinic database
dbs_normalitzada <- df_union_s_na %>%  
  # Find values starting with C_ or P_ followed by a single digit (1-9) and prepend a 0
  mutate(
    code_id = str_replace(code_id, "C_([1-9])$", "C_0\\1"),
    code_id = str_replace(code_id, "P_([1-9])$", "P_0\\1") 
  ) %>% 
  arrange(code_id)

# View(dbs_normalitzada)


db_join <- dbs_normalitzada %>% left_join(db_horm_normalitzada,by = "code_id") %>% collect()
file <- db_join %>% collect() %>% mutate(across(where(is.list), ~ sapply(., toString))) 
write.csv(file, file="results/db_join.csv", row.names = FALSE)
rm(file)


# Remove empty values
db_join_na <-db_join %>% filter(is.na(bss_total)) 
db_join_s_na <- db_join %>% filter(!is.na(bss_total )) %>% collect() %>% mutate(across(where(is.list), ~ sapply(., toString)))
db_join_s_na <- db_join_s_na %>% filter(!is.na(X4A))%>% filter(trimws(X4A) != "") %>% collect()
write.csv(db_join_s_na, file="results/db_join_s_na.csv", row.names = FALSE)

# Convert fields that should be numeric to numeric
db_join_s_na$PD.5a3a20a.20S <- as.numeric(db_join_s_na$PD.5a3a20a.20S)
db_join_s_na$bss_total     <- as.numeric(db_join_s_na$bss_total)
db_join_s_na$X5PT.S_2  <- as.numeric(db_join_s_na$X5PT.S_2)
db_join_s_na$X5PT.S_1 <- as.numeric(db_join_s_na$X5PT.S_1)
db_join_s_na$X5PD.3b20a.3S <- as.numeric(db_join_s_na$X5PD.3b20a.3S)
db_join_s_na$epiAN.S <- as.numeric(db_join_s_na$epiAN.S)
db_join_s_na$Cortisol..F. <- as.numeric(db_join_s_na$Cortisol..F.)
db_join_s_na$X20a.DHF <- as.numeric(db_join_s_na$X20a.DHF)
db_join_s_na$E1.S <- as.numeric(db_join_s_na$E1.S)
db_join_s_na$X5aTHE <- as.numeric(db_join_s_na$X5aTHE)
db_join_s_na$X20b.DHF <- as.numeric(db_join_s_na$X20b.DHF)
db_join_s_na$X11.CO.4A <- as.numeric(db_join_s_na$X11.CO.4A)
db_join_s_na$epds_total <- as.numeric(db_join_s_na$epds_total)
db_join_s_na$pbq_total <- as.numeric(db_join_s_na$pbq_total)
db_join_s_na$braz_total <- as.numeric(db_join_s_na$braz_total)
db_join_s_na$hamd_total <- as.numeric(db_join_s_na$hamd_total)
write.csv(db_join_s_na, file="results/db_join_s_na.csv", row.names = FALSE)

# Pause
Sys.sleep(1)



# Shapiro-Wilk Test for Hormones ----
Value_H <- db_join_s_na[, 943:983]
Value_H <- data.frame(lapply(Value_H, as.numeric))


# Shapiro-Wilk loop
shapiro_results <- lapply(Value_H, function(x) {
  x <- x[!is.na(x)]
  if (length(x) >= 3 && length(x) <= 5000) {
    shapiro.test(x)
  } else {
    NA
  }})

shapiro_table <- data.frame(
  hormone = names(Value_H),
  W = sapply(shapiro_results, function(x) if (is.list(x)) x$statistic else NA),
  p_value = sapply(shapiro_results, function(x) if (is.list(x)) x$p.value else NA)
)
# View(shapiro_table) # No normal


# Shapiro-wilk BSS-R ----
BSS_value <- db_join_s_na$bss_total
shapiro.test(BSS_value) # Normal

hist(BSS_value, probability = TRUE,
     main = "Distribution of total BSS-R score",
     xlab = "BSS-R",
     col = "lightblue", border = "white")
lines(density(BSS_value), col = "red", lwd = 2) # Add density curve


# Shapiro-wilk EPDS ----
EPDS_value <- db_join_s_na$epds_total
shapiro.test(EPDS_value) # No normal

hist(EPDS_value, probability = TRUE,
     main = "Distribution of total EPDS score",
     xlab = "EPDS",
     col = "lightblue", border = "white")
lines(density(EPDS_value), col = "red", lwd = 2)


# Shapiro-wilk HDRS ----
Hamd_value <- db_join_s_na$hamd_total
shapiro.test(Hamd_value) # No normal

hist(Hamd_value, probability = TRUE,
     main = "Distribution of total HDRS score",
     xlab = "HDRS",
     col = "lightblue", border = "white")
lines(density(Hamd_value), col = "red", lwd = 2)


# Shapiro-wilk PBQ-16 ----
PBQ_value <- db_join_s_na$pbq_total
shapiro.test(PBQ_value) # No normal

hist(PBQ_value, probability = TRUE,
     main = "Distribution of total PBQ-16 score",
     xlab = "PBQ-16",
     col = "lightblue", border = "white")
lines(density(PBQ_value), col = "red", lwd = 2)


# Shapiro-wilk Brazelton Scale ----
Braz_value <- db_join_s_na$braz_total
Braz_value <- Braz_value[!is.na(Braz_value)]
shapiro.test(Braz_value) # No normal

hist(Braz_value, probability = TRUE,
     main = "Distribution of total Brazelton score",
     xlab = "Brazelton Scale",
     col = "lightblue", border = "white")
lines(density(Braz_value), col = "red", lwd = 2)


# Shapiro-wilk maternal age ----
Age_value <- as.numeric(db_join_s_na$edat_1avis)
Age_value <- Age_value[!is.na(Age_value)]
shapiro.test(Age_value) # Normal

hist(Age_value, probability = TRUE,
     main = "Distribution of maternal age",
     xlab = "Maternal age",
     col = "lightblue", border = "white")
lines(density(Age_value), col = "red", lwd = 2)



# Shapiro-wilk parity  ----
Parity_value <- as.numeric(db_join_s_na$tpal_4_v2)
Parity_value <- Parity_value[!is.na(Parity_value)]
shapiro.test(Parity_value) # No normal

hist(Parity_value, probability = TRUE,
     main = "Distribution of parity",
     xlab = "Parity",
     col = "lightblue", border = "white")
lines(density(Parity_value), col = "red", lwd = 2)

# Pause
Sys.sleep(1)



# Create hormone variables ----
hormones <- colnames(db_join_s_na)[943:983]

CORT <- c(
  "Cortisol..F.", "X20a.DHF", "X20b.DHF",
  "X5a.THF", "X5b.THF", "OH.F", "b.cortol",
  "Cortisone", "X20a.DHE", "X20b.DHE", "X5aTHE",
  "X5bTHE", "aCortolone", "X11.Dehydrocorticosterone..A."
)

ANDR <- c(
  "X4A", "X11.CO.4A", "X11b.OH.4A", "X11.CO.T", "T",
  "X5AED.3b17b.SS", "AED.SS", "AD.5a3a17b.3S",
  "X16bOH.DHEA.3S", "AN.G.Etio.G", "epiAN.S",
  "AN.S", "DHEA.S", "X16aOH.DHEA.3S"
)

ESTR <- c(
  "E1.S", "E2.3S"
)

PROG <- c(
  "X21OH.4P..DOC.", "X20a.DH.17OH.4P", "X16OH.4P",
  "X5PD.3b20a.SS", "X5P.S", "alloP.S", "X5PD.3b20a.3S",
  "PD.5a3b20a.20S", "PD.5a3a20a.20S", "X5PT.S_1", "X5PT.S_2"
)

# Pause
Sys.sleep(1)



# CORT ----
## Spearman with BSS-R ----
spearman_bss_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na$bss_total),
    as.numeric(db_join_s_na[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_bss_CORT) <- CORT

spearman_bss_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_bss_CORT, function(x) x$estimate),
  P_value = sapply(spearman_bss_CORT, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_bss_CORT_table)
write.csv(spearman_bss_CORT_table, file="results/spearman_bss_CORT_table.csv")

# Create plot
spearman_bss_CORT_table_t <- spearman_bss_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_bss_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with BSS-R",
    title = "Spearman correlation - BSS-R with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_bss_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
db_join_s_na[CORT] <- lapply(db_join_s_na[CORT], as.numeric)

long_bss_CORT <- db_join_s_na %>%
  dplyr::select(bss_total, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(bss_total) & !is.na(corticosteroid_value))

pvals <- long_bss_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      bss_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_bss_CORT <- long_bss_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_bss_CORT <- ggplot(long_bss_CORT, aes(x = corticosteroid_value, y = bss_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", " 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "BSS-R",
    title = "Spearman correlation - BSS-R with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_bss_CORT.png", scatter_bss_CORT, width = 12, height = 8, dpi = 300)



## Spearman with EPDS ----
db_join_na_epds <-db_join_s_na %>% filter(is.na(epds_total))
db_join_s_na_epds <- db_join_s_na %>%
  filter(!is.na(epds_total )) %>%
  collect() %>% mutate(across(where(is.list), ~ sapply(., toString)))

spearman_EPds_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na_epds$epds_total),
    as.numeric(db_join_s_na_epds[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_EPds_CORT) <- CORT

spearman_epds_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_EPds_CORT, function(x) x$estimate),
  P_value = sapply(spearman_EPds_CORT, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>%
  collect()

# View(spearman_epds_CORT_table)
write.csv(spearman_epds_CORT_table, file="results/spearman_epds_CORT_table.csv")

# Create plot
spearman_epds_CORT_table_t <- spearman_epds_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_epds_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with EPDS",
    title = "Spearman correlation - EPDS with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_epds_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_epds_CORT <- db_join_s_na %>%
  dplyr::select(epds_total, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(epds_total) & !is.na(corticosteroid_value))

pvals <- long_epds_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      epds_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_epds_CORT <- long_epds_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_epds_CORT <- ggplot(long_epds_CORT, aes(x = corticosteroid_value, y = epds_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "EPDS",
    title = "Spearman correlation - EPDS with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_epds_CORT.png", scatter_epds_CORT, width = 12, height = 8, dpi = 300)



## Spearman with HDRS ----
spearman_hamd_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na$hamd_total),
    as.numeric(db_join_s_na[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_hamd_CORT) <- CORT

spearman_hamd_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_hamd_CORT, function(x) x$estimate),
  P_value = sapply(spearman_hamd_CORT, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_hamd_CORT_table)
write.csv(spearman_hamd_CORT_table, file="results/spearman_hamd_CORT_table.csv")

# Create plot
spearman_hamd_CORT_table_t <- spearman_hamd_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_hamd_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with HDRS",
    title = "Spearman correlation - HDRS with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_hamd_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_hamd_CORT <- db_join_s_na %>%
  dplyr::select(hamd_total, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(hamd_total) & !is.na(corticosteroid_value))

pvals <- long_hamd_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      hamd_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_hamd_CORT <- long_hamd_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_hamd_CORT <- ggplot(long_hamd_CORT, aes(x = corticosteroid_value, y = hamd_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "HDRS",
    title = "Spearman correlation - HDRS with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_hamd_CORT.png", scatter_hamd_CORT, width = 12, height = 8, dpi = 300)



## Spearman with PBQ-16 ----
db_join_na_pbq <-db_join_s_na %>% filter(is.na(pbq_total))
db_join_s_na_pbq <- db_join_s_na %>% filter(!is.na(pbq_total )) %>% collect() %>%
  mutate(across(where(is.list), ~ sapply(., toString)))


spearman_pbq_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na_pbq$pbq_total),
    as.numeric(db_join_s_na_pbq[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_pbq_CORT) <- CORT

spearman_pbq_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_pbq_CORT, function(x) x$estimate),
  P_value = sapply(spearman_pbq_CORT, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_pbq_CORT_table)
write.csv(spearman_pbq_CORT_table, file="results/spearman_pbq_CORT_table.csv")

# Create plot
spearman_pbq_CORT_table_t <- spearman_pbq_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_pbq_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with PBQ-16",
    title = "Spearman correlation - PBQ-16 with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_pbq_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_pbq_CORT <- db_join_s_na %>%
  dplyr::select(pbq_total, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(pbq_total) & !is.na(corticosteroid_value))

pvals <- long_pbq_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      pbq_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_pbq_CORT <- long_pbq_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_pbq_CORT <- ggplot(long_pbq_CORT, aes(x = corticosteroid_value, y = pbq_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "PBQ-16",
    title = "Spearman correlation - PBQ-16 with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_pbq_CORT.png", scatter_pbq_CORT, width = 12, height = 8, dpi = 300)



## Spearman with Brazelton Scale ----
db_join_na_braz <-db_join_s_na %>% filter(is.na(braz_total))
db_join_s_na_braz <- db_join_s_na %>% filter(!is.na(braz_total )) %>% collect() %>%
  mutate(across(where(is.list), ~ sapply(., toString)))

spearman_braz_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na_braz$braz_total),
    as.numeric(db_join_s_na_braz[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_braz_CORT) <- CORT

spearman_braz_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_braz_CORT, function(x) x$estimate),
  P_value = sapply(spearman_braz_CORT, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_braz_CORT_table)
write.csv(spearman_braz_CORT_table, file="results/spearman_braz_CORT_table.csv")

# Create plot
spearman_braz_CORT_table_t <- spearman_braz_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_braz_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_braz_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_braz_CORT <- db_join_s_na %>%
  dplyr::select(braz_total, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(braz_total) & !is.na(corticosteroid_value))

pvals <- long_braz_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      braz_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_braz_CORT <- long_braz_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_braz_CORT <- ggplot(long_braz_CORT, aes(x = corticosteroid_value, y = braz_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_braz_CORT.png", scatter_braz_CORT, width = 12, height = 8, dpi = 300)



## Spearman with maternal age ----
db_join_s_na$edat_1avis <- as.numeric(db_join_s_na$edat_1avis)

spearman_age_CORT <- lapply(CORT, function(c) {
  cor.test(
    as.numeric(db_join_s_na$edat_1avis),
    as.numeric(db_join_s_na[[c]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_age_CORT) <- CORT

spearman_age_CORT_table <- data.frame(
  Corticosteroids = CORT,
  Correlation = sapply(spearman_age_CORT, function(x) x$estimate),
  P_value = sapply(spearman_age_CORT, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_age_CORT_table)
write.csv(spearman_age_CORT_table, file="results/spearman_age_CORT_table.csv")

# Create plot
spearman_age_CORT_table_t <- spearman_age_CORT_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_age_CORT_table_t, aes(x = reorder(Corticosteroids, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Corticosteroids",
    y = "Correlation with maternal age",
    title = "Spearman correlation - Maternal age with each corticosteroid",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_age_CORT_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot maternal
long_age_CORT <- db_join_s_na %>%
  dplyr::select(edat_1avis, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroids",
    values_to = "corticosteroid_value"
  ) %>%
  filter(!is.na(edat_1avis) & !is.na(corticosteroid_value))

pvals <- long_age_CORT %>%
  group_by(corticosteroids) %>%
  summarise(
    p_value = cor.test(
      corticosteroid_value,
      edat_1avis,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_age_CORT <- long_age_CORT %>%
  left_join(pvals, by = "corticosteroids")


scatter_age_CORT <- ggplot(long_age_CORT, aes(x = corticosteroid_value, y = edat_1avis)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ corticosteroids, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Corticosteroid value",
    y = "Maternal age",
    title = "Spearman correlation - Maternal age with each corticosteroid"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_age_CORT.png", scatter_age_CORT, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - delivery ----
db_join_s_na$tipus_part_via <- dplyr::recode(
  db_join_s_na$tipus_part_via,
  "1" = "Vaginal",
  "2" = "Vaginal",
  "3" = "C-section",
  "4" = "C-section",
  "5" = "C-section"
)

results_mw_delivery_CORT <- lapply(CORT, function(c) {
  
  vaginal_group <- db_join_s_na %>%
    filter(tipus_part_via == "Vaginal") %>%
    pull(c) %>% na.omit()
  
  c_section_group <- db_join_s_na %>%
    filter(tipus_part_via == "C-section") %>%
    pull(c) %>% na.omit()
  
  test <- wilcox.test(vaginal_group, c_section_group, exact = FALSE)
  
  data.frame(
    Corticosteroids = c,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_delivery_CORT <- do.call(rbind, results_mw_delivery_CORT) %>% arrange(P_value)

# View(mw_table_delivery_CORT)
write.csv(mw_table_delivery_CORT, file="results/mw_table_delivery_CORT.csv")

# Create box plot
long_delivery <- db_join_s_na %>%
  dplyr::select(tipus_part_via, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroid",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_delivery %>%
  group_by(corticosteroid) %>%
  summarise(
    p_value = wilcox.test(value ~ tipus_part_via)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_delivery_t <- long_delivery %>%
  left_join(pvals, by = "corticosteroid")

boxplot_delivery <- ggplot(long_delivery_t, aes(x = tipus_part_via, y = value, fill = tipus_part_via)) +
  
  geom_boxplot(
  aes(color = Significance, size = Significance),
  alpha = 0.7
) +

scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey70")) +
scale_size_manual(values = c("< 0.05" = 0.7, "> 0.05" = 0.4)) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ corticosteroid, scales = "free_y") +
  
  scale_fill_manual(
    name = "Delivery type", 
    values = c("Vaginal" = "steelblue", "C-section" = "tomato")
  ) +
  
  labs(
    x = "Delivery type",
    y = "Corticosteroid value",
    title = "Corticosteroid comparison by delivery type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    axis.text.x = element_blank()
  )

ggsave("results/box_delivery_CORT.png", boxplot_delivery, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - breastfeeding ----
db_join_s_na$lact_0 <- dplyr::recode(
  db_join_s_na$lact_0,
  "0" = "Breastfeeding",
  "1" = "Formula",
  "2" = "Formula"
)

results_mw_breastfeeding_CORT <- lapply(CORT, function(c) {
  
  # Extract values for each group
  group_breastfeeding  <- db_join_s_na %>%
    filter(lact_0 == "Breastfeeding") %>%
    pull(c) %>% na.omit()
  
  group_formula <- db_join_s_na %>%
    filter(lact_0 == "Formula") %>%
    pull(c) %>% na.omit()
  
  test <- wilcox.test(group_breastfeeding , group_formula, exact = FALSE)
  
  data.frame(
    Corticosteroids = c,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_breastfeeding_CORT <- do.call(rbind, results_mw_breastfeeding_CORT) %>% arrange(P_value)

# View(mw_table_breastfeeding_CORT)
write.csv(mw_table_breastfeeding_CORT, file="results/mw_table_breastfeeding_CORT.csv")

# Create box plot
long_breastfeeding <- db_join_s_na %>%
  dplyr::select(lact_0, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroid",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_breastfeeding %>%
  group_by(corticosteroid) %>%
  summarise(
    p_value = wilcox.test(value ~ lact_0, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_breastfeeding_t <- long_breastfeeding %>%
  left_join(pvals, by = "corticosteroid")

boxplot_breastfeeding <- ggplot(
  long_breastfeeding_t,
  aes(x = lact_0, y = value, fill = lact_0)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ corticosteroid, scales = "free_y") +
  
  scale_fill_manual(
    name = "Breastfeeding type", 
    values = c("Breastfeeding" = "steelblue", "Formula" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Breastfeeding type",
    y = "Corticosteroid value",
    title = "Corticosteroid comparison by breastfeeding type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_breastfeeding_CORT.png", boxplot_breastfeeding, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - ATD ----
db_join_s_na$fcodic.factor <- dplyr::recode(
  db_join_s_na$fcodic.factor,
  "no antidepresivo nunca" = "No ATD use",
  "ATD en algun periodo" = "ATD use at some point"
)

results_mw_ATD_CORT <- lapply(CORT, function(c) {

  group_yes <- db_join_s_na %>%
    filter(fcodic == 1) %>%   
    pull(all_of(c)) %>%       
    na.omit()
  
  grupo_no <- db_join_s_na %>%
    filter(fcodic == 0) %>%
    pull(all_of(c)) %>%
    na.omit()
  
  test <- wilcox.test(group_yes, grupo_no, exact = FALSE)
  
  data.frame(
    Corticosteroids = c,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_ATD_CORT <- do.call(rbind, results_mw_ATD_CORT) %>% arrange(P_value)

# View(mw_table_ATD_CORT)
write.csv(mw_table_ATD_CORT, file="results/mw_table_ATD_CORT.csv")

# Create box plot
long_ATD <- db_join_s_na %>%
  dplyr::select(fcodic.factor, dplyr::all_of(CORT)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(CORT),
    names_to = "corticosteroid",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_ATD %>%
  group_by(corticosteroid) %>%
  summarise(
    p_value = wilcox.test(value ~ fcodic.factor, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_ATD_t <- long_ATD %>%
  left_join(pvals, by = "corticosteroid")

boxplot_ATD <- ggplot(
  long_ATD_t,
  aes(x = fcodic.factor, y = value, fill = fcodic.factor)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ corticosteroid, scales = "free_y") +
  
  scale_fill_manual(
    name = "Treatment", 
    values = c("No ATD use" = "steelblue", "ATD use at some point" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Treatment",
    y = "Corticosteroid value",
    title = "Corticosteroid comparison by treatment (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_ATD_CORT.png", boxplot_ATD, width = 12, height = 8, dpi = 300)

# Pause
Sys.sleep(1)





# ANDR ----
## Spearman with BSS-R ----
spearman_bss_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na$bss_total),
    as.numeric(db_join_s_na[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_bss_ANDR) <- ANDR

spearman_bss_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_bss_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_bss_ANDR, function(x) x$p.value)) %>% 
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_bss_ANDR_table)
write.csv(spearman_bss_ANDR_table, file="results/spearman_bss_ANDR_table.csv")

# Create plot
spearman_bss_ANDR_table_t <- spearman_bss_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_bss_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with BSS-R",
    title = "Spearman correlation - BSS-R with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_bss_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
db_join_s_na[ANDR] <- lapply(db_join_s_na[ANDR], as.numeric)

long_bss_ANDR <- db_join_s_na %>%
  dplyr::select(bss_total, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(bss_total) & !is.na(androgen_value))

pvals <- long_bss_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      bss_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_bss_ANDR <- long_bss_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_bss_ANDR <- ggplot(long_bss_ANDR, aes(x = androgen_value, y = bss_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "BSS-R",
    title = "Spearman correlation - BSS-R with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_bss_ANDR.png", scatter_bss_ANDR, width = 12, height = 8, dpi = 300)



## Spearman with EPDS ----
spearman_EPds_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na_epds$epds_total),
    as.numeric(db_join_s_na_epds[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_EPds_ANDR) <- ANDR

spearman_epds_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_EPds_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_EPds_ANDR, function(x) x$p.value)
)%>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_epds_ANDR_table)
write.csv(spearman_epds_ANDR_table, file="results/spearman_epds_ANDR_table.csv")

# Create plot
spearman_epds_ANDR_table_t <- spearman_epds_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_epds_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with EPDS",
    title = "Spearman correlation - EPDS with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_epds_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_epds_ANDR <- db_join_s_na %>%
  dplyr::select(epds_total, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(epds_total) & !is.na(androgen_value))

pvals <- long_epds_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      epds_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_epds_ANDR <- long_epds_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_epds_ANDR <- ggplot(long_epds_ANDR, aes(x = androgen_value, y = epds_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "EPDS",
    title = "Spearman correlation - EPDS with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_epds_ANDR.png", scatter_epds_ANDR, width = 12, height = 8, dpi = 300)



## Spearman with HDRS ----
spearman_hamd_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na$hamd_total),
    as.numeric(db_join_s_na[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_hamd_ANDR) <- ANDR

spearman_hamd_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_hamd_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_hamd_ANDR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_hamd_ANDR_table)
write.csv(spearman_hamd_ANDR_table, file="results/spearman_hamd_ANDR_table.csv")

# Create plot
spearman_hamd_ANDR_table_t <- spearman_hamd_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_hamd_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with HDRS",
    title = "Spearman correlation - HDRS with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_hamd_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_hamd_ANDR <- db_join_s_na %>%
  dplyr::select(hamd_total, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(hamd_total) & !is.na(androgen_value))

pvals <- long_hamd_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      hamd_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_hamd_ANDR <- long_hamd_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_hamd_ANDR <- ggplot(long_hamd_ANDR, aes(x = androgen_value, y = hamd_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "HDRS",
    title = "Spearman correlation - HDRS with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_hamd_ANDR.png", scatter_hamd_ANDR, width = 12, height = 8, dpi = 300)



## Spearman with PBQ-16 ----
spearman_pbq_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na_pbq$pbq_total),
    as.numeric(db_join_s_na_pbq[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_pbq_ANDR) <- ANDR

spearman_pbq_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_pbq_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_pbq_ANDR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_pbq_ANDR_table)
write.csv(spearman_pbq_ANDR_table, file="results/spearman_pbq_ANDR_table.csv")

# Create plot
spearman_pbq_ANDR_table_t <- spearman_pbq_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_pbq_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with PBQ-16",
    title = "Spearman correlation - PBQ-16 with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_pbq_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

# Create scatter plot
long_pbq_ANDR <- db_join_s_na %>%
  dplyr::select(pbq_total, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(pbq_total) & !is.na(androgen_value))

pvals <- long_pbq_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      pbq_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_pbq_ANDR <- long_pbq_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_pbq_ANDR <- ggplot(long_pbq_ANDR, aes(x = androgen_value, y = pbq_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "PBQ-16",
    title = "Spearman correlation - PBQ-16 with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_pbq_ANDR.png", scatter_pbq_ANDR, width = 12, height = 8, dpi = 300)



## Spearman with Brazelton Scale ----
spearman_braz_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na_braz$braz_total),
    as.numeric(db_join_s_na_braz[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_braz_ANDR) <- ANDR

spearman_braz_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_braz_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_braz_ANDR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_braz_ANDR_table)
write.csv(spearman_braz_ANDR_table, file="results/spearman_braz_ANDR_table.csv")

# Create plot
spearman_braz_ANDR_table_t <- spearman_braz_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_braz_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_braz_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_braz_ANDR <- db_join_s_na %>%
  dplyr::select(braz_total, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(braz_total) & !is.na(androgen_value))

pvals <- long_braz_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      braz_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_braz_ANDR <- long_braz_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_braz_ANDR <- ggplot(long_braz_ANDR, aes(x = androgen_value, y = braz_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_braz_ANDR.png", scatter_braz_ANDR, width = 12, height = 8, dpi = 300)



## Spearman with maternal age ----
spearman_age_ANDR <- lapply(ANDR, function(a) {
  cor.test(
    as.numeric(db_join_s_na$edat_1avis),
    as.numeric(db_join_s_na[[a]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_age_ANDR) <- ANDR

spearman_age_ANDR_table <- data.frame(
  Androgens = ANDR,
  Correlation = sapply(spearman_age_ANDR, function(x) x$estimate),
  P_value = sapply(spearman_age_ANDR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_age_ANDR_table)
write.csv(spearman_age_ANDR_table, file="results/spearman_age_ANDR_table.csv")

# Create plot
spearman_age_ANDR_table_t <- spearman_age_ANDR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_age_ANDR_table_t, aes(x = reorder(Androgens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Androgens",
    y = "Correlation with maternal age",
    title = "Spearman correlation - Maternal age with each androgen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8))  

ggsave("results/spearman_age_ANDR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_age_ANDR <- db_join_s_na %>%
  dplyr::select(edat_1avis, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgens",
    values_to = "androgen_value"
  ) %>%
  filter(!is.na(edat_1avis) & !is.na(androgen_value))

pvals <- long_age_ANDR %>%
  group_by(androgens) %>%
  summarise(
    p_value = cor.test(
      androgen_value,
      edat_1avis,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_age_ANDR <- long_age_ANDR %>%
  left_join(pvals, by = "androgens")


scatter_age_ANDR <- ggplot(long_age_ANDR, aes(x = androgen_value, y = edat_1avis)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ androgens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Androgen value",
    y = "Maternal age",
    title = "Spearman correlation - Maternal age with each androgen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_age_ANDR.png", scatter_age_ANDR, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - delivery ----
results_mw_delivery_ANDR <- lapply(ANDR, function(a) {
  
  # Extract values for each group
  vaginal_group <- db_join_s_na %>%
    filter(tipus_part_via == "Vaginal") %>%
    pull(a) %>% na.omit()
  
  c_section_group <- db_join_s_na %>%
    filter(tipus_part_via == "C-section") %>%
    pull(a) %>% na.omit()
  
  test <- wilcox.test(vaginal_group, c_section_group, exact = FALSE)
  
  data.frame(
    Androgens = a,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_delivery_ANDR <- do.call(rbind, results_mw_delivery_ANDR) %>% arrange(P_value)

# View(mw_table_delivery_ANDR)
write.csv(mw_table_delivery_ANDR, file="results/mw_table_delivery_ANDR.csv")

# Create box plot
long_delivery <- db_join_s_na %>%
  dplyr::select(tipus_part_via, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_delivery %>%
  group_by(androgen) %>%
  summarise(
    p_value = wilcox.test(value ~ tipus_part_via)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_delivery_t <- long_delivery %>%
  left_join(pvals, by = "androgen")

boxplot_delivery <- ggplot(long_delivery_t, aes(x = tipus_part_via, y = value, fill = tipus_part_via)) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey70")) +
  scale_size_manual(values = c("< 0.05" = 0.7, "> 0.05" = 0.4)) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ androgen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Delivery type", 
    values = c("Vaginal" = "steelblue", "C-section" = "tomato")
  ) +
  
  labs(
    x = "Delivery type",
    y = "Androgen value",
    title = "Androgen comparison by delivery type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    axis.text.x = element_blank()
  )

ggsave("results/box_delivery_ANDR.png", boxplot_delivery, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - breastfeeding ----
results_mw_breastfeeding_ANDR <- lapply(ANDR, function(a) {
  
  # Extract values for each group
  group_breastfeeding  <- db_join_s_na %>%
    filter(lact_0 == "Breastfeeding") %>%
    pull(a) %>% na.omit()
  
  group_formula <- db_join_s_na %>%
    filter(lact_0 == "Formula") %>%
    pull(a) %>% na.omit()
  
  test <- wilcox.test(group_breastfeeding , group_formula, exact = FALSE)
  
  data.frame(
    Androgens = a,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_breastfeeding_ANDR <- do.call(rbind, results_mw_breastfeeding_ANDR) %>% arrange(P_value)

# View(mw_table_breastfeeding_ANDR)
write.csv(mw_table_breastfeeding_ANDR, file="results/mw_table_breastfeeding_ANDR.csv")

# Create box plot
long_breastfeeding <- db_join_s_na %>%
  dplyr::select(lact_0, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_breastfeeding %>%
  group_by(androgen) %>%
  summarise(
    p_value = wilcox.test(value ~ lact_0, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_breastfeeding_t <- long_breastfeeding %>%
  left_join(pvals, by = "androgen")

boxplot_breastfeeding <- ggplot(
  long_breastfeeding_t,
  aes(x = lact_0, y = value, fill = lact_0)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",  
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ androgen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Breastfeeding type",
    values = c("Breastfeeding" = "steelblue", "Formula" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Breastfeeding type",
    y = "Androgen value",
    title = "Androgen comparison by breastfeeding type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_breastfeeding_ANDR.png", boxplot_breastfeeding, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - ATD ----
results_mw_ATD_ANDR <- lapply(ANDR, function(a) {
  
  group_yes <- db_join_s_na %>%
    filter(fcodic == 1) %>%   
    pull(all_of(a)) %>%       
    na.omit()
  
  grupo_no <- db_join_s_na %>%
    filter(fcodic == 0) %>%
    pull(all_of(a)) %>%
    na.omit()
  
  test <- wilcox.test(group_yes, grupo_no, exact = FALSE)
  
  data.frame(
    Androgens = a,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_ATD_ANDR <- do.call(rbind, results_mw_ATD_ANDR) %>% arrange(P_value)

# View(mw_table_ATD_ANDR)
write.csv(mw_table_ATD_ANDR, file="results/mw_table_ATD_ANDR.csv")

# Create box plot
long_ATD <- db_join_s_na %>%
  dplyr::select(fcodic.factor, dplyr::all_of(ANDR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ANDR),
    names_to = "androgen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_ATD %>%
  group_by(androgen) %>%
  summarise(
    p_value = wilcox.test(value ~ fcodic.factor, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_ATD_t <- long_ATD %>%
  left_join(pvals, by = "androgen")

boxplot_ATD <- ggplot(
  long_ATD_t,
  aes(x = fcodic.factor, y = value, fill = fcodic.factor)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ androgen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Treatment", 
    values = c("No ATD use" = "steelblue", "ATD use at some point" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Treatment",
    y = "Androgen value",
    title = "Androgen comparison by treatment (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_ATD_ANDR.png", boxplot_ATD, width = 12, height = 8, dpi = 300)

# Pause
Sys.sleep(1)








# ESTR ----
## Spearman with BSS-R ----
spearman_bss_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na$bss_total),
    as.numeric(db_join_s_na[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_bss_ESTR) <- ESTR

spearman_bss_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_bss_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_bss_ESTR, function(x) x$p.value)) %>%
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_bss_ESTR_table)
write.csv(spearman_bss_ESTR_table, file="results/spearman_bss_ESTR_table.csv")

# Create plot
spearman_bss_ESTR_table_t <- spearman_bss_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_bss_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with BSS-R",
    title = "Spearman correlation - BSS-R with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_bss_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
db_join_s_na[ESTR] <- lapply(db_join_s_na[ESTR], as.numeric)

long_bss_ESTR <- db_join_s_na %>%
  dplyr::select(bss_total, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(bss_total) & !is.na(estrogen_value))

pvals <- long_bss_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      bss_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_bss_ESTR <- long_bss_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_bss_ESTR <- ggplot(long_bss_ESTR, aes(x = estrogen_value, y = bss_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "BSS-R",
    title = "Spearman correlation - BSS-R with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_bss_ESTR.png", scatter_bss_ESTR, width = 12, height = 8, dpi = 300)



## Spearman with EPDS ----
spearman_EPds_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na_epds$epds_total),
    as.numeric(db_join_s_na_epds[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_EPds_ESTR) <- ESTR

spearman_epds_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_EPds_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_EPds_ESTR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_epds_ESTR_table)
write.csv(spearman_epds_ESTR_table, file="results/spearman_epds_ESTR_table.csv")

# Create plot
spearman_epds_ESTR_table_t <- spearman_epds_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_epds_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with EPDS",
    title = "Spearman correlation - EPDS with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_epds_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_epds_ESTR <- db_join_s_na %>%
  dplyr::select(epds_total, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(epds_total) & !is.na(estrogen_value))

pvals <- long_epds_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      epds_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_epds_ESTR <- long_epds_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_epds_ESTR <- ggplot(long_epds_ESTR, aes(x = estrogen_value, y = epds_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "EPDS",
    title = "Spearman correlation - EPDS with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_epds_ESTR.png", scatter_epds_ESTR, width = 12, height = 8, dpi = 300)



## Spearman with HDRS ----
spearman_hamd_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na$hamd_total),
    as.numeric(db_join_s_na[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_hamd_ESTR) <- ESTR

spearman_hamd_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_hamd_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_hamd_ESTR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_hamd_ESTR_table)
write.csv(spearman_hamd_ESTR_table, file="results/spearman_hamd_ESTR_table.csv")

# Create plot
spearman_hamd_ESTR_table_t <- spearman_hamd_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_hamd_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with HDRS",
    title = "Spearman correlation - HDRS with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_hamd_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_hamd_ESTR <- db_join_s_na %>%
  dplyr::select(hamd_total, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(hamd_total) & !is.na(estrogen_value))

pvals <- long_hamd_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      hamd_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_hamd_ESTR <- long_hamd_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_hamd_ESTR <- ggplot(long_hamd_ESTR, aes(x = estrogen_value, y = hamd_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "HDRS",
    title = "Spearman correlation - HDRS with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_hamd_ESTR.png", scatter_hamd_ESTR, width = 12, height = 8, dpi = 300)



## Spearman with PBQ-16 ----
spearman_pbq_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na_pbq$pbq_total),
    as.numeric(db_join_s_na_pbq[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_pbq_ESTR) <- ESTR

spearman_pbq_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_pbq_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_pbq_ESTR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_pbq_ESTR_table)
write.csv(spearman_pbq_ESTR_table, file="results/spearman_pbq_ESTR_table.csv")

# Create plot
spearman_pbq_ESTR_table_t <- spearman_pbq_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_pbq_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with PBQ-16",
    title = "Spearman correlation - PBQ-16 with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_pbq_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_pbq_ESTR <- db_join_s_na %>%
  dplyr::select(pbq_total, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(pbq_total) & !is.na(estrogen_value))

pvals <- long_pbq_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      pbq_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_pbq_ESTR <- long_pbq_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_pbq_ESTR <- ggplot(long_pbq_ESTR, aes(x = estrogen_value, y = pbq_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "PBQ-16",
    title = "Spearman correlation - PBQ-16 with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_pbq_ESTR.png", scatter_pbq_ESTR, width = 12, height = 8, dpi = 300)



## Spearman with Brazelton Scale ----
spearman_braz_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na_braz$braz_total),
    as.numeric(db_join_s_na_braz[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_braz_ESTR) <- ESTR

spearman_braz_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_braz_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_braz_ESTR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_braz_ESTR_table)
write.csv(spearman_braz_ESTR_table, file="results/spearman_braz_ESTR_table.csv")

# Create plot
spearman_braz_ESTR_table_t <- spearman_braz_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_braz_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_braz_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_braz_ESTR <- db_join_s_na %>%
  dplyr::select(braz_total, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(braz_total) & !is.na(estrogen_value))

pvals <- long_braz_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      braz_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_braz_ESTR <- long_braz_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_braz_ESTR <- ggplot(long_braz_ESTR, aes(x = estrogen_value, y = braz_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_braz_ESTR.png", scatter_braz_ESTR, width = 12, height = 8, dpi = 300)



## Spearman with maternal age ----
spearman_age_ESTR <- lapply(ESTR, function(e) {
  cor.test(
    as.numeric(db_join_s_na$edat_1avis),
    as.numeric(db_join_s_na[[e]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_age_ESTR) <- ESTR

spearman_age_ESTR_table <- data.frame(
  Estrogens = ESTR,
  Correlation = sapply(spearman_age_ESTR, function(x) x$estimate),
  P_value = sapply(spearman_age_ESTR, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_age_ESTR_table)
write.csv(spearman_age_ESTR_table, file="results/spearman_age_ESTR_table.csv")

# Create plot
spearman_age_ESTR_table_t <- spearman_age_ESTR_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_age_ESTR_table_t, aes(x = reorder(Estrogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30")
  ) +
  labs(
    x = "Estrogens",
    y = "Correlation with maternal age",
    title = "Spearman correlation - Maternal age with each estrogen",
    fill = "P-value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/spearman_age_ESTR_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_age_ESTR <- db_join_s_na %>%
  dplyr::select(edat_1avis, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogens",
    values_to = "estrogen_value"
  ) %>%
  filter(!is.na(edat_1avis) & !is.na(estrogen_value))

pvals <- long_age_ESTR %>%
  group_by(estrogens) %>%
  summarise(
    p_value = cor.test(
      estrogen_value,
      edat_1avis,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_age_ESTR <- long_age_ESTR %>%
  left_join(pvals, by = "estrogens")


scatter_age_ESTR <- ggplot(long_age_ESTR, aes(x = estrogen_value, y = edat_1avis)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ estrogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Estrogen value",
    y = "Maternal age",
    title = "Spearman correlation - Maternal age with each estrogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_age_ESTR.png", scatter_age_ESTR, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - delivery ----
results_mw_delivery_ESTR <- lapply(ESTR, function(e) {
  
  # Extract values for each group
  vaginal_group <- db_join_s_na %>%
    filter(tipus_part_via == "Vaginal") %>%
    pull(e) %>% na.omit()
  
  c_section_group <- db_join_s_na %>%
    filter(tipus_part_via == "C-section") %>%
    pull(e) %>% na.omit()
  
  test <- wilcox.test(vaginal_group, c_section_group, exact = FALSE)
  
  data.frame(
    Estrogens = e,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_delivery_ESTR <- do.call(rbind, results_mw_delivery_ESTR) %>% arrange(P_value)

# View(mw_table_delivery_ESTR)
write.csv(mw_table_delivery_ESTR, file="results/mw_table_delivery_ESTR.csv")

# Create box plot
long_delivery <- db_join_s_na %>%
  dplyr::select(tipus_part_via, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_delivery %>%
  group_by(estrogen) %>%
  summarise(
    p_value = wilcox.test(value ~ tipus_part_via)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_delivery_t <- long_delivery %>%
  left_join(pvals, by = "estrogen")

boxplot_delivery <- ggplot(long_delivery_t, aes(x = tipus_part_via, y = value, fill = tipus_part_via)) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey70")) +
  scale_size_manual(values = c("< 0.05" = 0.7, "> 0.05" = 0.4)) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ estrogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Delivery type", 
    values = c("Vaginal" = "steelblue", "C-section" = "tomato")
  ) +
  
  labs(
    x = "Delivery type",
    y = "Estrogen value",
    title = "Estrogen comparison by delivery type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    axis.text.x = element_blank()
  )

ggsave("results/box_delivery_ESTR.png", boxplot_delivery, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - breastfeeding ----
results_mw_breastfeeding_ESTR <- lapply(ESTR, function(e) {
  
  # Extract values for each group
  group_breastfeeding  <- db_join_s_na %>%
    filter(lact_0 == "Breastfeeding") %>%
    pull(e) %>% na.omit()
  
  group_formula <- db_join_s_na %>%
    filter(lact_0 == "Formula") %>%
    pull(e) %>% na.omit()
  
  test <- wilcox.test(group_breastfeeding , group_formula, exact = FALSE)
  
  data.frame(
    Estrogens = e,
    U = test$statistic,
    P_value = test$p.value
  )
})


mw_table_breastfeeding_ESTR <- do.call(rbind, results_mw_breastfeeding_ESTR) %>% arrange(P_value)

# View(mw_table_breastfeeding_ESTR)
write.csv(mw_table_breastfeeding_ESTR, file="results/mw_table_breastfeeding_ESTR.csv")

# Create box plot
long_breastfeeding <- db_join_s_na %>%
  dplyr::select(lact_0, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_breastfeeding %>%
  group_by(estrogen) %>%
  summarise(
    p_value = wilcox.test(value ~ lact_0, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_breastfeeding_t <- long_breastfeeding %>%
  left_join(pvals, by = "estrogen")

boxplot_breastfeeding <- ggplot(
  long_breastfeeding_t,
  aes(x = lact_0, y = value, fill = lact_0)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",  
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ estrogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Breastfeeding type",
    values = c("Breastfeeding" = "steelblue", "Formula" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Breastfeeding type",
    y = "Estrogen value",
    title = "Estrogen comparison by breastfeeding type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_breastfeeding_ESTR.png", boxplot_breastfeeding, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - ATD ----
results_mw_ATD_ESTR <- lapply(ESTR, function(e) {
  
  group_yes <- db_join_s_na %>%
    filter(fcodic == 1) %>%   
    pull(all_of(e)) %>%       
    na.omit()
  
  grupo_no <- db_join_s_na %>%
    filter(fcodic == 0) %>%
    pull(all_of(e)) %>%
    na.omit()
  
  test <- wilcox.test(group_yes, grupo_no, exact = FALSE)
  
  data.frame(
    Estrogens = e,
    U = test$statistic,
    P_value = test$p.value
  )
})


mw_table_ATD_ESTR <- do.call(rbind, results_mw_ATD_ESTR) %>% arrange(P_value)

# View(mw_table_ATD_ESTR)
write.csv(mw_table_ATD_ESTR, file="results/mw_table_ATD_ESTR.csv")

# Create box plot
long_ATD <- db_join_s_na %>%
  dplyr::select(fcodic.factor, dplyr::all_of(ESTR)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(ESTR),
    names_to = "estrogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_ATD %>%
  group_by(estrogen) %>%
  summarise(
    p_value = wilcox.test(value ~ fcodic.factor, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_ATD_t <- long_ATD %>%
  left_join(pvals, by = "estrogen")

boxplot_ATD <- ggplot(
  long_ATD_t,
  aes(x = fcodic.factor, y = value, fill = fcodic.factor)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ estrogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Treatment", 
    values = c("No ATD use" = "steelblue", "ATD use at some point" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Treatment",
    y = "Estrogen value",
    title = "Estrogen comparison by treatment (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_ATD_ESTR.png", boxplot_ATD, width = 12, height = 8, dpi = 300)

# Pause
Sys.sleep(1)







# PROG ----
## Spearman with BSS-R ----
spearman_bss_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na$bss_total),
    as.numeric(db_join_s_na[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_bss_PROG) <- PROG

spearman_bss_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_bss_PROG, function(x) x$estimate),
  P_value = sapply(spearman_bss_PROG, function(x) x$p.value) ) %>% 
  arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_bss_PROG_table)
write.csv(spearman_bss_PROG_table, file="results/spearman_bss_PROG_table.csv")

# Create plot
spearman_bss_PROG_table_t <- spearman_bss_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_bss_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with BSS-R",
    title = "Spearman correlation - BSS-R with each progestogen",
    fill = "P-value"
    ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_bss_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
db_join_s_na[PROG] <- lapply(db_join_s_na[PROG], as.numeric)

long_bss_PROG <- db_join_s_na %>%
  dplyr::select(bss_total, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(bss_total) & !is.na(progestogen_value))

pvals <- long_bss_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      bss_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_bss_PROG <- long_bss_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_bss_PROG <- ggplot(long_bss_PROG, aes(x = progestogen_value, y = bss_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "BSS-R",
    title = "Spearman correlation - BSS-R with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_bss_PROG.png", scatter_bss_PROG, width = 12, height = 8, dpi = 300)



## Spearman with EPDS ----
spearman_EPds_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na_epds$epds_total),
    as.numeric(db_join_s_na_epds[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_EPds_PROG) <- PROG

spearman_epds_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_EPds_PROG, function(x) x$estimate),
  P_value = sapply(spearman_EPds_PROG, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_epds_PROG_table)
write.csv(spearman_epds_PROG_table, file="results/spearman_epds_PROG_table.csv")

# Create plot
spearman_epds_PROG_table_t <- spearman_epds_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_epds_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with EPDS",
    title = "Spearman correlation - EPDS with each progestogen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_epds_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_epds_PROG <- db_join_s_na %>%
  dplyr::select(epds_total, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(epds_total) & !is.na(progestogen_value))

pvals <- long_epds_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      epds_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_epds_PROG <- long_epds_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_epds_PROG <- ggplot(long_epds_PROG, aes(x = progestogen_value, y = epds_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "EPDS",
    title = "Spearman correlation - EPDS with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_epds_PROG.png", scatter_epds_PROG, width = 12, height = 8, dpi = 300)



## Spearman with HDRS ----
spearman_hamd_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na$hamd_total),
    as.numeric(db_join_s_na[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_hamd_PROG) <- PROG

spearman_hamd_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_hamd_PROG, function(x) x$estimate),
  P_value = sapply(spearman_hamd_PROG, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_hamd_PROG_table)
write.csv(spearman_hamd_PROG_table, file="results/spearman_hamd_PROG_table.csv")

# Create plot
spearman_hamd_PROG_table_t <- spearman_hamd_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_hamd_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with HDRS",
    title = "Spearman correlation - HDRS with each progestogen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_hamd_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_hamd_PROG <- db_join_s_na %>%
  dplyr::select(hamd_total, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(hamd_total) & !is.na(progestogen_value))

pvals <- long_hamd_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      hamd_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_hamd_PROG <- long_hamd_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_hamd_PROG <- ggplot(long_hamd_PROG, aes(x = progestogen_value, y = hamd_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "HDRS",
    title = "Spearman correlation - HDRS with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_hamd_PROG.png", scatter_hamd_PROG, width = 12, height = 8, dpi = 300)



## Spearman with PBQ-16 ----
spearman_pbq_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na_pbq$pbq_total),
    as.numeric(db_join_s_na_pbq[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_pbq_PROG) <- PROG

spearman_pbq_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_pbq_PROG, function(x) x$estimate),
  P_value = sapply(spearman_pbq_PROG, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_pbq_PROG_table)
write.csv(spearman_pbq_PROG_table, file="results/spearman_pbq_PROG_table.csv")

# Create plot
spearman_pbq_PROG_table_t <- spearman_pbq_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_pbq_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with PBQ-16",
    title = "Spearman correlation - PBQ-16 with each progestogen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_pbq_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_pbq_PROG <- db_join_s_na %>%
  dplyr::select(pbq_total, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(pbq_total) & !is.na(progestogen_value))

pvals <- long_pbq_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      pbq_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_pbq_PROG <- long_pbq_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_pbq_PROG <- ggplot(long_pbq_PROG, aes(x = progestogen_value, y = pbq_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "PBQ-16",
    title = "Spearman correlation - PBQ-16 with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_pbq_PROG.png", scatter_pbq_PROG, width = 12, height = 8, dpi = 300)



## Spearman with Brazelton Scale ----
spearman_braz_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na_braz$braz_total),
    as.numeric(db_join_s_na_braz[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_braz_PROG) <- PROG

spearman_braz_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_braz_PROG, function(x) x$estimate),
  P_value = sapply(spearman_braz_PROG, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_braz_PROG_table)
write.csv(spearman_braz_PROG_table, file="results/spearman_braz_PROG_table.csv")

# Create plot
spearman_braz_PROG_table_t <- spearman_braz_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_braz_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each progestogen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_braz_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_braz_PROG <- db_join_s_na %>%
  dplyr::select(braz_total, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(braz_total) & !is.na(progestogen_value))

pvals <- long_braz_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      braz_total,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_braz_PROG <- long_braz_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_braz_PROG <- ggplot(long_braz_PROG, aes(x = progestogen_value, y = braz_total)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "Brazelton Scale",
    title = "Spearman correlation - Brazelton Scale with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_braz_PROG.png", scatter_braz_PROG, width = 12, height = 8, dpi = 300)



## Spearman with maternal age ----
spearman_age_PROG <- lapply(PROG, function(p) {
  cor.test(
    as.numeric(db_join_s_na$edat_1avis),
    as.numeric(db_join_s_na[[p]]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

names(spearman_age_PROG) <- PROG

spearman_age_PROG_table <- data.frame(
  Progestogens = PROG,
  Correlation = sapply(spearman_age_PROG, function(x) x$estimate),
  P_value = sapply(spearman_age_PROG, function(x) x$p.value)
) %>% arrange(desc(abs(Correlation))) %>% collect()

# View(spearman_age_PROG_table)
write.csv(spearman_age_PROG_table, file="results/spearman_age_PROG_table.csv")

# Create plot
spearman_age_PROG_table_t <- spearman_age_PROG_table %>%
  mutate(
    signif = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

b <- ggplot(spearman_age_PROG_table_t, aes(x = reorder(Progestogens, Correlation), y = Correlation, fill = signif)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("< 0.05" = "steelblue", "> 0.05" = "grey30") 
  ) +
  labs(
    x = "Progestogens",
    y = "Correlation with maternal age",
    title = "Spearman correlation - Maternal age with each progestogen",
    fill = "P-value"
  ) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8)) 

ggsave("results/spearman_age_PROG_graphic.png", plot = b, width = 6, height = 4, dpi = 300)

rm(b)

# Create scatter plot
long_age_PROG <- db_join_s_na %>%
  dplyr::select(edat_1avis, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogens",
    values_to = "progestogen_value"
  ) %>%
  filter(!is.na(edat_1avis) & !is.na(progestogen_value))

pvals <- long_age_PROG %>%
  group_by(progestogens) %>%
  summarise(
    p_value = cor.test(
      progestogen_value,
      edat_1avis,
      method = "spearman",
      exact = FALSE
    )$p.value
  ) %>%
  mutate(
    signif = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_age_PROG <- long_age_PROG %>%
  left_join(pvals, by = "progestogens")


scatter_age_PROG <- ggplot(long_age_PROG, aes(x = progestogen_value, y = edat_1avis)) +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(
    aes(color = signif),
    method = "lm",
    se = FALSE
  ) +
  
  stat_cor(
    method = "spearman",
    color = "blue",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    size = 3
  ) +
  
  facet_wrap(~ progestogens, scales = "free_x") +
  
  scale_color_manual(
    values = c("< 0.05" = "red", "> 0.05" = "grey60")
  ) +
  
  labs(
    x = "Progestogen value",
    y = "Maternal age",
    title = "Spearman correlation - Maternal age with each progestogen"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

ggsave("results/scatter_age_PROG.png", scatter_age_PROG, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - delivery ----
results_mw_delivery_PROG <- lapply(PROG, function(p) {
  
  vaginal_group <- db_join_s_na %>%
    filter(tipus_part_via == "Vaginal") %>%
    pull(p) %>% na.omit()
  
  c_section_group <- db_join_s_na %>%
    filter(tipus_part_via == "C-section") %>%
    pull(p) %>% na.omit()
  
  test <- wilcox.test(vaginal_group, c_section_group, exact = FALSE)
  
  data.frame(
    Progestogens = p,
    U = test$statistic,
    P_value = test$p.value
  )
})


mw_table_delivery_PROG <- do.call(rbind, results_mw_delivery_PROG) %>% arrange(P_value)

# View(mw_table_delivery_PROG)
write.csv(mw_table_delivery_PROG, file="results/mw_table_delivery_PROG.csv")

# Create box plot
long_delivery <- db_join_s_na %>%
  dplyr::select(tipus_part_via, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_delivery %>%
  group_by(progestogen) %>%
  summarise(
    p_value = wilcox.test(value ~ tipus_part_via)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_delivery_t <- long_delivery %>%
  left_join(pvals, by = "progestogen")

boxplot_delivery <- ggplot(long_delivery_t, aes(x = tipus_part_via, y = value, fill = tipus_part_via)) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey70")) +
  scale_size_manual(values = c("< 0.05" = 0.7, "> 0.05" = 0.4)) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ progestogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Delivery type", 
    values = c("Vaginal" = "steelblue", "C-section" = "tomato")
  ) +
  
  labs(
    x = "Delivery type",
    y = "Progestogen value",
    title = "Progestogen comparison by delivery type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    axis.text.x = element_blank()
  )

ggsave("results/box_delivery_PROG.png", boxplot_delivery, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test - breastfeeding ----
results_mw_breastfeeding_PROG <- lapply(PROG, function(p) {
  
  # Extract values for each group
  group_breastfeeding  <- db_join_s_na %>%
    filter(lact_0 == "Breastfeeding") %>%
    pull(p) %>% na.omit()
  
  group_formula <- db_join_s_na %>%
    filter(lact_0 == "Formula") %>%
    pull(p) %>% na.omit()
  
  test <- wilcox.test(group_breastfeeding , group_formula, exact = FALSE)
  
  data.frame(
    Progestogens = p,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_breastfeeding_PROG <- do.call(rbind, results_mw_breastfeeding_PROG) %>% arrange(P_value)

# View(mw_table_breastfeeding_PROG)
write.csv(mw_table_breastfeeding_PROG, file="results/mw_table_breastfeeding_PROG.csv")

# Create box plot
long_breastfeeding <- db_join_s_na %>%
  dplyr::select(lact_0, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_breastfeeding %>%
  group_by(progestogen) %>%
  summarise(
    p_value = wilcox.test(value ~ lact_0, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_breastfeeding_t <- long_breastfeeding %>%
  left_join(pvals, by = "progestogen")

boxplot_breastfeeding <- ggplot(
  long_breastfeeding_t,
  aes(x = lact_0, y = value, fill = lact_0)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",  
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ progestogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Breastfeeding type",
    values = c("Breastfeeding" = "steelblue", "Formula" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Breastfeeding type",
    y = "Progestogen value",
    title = "Progestogen comparison by breastfeeding type (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_breastfeeding_PROG.png", boxplot_breastfeeding, width = 12, height = 8, dpi = 300)



## Mann-Whitney U test -  ATD ----
results_mw_ATD_PROG <- lapply(PROG, function(p) {
  
  group_yes <- db_join_s_na %>%
    filter(fcodic == 1) %>%   
    pull(all_of(p)) %>%       
    na.omit()
  
  grupo_no <- db_join_s_na %>%
    filter(fcodic == 0) %>%
    pull(all_of(p)) %>%
    na.omit()
  
  test <- wilcox.test(group_yes, grupo_no, exact = FALSE)
  
  data.frame(
    Progestogens = p,
    U = test$statistic,
    P_value = test$p.value
  )
})

mw_table_ATD_PROG <- do.call(rbind, results_mw_ATD_PROG) %>% arrange(P_value)

# View(mw_table_ATD_PROG)
write.csv(mw_table_ATD_PROG, file="results/mw_table_ATD_PROG.csv")

# Create box plot
long_ATD <- db_join_s_na %>%
  dplyr::select(fcodic.factor, dplyr::all_of(PROG)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(PROG),
    names_to = "progestogen",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

pvals <- long_ATD %>%
  group_by(progestogen) %>%
  summarise(
    p_value = wilcox.test(value ~ fcodic.factor, exact = FALSE)$p.value
  ) %>%
  mutate(
    Significance = ifelse(p_value < 0.05, "< 0.05", "> 0.05")
  )

long_ATD_t <- long_ATD %>%
  left_join(pvals, by = "progestogen")

boxplot_ATD <- ggplot(
  long_ATD_t,
  aes(x = fcodic.factor, y = value, fill = fcodic.factor)
) +
  
  geom_boxplot(
    aes(color = Significance, size = Significance),
    alpha = 0.7
  ) +
  
  stat_compare_means(
    method = "wilcox.test",   
    label = "p.format",  
    color = "blue",
    label.y.npc = 0.90,
    label.x.npc = 0.5
  ) +
  
  facet_wrap(~ progestogen, scales = "free_y") +
  
  scale_fill_manual(
    name = "Treatment", 
    values = c("No ATD use" = "steelblue", "ATD use at some point" = "tomato")
  ) +
  
  scale_color_manual(
    name = "Significance",
    values = c("< 0.05" = "red", "> 0.05" = "grey70")
  ) +
  
  scale_size_manual(
    name = "Significance",
    values = c("< 0.05" = 0.7, "> 0.05" = 0.4)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  
  labs(
    x = "Treatment",
    y = "Progestogen value",
    title = "Progestogen comparison by treatment (Mann-Whitney U test)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(color = "black", face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_blank()
  )

ggsave("results/box_ATD_PROG.png", boxplot_ATD, width = 12, height = 8, dpi = 300)

# Pause
Sys.sleep(1)






# FDR corrections - Spearman / Mann-Whitney U test ----
spearman_bss_ANDR_table_FDR <- spearman_bss_ANDR_table
spearman_bss_ANDR_table_FDR$FDR <- p.adjust(spearman_bss_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_bss_ANDR_table_FDR)
write.csv(spearman_bss_ANDR_table_FDR, file="results/spearman_bss_ANDR_table_FDR.csv")

spearman_bss_CORT_table_FDR <- spearman_bss_CORT_table
spearman_bss_CORT_table_FDR$FDR <- p.adjust(spearman_bss_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_bss_CORT_table_FDR)
write.csv(spearman_bss_CORT_table_FDR, file="results/spearman_bss_CORT_table_FDR.csv")

spearman_bss_ESTR_table_FDR <- spearman_bss_ESTR_table
spearman_bss_ESTR_table_FDR$FDR <- p.adjust(spearman_bss_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_bss_ESTR_table_FDR)
write.csv(spearman_bss_ESTR_table_FDR, file="results/spearman_bss_ESTR_table_FDR.csv")

spearman_bss_PROG_table_FDR <- spearman_bss_PROG_table
spearman_bss_PROG_table_FDR$FDR <- p.adjust(spearman_bss_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_bss_PROG_table_FDR)
write.csv(spearman_bss_PROG_table_FDR, file="results/spearman_bss_PROG_table_FDR.csv")

spearman_epds_ANDR_table_FDR <- spearman_epds_ANDR_table
spearman_epds_ANDR_table_FDR$FDR <- p.adjust(spearman_epds_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_epds_ANDR_table_FDR)
write.csv(spearman_epds_ANDR_table_FDR, file="results/spearman_epds_ANDR_table_FDR.csv")

spearman_epds_CORT_table_FDR <- spearman_epds_CORT_table
spearman_epds_CORT_table_FDR$FDR <- p.adjust(spearman_epds_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_epds_CORT_table_FDR)
write.csv(spearman_epds_CORT_table_FDR, file="results/spearman_epds_CORT_table_FDR.csv")

spearman_epds_ESTR_table_FDR <- spearman_epds_ESTR_table
spearman_epds_ESTR_table_FDR$FDR <- p.adjust(spearman_epds_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_epds_ESTR_table_FDR)
write.csv(spearman_epds_ESTR_table_FDR, file="results/spearman_epds_ESTR_table_FDR.csv")

spearman_epds_PROG_table_FDR <- spearman_epds_PROG_table
spearman_epds_PROG_table_FDR$FDR <- p.adjust(spearman_epds_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_epds_PROG_table_FDR)
write.csv(spearman_epds_PROG_table_FDR, file="results/spearman_epds_PROG_table_FDR.csv")

spearman_hamd_ANDR_table_FDR <- spearman_hamd_ANDR_table
spearman_hamd_ANDR_table_FDR$FDR <- p.adjust(spearman_hamd_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_hamd_ANDR_table_FDR)
write.csv(spearman_hamd_ANDR_table_FDR, file="results/spearman_hamd_ANDR_table_FDR.csv")

spearman_hamd_CORT_table_FDR <- spearman_hamd_CORT_table
spearman_hamd_CORT_table_FDR$FDR <- p.adjust(spearman_hamd_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_hamd_CORT_table_FDR)
write.csv(spearman_hamd_CORT_table_FDR, file="results/spearman_hamd_CORT_table_FDR.csv")

spearman_hamd_ESTR_table_FDR <- spearman_hamd_ESTR_table
spearman_hamd_ESTR_table_FDR$FDR <- p.adjust(spearman_hamd_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_hamd_ESTR_table_FDR)
write.csv(spearman_hamd_ESTR_table_FDR, file="results/spearman_hamd_ESTR_table_FDR.csv")

spearman_hamd_PROG_table_FDR <- spearman_hamd_PROG_table
spearman_hamd_PROG_table_FDR$FDR <- p.adjust(spearman_hamd_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_hamd_PROG_table_FDR)
write.csv(spearman_hamd_PROG_table_FDR, file="results/spearman_hamd_PROG_table_FDR.csv")

spearman_pbq_ANDR_table_FDR <- spearman_pbq_ANDR_table
spearman_pbq_ANDR_table_FDR$FDR <- p.adjust(spearman_pbq_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_pbq_ANDR_table_FDR)
write.csv(spearman_pbq_ANDR_table_FDR, file="results/spearman_pbq_ANDR_table_FDR.csv")

spearman_pbq_CORT_table_FDR <- spearman_pbq_CORT_table
spearman_pbq_CORT_table_FDR$FDR <- p.adjust(spearman_pbq_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_pbq_CORT_table_FDR)
write.csv(spearman_pbq_CORT_table_FDR, file="results/spearman_pbq_CORT_table_FDR.csv")

spearman_pbq_ESTR_table_FDR <- spearman_pbq_ESTR_table
spearman_pbq_ESTR_table_FDR$FDR <- p.adjust(spearman_pbq_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_pbq_ESTR_table_FDR)
write.csv(spearman_pbq_ESTR_table_FDR, file="results/spearman_pbq_ESTR_table_FDR.csv")

spearman_pbq_PROG_table_FDR <- spearman_pbq_PROG_table
spearman_pbq_PROG_table_FDR$FDR <- p.adjust(spearman_pbq_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_pbq_PROG_table_FDR)
write.csv(spearman_pbq_PROG_table_FDR, file="results/spearman_pbq_PROG_table_FDR.csv")

spearman_braz_ANDR_table_FDR <- spearman_braz_ANDR_table
spearman_braz_ANDR_table_FDR$FDR <- p.adjust(spearman_braz_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_braz_ANDR_table_FDR)
write.csv(spearman_braz_ANDR_table_FDR, file="results/spearman_braz_ANDR_table_FDR.csv")

spearman_braz_CORT_table_FDR <- spearman_braz_CORT_table
spearman_braz_CORT_table_FDR$FDR <- p.adjust(spearman_braz_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_braz_CORT_table_FDR)
write.csv(spearman_braz_CORT_table_FDR, file="results/spearman_braz_CORT_table_FDR.csv")

spearman_braz_ESTR_table_FDR <- spearman_braz_ESTR_table
spearman_braz_ESTR_table_FDR$FDR <- p.adjust(spearman_braz_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_braz_ESTR_table_FDR)
write.csv(spearman_braz_ESTR_table_FDR, file="results/spearman_braz_ESTR_table_FDR.csv")

spearman_braz_PROG_table_FDR <- spearman_braz_PROG_table
spearman_braz_PROG_table_FDR$FDR <- p.adjust(spearman_braz_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_braz_PROG_table_FDR)
write.csv(spearman_braz_PROG_table_FDR, file="results/spearman_braz_PROG_table_FDR.csv")

spearman_age_ANDR_table_FDR <- spearman_age_ANDR_table
spearman_age_ANDR_table_FDR$FDR <- p.adjust(spearman_age_ANDR_table_FDR$P_value, method = "fdr")
# View(spearman_age_ANDR_table_FDR)
write.csv(spearman_age_ANDR_table_FDR, file="results/spearman_age_ANDR_table_FDR.csv")

spearman_age_CORT_table_FDR <- spearman_age_CORT_table
spearman_age_CORT_table_FDR$FDR <- p.adjust(spearman_age_CORT_table_FDR$P_value, method = "fdr")
# View(spearman_age_CORT_table_FDR)
write.csv(spearman_age_CORT_table_FDR, file="results/spearman_age_CORT_table_FDR.csv")

spearman_age_ESTR_table_FDR <- spearman_age_ESTR_table
spearman_age_ESTR_table_FDR$FDR <- p.adjust(spearman_age_ESTR_table_FDR$P_value, method = "fdr")
# View(spearman_age_ESTR_table_FDR)
write.csv(spearman_age_ESTR_table_FDR, file="results/spearman_age_ESTR_table_FDR.csv")

spearman_age_PROG_table_FDR <- spearman_age_PROG_table
spearman_age_PROG_table_FDR$FDR <- p.adjust(spearman_age_PROG_table_FDR$P_value, method = "fdr")
# View(spearman_age_PROG_table_FDR)
write.csv(spearman_age_PROG_table_FDR, file="results/spearman_age_PROG_table_FDR.csv")

mw_table_delivery_ANDR_FDR <- mw_table_delivery_ANDR
mw_table_delivery_ANDR_FDR$FDR <- p.adjust(mw_table_delivery_ANDR_FDR$P_value, method = "fdr")
# View(mw_table_delivery_ANDR_FDR)
write.csv(mw_table_delivery_ANDR_FDR, file="results/mw_table_delivery_ANDR_FDR.csv")

mw_table_delivery_CORT_FDR <- mw_table_delivery_CORT
mw_table_delivery_CORT_FDR$FDR <- p.adjust(mw_table_delivery_CORT_FDR$P_value, method = "fdr")
# View(mw_table_delivery_CORT_FDR)
write.csv(mw_table_delivery_CORT_FDR, file="results/mw_table_delivery_CORT_FDR.csv")

mw_table_delivery_ESTR_FDR <- mw_table_delivery_ESTR
mw_table_delivery_ESTR_FDR$FDR <- p.adjust(mw_table_delivery_ESTR_FDR$P_value, method = "fdr")
# View(mw_table_delivery_ESTR_FDR)
write.csv(mw_table_delivery_ESTR_FDR, file="results/mw_table_delivery_ESTR_FDR.csv")

mw_table_delivery_PROG_FDR <- mw_table_delivery_PROG
mw_table_delivery_PROG_FDR$FDR <- p.adjust(mw_table_delivery_PROG_FDR$P_value, method = "fdr")
# View(mw_table_delivery_PROG_FDR)
write.csv(mw_table_delivery_PROG_FDR, file="results/mw_table_delivery_PROG_FDR.csv")

mw_table_breastfeeding_ANDR_FDR <- mw_table_breastfeeding_ANDR
mw_table_breastfeeding_ANDR_FDR$FDR <- p.adjust(mw_table_breastfeeding_ANDR_FDR$P_value, method = "fdr")
# View(mw_table_breastfeeding_ANDR_FDR)
write.csv(mw_table_breastfeeding_ANDR_FDR, file="results/mw_table_breastfeeding_ANDR_FDR.csv")

mw_table_breastfeeding_PROG_FDR <- mw_table_breastfeeding_PROG
mw_table_breastfeeding_PROG_FDR$FDR <- p.adjust(mw_table_breastfeeding_PROG_FDR$P_value, method = "fdr")
# View(mw_table_breastfeeding_PROG_FDR)
write.csv(mw_table_breastfeeding_PROG_FDR, file="results/mw_table_breastfeeding_PROG_FDR.csv")

mw_table_breastfeeding_ESTR_FDR <- mw_table_breastfeeding_ESTR
mw_table_breastfeeding_ESTR_FDR$FDR <- p.adjust(mw_table_breastfeeding_ESTR_FDR$P_value, method = "fdr")
# View(mw_table_breastfeeding_ESTR_FDR)
write.csv(mw_table_breastfeeding_ESTR_FDR, file="results/mw_table_breastfeeding_ESTR_FDR.csv")

mw_table_breastfeeding_PROG_FDR <- mw_table_breastfeeding_PROG
mw_table_breastfeeding_PROG_FDR$FDR <- p.adjust(mw_table_breastfeeding_PROG_FDR$P_value, method = "fdr")
# View(mw_table_breastfeeding_PROG_FDR)
write.csv(mw_table_breastfeeding_PROG_FDR, file="results/mw_table_breastfeeding_PROG_FDR.csv")

mw_table_ATD_ANDR_FDR <- mw_table_ATD_ANDR
mw_table_ATD_ANDR_FDR$FDR <- p.adjust(mw_table_ATD_ANDR_FDR$P_value, method = "fdr")
# View(mw_table_ATD_ANDR_FDR)
write.csv(mw_table_ATD_ANDR_FDR, file="results/mw_table_ATD_ANDR_FDR.csv")

mw_table_ATD_CORT_FDR <- mw_table_ATD_CORT
mw_table_ATD_CORT_FDR$FDR <- p.adjust(mw_table_ATD_CORT_FDR$P_value, method = "fdr")
# View(mw_table_ATD_CORT_FDR)
write.csv(mw_table_ATD_CORT_FDR, file="results/mw_table_ATD_CORT_FDR.csv")

mw_table_ATD_ESTR_FDR <- mw_table_ATD_ESTR
mw_table_ATD_ESTR_FDR$FDR <- p.adjust(mw_table_ATD_ESTR_FDR$P_value, method = "fdr")
# View(mw_table_ATD_ESTR_FDR)
write.csv(mw_table_ATD_ESTR_FDR, file="results/mw_table_ATD_ESTR_FDR.csv")

mw_table_ATD_PROG_FDR <- mw_table_ATD_PROG
mw_table_ATD_PROG_FDR$FDR <- p.adjust(mw_table_ATD_PROG_FDR$P_value, method = "fdr")
# View(mw_table_ATD_PROG_FDR)
write.csv(mw_table_ATD_PROG_FDR, file="results/mw_table_ATD_PROG_FDR.csv")

# Pause
Sys.sleep(1)








# Robust lineal regression CORT ----

## BSS-R + maternal age + delivery ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_part = factor(tipus_part_via,
                        levels = c("Vaginal", "C-section"))
  )

results_rlm_bss_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Corticosteroids = c,
    N = nrow(db_join_s_na),
    Beta_BSS = s$coefficients["bss_total", "Estimate"],
    CI_low = ic["bss_total", 1],
    CI_upp = ic["bss_total", 2],
    P_value = s$coefficients["bss_total", "Pr(>|t|)"]
  )
})

rlm_bss_table_CORT <- do.call(rbind, results_rlm_bss_CORT)

rlm_bss_table_CORT_arranged <- rlm_bss_table_CORT %>%
  filter(!is.na(P_value)) %>%  
  arrange(P_value)

write.csv(rlm_bss_table_CORT_arranged, file="results/rlm_bss_table_CORT_arranged.csv")

rlm_bss_table_CORT_FDR <- rlm_bss_table_CORT_arranged
rlm_bss_table_CORT_FDR$FDR <- p.adjust(rlm_bss_table_CORT_FDR$P_value, method = "fdr")
# View(rlm_bss_table_CORT_FDR)
write.csv(rlm_bss_table_CORT_FDR, file="results/rlm_bss_table_CORT_FDR.csv")

# Forest plot
df_plot <- rlm_bss_table_CORT_arranged %>%
  filter(!is.na(Beta_BSS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_BSS) %>%
  mutate(
    Corticosteroids = factor(Corticosteroids, levels = Corticosteroids),
    
    label = paste0(
      sprintf("%.6f", Beta_BSS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_bss_CORT <- ggplot(df_plot, aes(x = Beta_BSS, y = Corticosteroids)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Corticosteroids) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_BSS-R (95% CI)",
    y = "Corticosteroids",
    title = "Association between BSS-R and corticosteroids"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_bss_CORT.png", forest_plot_bss_CORT, width = 18, height = 7, dpi = 300)

# VIF
vif_list_bss_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Corticosteroids = c,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_bss_table_CORT <- do.call(rbind, vif_list_bss_CORT)

vif_bss_table_CORT_arranged <- vif_bss_table_CORT %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## Brazelton Scale + maternal age + type of delivery ----
results_rlm_braz_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Corticosteroids = c,
    N = nrow(db_join_s_na),
    Beta_braz = s$coefficients["braz_total", "Estimate"],
    CI_low = ic["braz_total", 1],
    CI_upp = ic["braz_total", 2],
    P_value = s$coefficients["braz_total", "Pr(>|t|)"]
  )
})

rlm_braz_table_CORT <- do.call(rbind, results_rlm_braz_CORT)

rlm_braz_table_CORT_arranged <- rlm_braz_table_CORT %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(rlm_braz_table_CORT_arranged, file="results/rlm_braz_table_CORT_arranged.csv")

rlm_braz_table_CORT_FDR <- rlm_braz_table_CORT_arranged
rlm_braz_table_CORT_FDR$FDR <- p.adjust(rlm_braz_table_CORT_FDR$P_value, method = "fdr")
# View(rlm_braz_table_CORT_FDR)
write.csv(rlm_braz_table_CORT_FDR, file="results/rlm_braz_table_CORT_FDR.csv")

# Forest plot
df_plot <- rlm_braz_table_CORT_arranged %>%
  filter(!is.na(Beta_braz), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_braz) %>%
  mutate(
    Corticosteroids = factor(Corticosteroids, levels = Corticosteroids),
    
    label = paste0(
      sprintf("%.6f", Beta_braz), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_braz_CORT <- ggplot(df_plot, aes(x = Beta_braz, y = Corticosteroids)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Corticosteroids) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_Brazelton Scale (95% CI)",
    y = "Corticosteroids",
    title = "Association between Brazelton Scale and corticosteroids"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_braz_CORT.png", forest_plot_braz_CORT, width = 18, height = 7, dpi = 300)


# VIF
vif_list_braz_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Corticosteroids = c,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_braz_CORT <- do.call(rbind, vif_list_braz_CORT)

table_vif_braz_CORT_arranged <- table_vif_braz_CORT %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## PBQ-16 + maternal age + breastfeeding ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_lact = factor(lact_0,
                        levels = c("Breastfeeding", "Formula"))
  )

results_rlm_pbq_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ pbq_total + edat_1avis + grupo_lact")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Corticosteroids = c,
    N = nrow(db_join_s_na),
    Beta_PBQ = s$coefficients["pbq_total", "Estimate"],
    CI_low = ic["pbq_total", 1],
    CI_upp = ic["pbq_total", 2],
    P_value = s$coefficients["pbq_total", "Pr(>|t|)"]
  )
})

table_rlm_pbq_CORT <- do.call(rbind, results_rlm_pbq_CORT)

table_rlm_pbq_CORT_arranged <- table_rlm_pbq_CORT %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_pbq_CORT_arranged, file="results/table_rlm_pbq_CORT_arranged.csv")

table_rlm_pbq_CORT_FDR <- table_rlm_pbq_CORT_arranged
table_rlm_pbq_CORT_FDR$FDR <- p.adjust(table_rlm_pbq_CORT_FDR$P_value, method = "fdr")
# View(table_rlm_pbq_CORT_FDR)
write.csv(table_rlm_pbq_CORT_FDR, file="results/table_rlm_pbq_CORT_FDR.csv")

# Forest plot
df_plot <- table_rlm_pbq_CORT_arranged %>%
  filter(!is.na(Beta_PBQ), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_PBQ) %>%
  mutate(
    Corticosteroids = factor(Corticosteroids, levels = Corticosteroids),
    
    label = paste0(
      sprintf("%.6f", Beta_PBQ), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_pbq_CORT <- ggplot(df_plot, aes(x = Beta_PBQ, y = Corticosteroids)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Corticosteroids) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_PBQ-16 (95% CI)",
    y = "Corticosteroids",
    title = "Association between PBQ-16 and corticosteroids"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_pbq_CORT.png", forest_plot_pbq_CORT, width = 18, height = 7, dpi = 300)

# VIF
vif_list_pbq_CORT <- lapply(CORT, function(c) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(c, "~ pbq_total + edat_1avis + grupo_lact")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  
  if(is.null(modelo)) {
    return(data.frame(
      Corticosteroids = c,
      Variable = c("pbq_total", "edat_1avis", "pbq_total"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("pbq_total", "edat_1avis", "pbq_total")
  
  
  data.frame(
    corticosteroids = c,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_pbq_CORT <- do.call(rbind, vif_list_pbq_CORT)

table_vif_pbq_CORT_arranged <- table_vif_pbq_CORT %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## EPDS + maternal age + ATD ----
ctrl <- lmrob.control(
  maxit.scale = 500,
  max.it = 500,
  k.max = 2000
)

results_rlm_epds_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ epds_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na,
    control = ctrl
  )

  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Corticosteroids = c,
    N = nrow(db_join_s_na),
    Beta_EPDS = s$coefficients["epds_total", "Estimate"],
    CI_low = ic["epds_total", 1],
    CI_upp = ic["epds_total", 2],
    P_value = s$coefficients["epds_total", "Pr(>|t|)"]
  )
})

table_rlm_epds_CORT <- do.call(rbind, results_rlm_epds_CORT)

table_rlm_epds_CORT_arranged <- table_rlm_epds_CORT %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_epds_CORT_arranged, file="results/table_rlm_epds_CORT_arranged.csv")

table_rlm_epds_CORT_FDR <- table_rlm_epds_CORT_arranged
table_rlm_epds_CORT_FDR$FDR <- p.adjust(table_rlm_epds_CORT_FDR$P_value, method = "fdr")
# View(table_rlm_epds_CORT_FDR)
write.csv(table_rlm_epds_CORT_FDR, file="results/table_rlm_epds_CORT_FDR.csv")

# Forest plot
df_plot <- table_rlm_epds_CORT_arranged %>%
  filter(!is.na(Beta_EPDS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_EPDS) %>%
  mutate(
    Corticosteroids = factor(Corticosteroids, levels = Corticosteroids),
    
    label = paste0(
      sprintf("%.6f", Beta_EPDS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_epds_CORT <- ggplot(df_plot, aes(x = Beta_EPDS, y = Corticosteroids)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Corticosteroids) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_EPDS (95% CI)",
    y = "Corticosteroids",
    title = "Association between EPDS and corticosteroids"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_epds_CORT.png", forest_plot_epds_CORT, width = 18, height = 7, dpi = 300)

# VIF
vif_list_epds_CORT <- lapply(CORT, function(c) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(c, "~ epds_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Corticosteroids = c,
      Variable = c("epds_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("epds_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    Corticosteroids = c,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_epds_CORT <- do.call(rbind, vif_list_epds_CORT)

table_vif_epds_CORT_arranged <- table_vif_epds_CORT %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)

## HDRS + maternal age + ATD ----
results_rlm_hamd_CORT <- lapply(CORT, function(c) {
  
  modelo <- lmrob(
    as.formula(paste(c, "~ hamd_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Corticosteroids = c,
    N = nrow(db_join_s_na),
    Beta_HDRS = s$coefficients["hamd_total", "Estimate"],
    CI_low = ic["hamd_total", 1],
    CI_upp = ic["hamd_total", 2],
    P_value = s$coefficients["hamd_total", "Pr(>|t|)"]
  )
})

table_rlm_hamd_CORT <- do.call(rbind, results_rlm_hamd_CORT)

table_rlm_hamd_CORT_arranged <- table_rlm_hamd_CORT %>%
  filter(!is.na(P_value)) %>%        
  arrange(P_value)

write.csv(table_rlm_hamd_CORT_arranged, file="results/table_rlm_hamd_CORT_arranged.csv")

table_rlm_hamd_CORT_FDR <- table_rlm_hamd_CORT_arranged
table_rlm_hamd_CORT_FDR$FDR <- p.adjust(table_rlm_hamd_CORT_FDR$P_value, method = "fdr")
# View(table_rlm_hamd_CORT_FDR)
write.csv(table_rlm_hamd_CORT_FDR, file="results/table_rlm_hamd_CORT_FDR.csv")

# Forest plot
df_plot <- table_rlm_hamd_CORT_arranged %>%
  filter(!is.na(Beta_HDRS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_HDRS) %>%
  mutate(
    Corticosteroids = factor(Corticosteroids, levels = Corticosteroids),
    
    label = paste0(
      sprintf("%.6f", Beta_HDRS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_hamd_CORT <- ggplot(df_plot, aes(x = Beta_HDRS, y = Corticosteroids)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Corticosteroids) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_HDRS (95% CI)",
    y = "Corticosteroids",
    title = "Association between HDRS and corticosteroids"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_hamd_CORT.png", forest_plot_hamd_CORT, width = 18, height = 7, dpi = 300)

# VIF
vif_list_hamd_CORT <- lapply(CORT, function(c) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(c, "~ hamd_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Corticosteroids = c,
      Variable = c("hamd_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("hamd_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    corticosteroids = c,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_hamd_CORT <- do.call(rbind, vif_list_hamd_CORT)

table_vif_hamd_CORT_arranged <- table_vif_hamd_CORT %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)






# Robust lineal regression ANDR ----

## BSS-R + maternal age + delivery ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_part = factor(tipus_part_via,
                        levels = c("Vaginal", "C-section"))
  )

results_rlm_bss_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Androgens = a,
    N = nrow(db_join_s_na),
    Beta_BSS = s$coefficients["bss_total", "Estimate"],
    CI_low = ic["bss_total", 1],
    CI_upp = ic["bss_total", 2],
    P_value = s$coefficients["bss_total", "Pr(>|t|)"]
  )
})

rlm_bss_table_ANDR <- do.call(rbind, results_rlm_bss_ANDR)

rlm_bss_table_ANDR_arranged <- rlm_bss_table_ANDR %>%
  filter(!is.na(P_value)) %>%  
  arrange(P_value)

write.csv(rlm_bss_table_ANDR_arranged, file="results/rlm_bss_table_ANDR_arranged.csv")

rlm_bss_table_ANDR_FDR <- rlm_bss_table_ANDR_arranged
rlm_bss_table_ANDR_FDR$FDR <- p.adjust(rlm_bss_table_ANDR_FDR$P_value, method = "fdr")
# View(rlm_bss_table_ANDR_FDR)
write.csv(rlm_bss_table_ANDR_FDR, file="results/rlm_bss_table_ANDR_FDR.csv")

# Forest plot
df_plot <- rlm_bss_table_ANDR_arranged %>%
  filter(!is.na(Beta_BSS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_BSS) %>%
  mutate(
    Androgens = factor(Androgens, levels = Androgens),
    
    label = paste0(
      sprintf("%.6f", Beta_BSS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_bss_ANDR <- ggplot(df_plot, aes(x = Beta_BSS, y = Androgens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Androgens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.5) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_BSS-R (95% CI)",
    y = "Androgens",
    title = "Association between BSS-R and androgens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_bss_ANDR.png", forest_plot_bss_ANDR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_bss_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Androgens = a,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_bss_table_ANDR <- do.call(rbind, vif_list_bss_ANDR)

vif_bss_table_ANDR_arranged <- vif_bss_table_ANDR %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## Brazelton Scale + maternal age + type of delivery ----
results_rlm_braz_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Androgens = a,
    N = nrow(db_join_s_na),
    Beta_braz = s$coefficients["braz_total", "Estimate"],
    CI_low = ic["braz_total", 1],
    CI_upp = ic["braz_total", 2],
    P_value = s$coefficients["braz_total", "Pr(>|t|)"]
  )
})

rlm_braz_table_ANDR <- do.call(rbind, results_rlm_braz_ANDR)

rlm_braz_table_ANDR_arranged <- rlm_braz_table_ANDR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(rlm_braz_table_ANDR_arranged, file="results/rlm_braz_table_ANDR_arranged.csv")

rlm_braz_table_ANDR_FDR <- rlm_braz_table_ANDR_arranged
rlm_braz_table_ANDR_FDR$FDR <- p.adjust(rlm_braz_table_ANDR_FDR$P_value, method = "fdr")
# View(rlm_braz_table_ANDR_FDR)
write.csv(rlm_braz_table_ANDR_FDR, file="results/rlm_braz_table_ANDR_FDR.csv")

# Forest plot
df_plot <- rlm_braz_table_ANDR_arranged %>%
  filter(!is.na(Beta_braz), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_braz) %>%
  mutate(
    Androgens = factor(Androgens, levels = Androgens),
    
    label = paste0(
      sprintf("%.6f", Beta_braz), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_braz_ANDR <- ggplot(df_plot, aes(x = Beta_braz, y = Androgens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Androgens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_Brazelton Scale (95% CI)",
    y = "Androgens",
    title = "Association between Brazelton Scale and androgens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_braz_ANDR.png", forest_plot_braz_ANDR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_braz_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Androgens = a,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_braz_ANDR <- do.call(rbind, vif_list_braz_ANDR)

table_vif_braz_ANDR_arranged <- table_vif_braz_ANDR %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## PBQ-16 + maternal age + breastfeeding ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_lact = factor(lact_0,
                        levels = c("Breastfeeding", "Formula"))
  )

results_rlm_pbq_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ pbq_total + edat_1avis + grupo_lact")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Androgens = a,
    N = nrow(db_join_s_na),
    Beta_PBQ = s$coefficients["pbq_total", "Estimate"],
    CI_low = ic["pbq_total", 1],
    CI_upp = ic["pbq_total", 2],
    P_value = s$coefficients["pbq_total", "Pr(>|t|)"]
  )
})

table_rlm_pbq_ANDR <- do.call(rbind, results_rlm_pbq_ANDR)

table_rlm_pbq_ANDR_arranged <- table_rlm_pbq_ANDR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_pbq_ANDR_arranged, file="results/table_rlm_pbq_ANDR_arranged.csv")

table_rlm_pbq_ANDR_FDR <- table_rlm_pbq_ANDR_arranged
table_rlm_pbq_ANDR_FDR$FDR <- p.adjust(table_rlm_pbq_ANDR_FDR$P_value, method = "fdr")
# View(table_rlm_pbq_ANDR_FDR)
write.csv(table_rlm_pbq_ANDR_FDR, file="results/table_rlm_pbq_ANDR_FDR.csv")

# Forest plot
df_plot <- table_rlm_pbq_ANDR_arranged %>%
  filter(!is.na(Beta_PBQ), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_PBQ) %>%
  mutate(
    Androgens = factor(Androgens, levels = Androgens),
    
    label = paste0(
      sprintf("%.6f", Beta_PBQ), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_pbq_ANDR <- ggplot(df_plot, aes(x = Beta_PBQ, y = Androgens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Androgens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_PBQ-16 (95% CI)",
    y = "Androgens",
    title = "Association between PBQ-16 and androgens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_pbq_ANDR.png", forest_plot_pbq_ANDR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_pbq_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ pbq_total + edat_1avis + grupo_lact")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Androgens = a,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_pbq_ANDR <- do.call(rbind, vif_list_pbq_ANDR)

table_vif_pbq_ANDR_arranged <- table_vif_pbq_ANDR %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## EPDS + maternal age + ATD ----
ctrl <- lmrob.control(
  maxit.scale = 500,
  max.it = 500,
  k.max = 2000
)

results_rlm_epds_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ epds_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na,
    control = ctrl
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Androgens = a,
    N = nrow(db_join_s_na),
    Beta_EPDS = s$coefficients["epds_total", "Estimate"],
    CI_low = ic["epds_total", 1],
    CI_upp = ic["epds_total", 2],
    P_value = s$coefficients["epds_total", "Pr(>|t|)"]
  )
})

table_rlm_epds_ANDR <- do.call(rbind, results_rlm_epds_ANDR)

table_rlm_epds_ANDR_arranged <- table_rlm_epds_ANDR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_epds_ANDR_arranged, file="results/table_rlm_epds_ANDR_arranged.csv")

table_rlm_epds_ANDR_FDR <- table_rlm_epds_ANDR_arranged
table_rlm_epds_ANDR_FDR$FDR <- p.adjust(table_rlm_epds_ANDR_FDR$P_value, method = "fdr")
# View(table_rlm_epds_ANDR_FDR)
write.csv(table_rlm_epds_ANDR_FDR, file="results/table_rlm_epds_ANDR_FDR.csv")

# Forest plot
df_plot <- table_rlm_epds_ANDR_arranged %>%
  filter(!is.na(Beta_EPDS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_EPDS) %>%
  mutate(
    Androgens = factor(Androgens, levels = Androgens),
    
    label = paste0(
      sprintf("%.6f", Beta_EPDS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_epds_ANDR <- ggplot(df_plot, aes(x = Beta_EPDS, y = Androgens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Androgens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_EPDS (95% CI)",
    y = "Androgens",
    title = "Association between EPDS and androgens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_epds_ANDR.png", forest_plot_epds_ANDR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_epds_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(a, "~ epds_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Androgens = a,
      Variable = c("epds_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("epds_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    androgens = a,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_epds_ANDR <- do.call(rbind, vif_list_epds_ANDR)

table_vif_epds_ANDR_arranged <- table_vif_epds_ANDR %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)

## HDRS + maternal age + ATD ----
results_rlm_hamd_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- lmrob(
    as.formula(paste(a, "~ hamd_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Androgens = a,
    N = nrow(db_join_s_na),
    Beta_HDRS = s$coefficients["hamd_total", "Estimate"],
    CI_low = ic["hamd_total", 1],
    CI_upp = ic["hamd_total", 2],
    P_value = s$coefficients["hamd_total", "Pr(>|t|)"]
  )
})

table_rlm_hamd_ANDR <- do.call(rbind, results_rlm_hamd_ANDR)

table_rlm_hamd_ANDR_arranged <- table_rlm_hamd_ANDR %>%
  filter(!is.na(P_value)) %>%        
  arrange(P_value)

write.csv(table_rlm_hamd_ANDR_arranged, file="results/table_rlm_hamd_ANDR_arranged.csv")

table_rlm_hamd_ANDR_FDR <- table_rlm_hamd_ANDR_arranged
table_rlm_hamd_ANDR_FDR$FDR <- p.adjust(table_rlm_hamd_ANDR_FDR$P_value, method = "fdr")
# View(table_rlm_hamd_ANDR_FDR)
write.csv(table_rlm_hamd_ANDR_FDR, file="results/table_rlm_hamd_ANDR_FDR.csv")

# Forest plot
df_plot <- table_rlm_hamd_ANDR_arranged %>%
  filter(!is.na(Beta_HDRS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_HDRS) %>%
  mutate(
    Androgens = factor(Androgens, levels = Androgens),
    
    label = paste0(
      sprintf("%.6f", Beta_HDRS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_hamd_ANDR <- ggplot(df_plot, aes(x = Beta_HDRS, y = Androgens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Androgens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_HDRS (95% CI)",
    y = "Androgens",
    title = "Association between HDRS and androgens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_hamd_ANDR.png", forest_plot_hamd_ANDR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_hamd_ANDR <- lapply(ANDR, function(a) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(a, "~ hamd_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Androgens = a,
      Variable = c("hamd_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("hamd_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    androgens = a,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_hamd_ANDR <- do.call(rbind, vif_list_hamd_ANDR)

table_vif_hamd_ANDR_arranged <- table_vif_hamd_ANDR %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)






# Robust lineal regression ESTR ----

## BSS-R + maternal age + delivery ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_part = factor(tipus_part_via,
                        levels = c("Vaginal", "C-section"))
  )

results_rlm_bss_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Estrogens = e,
    N = nrow(db_join_s_na),
    Beta_BSS = s$coefficients["bss_total", "Estimate"],
    CI_low = ic["bss_total", 1],
    CI_upp = ic["bss_total", 2],
    P_value = s$coefficients["bss_total", "Pr(>|t|)"]
  )
})

rlm_bss_table_ESTR <- do.call(rbind, results_rlm_bss_ESTR)

rlm_bss_table_ESTR_arranged <- rlm_bss_table_ESTR %>%
  filter(!is.na(P_value)) %>%  
  arrange(P_value)

write.csv(rlm_bss_table_ESTR_arranged, file="results/rlm_bss_table_ESTR_arranged.csv")

rlm_bss_table_ESTR_FDR <- rlm_bss_table_ESTR_arranged
rlm_bss_table_ESTR_FDR$FDR <- p.adjust(rlm_bss_table_ESTR_FDR$P_value, method = "fdr")
# View(rlm_bss_table_ESTR_FDR)
write.csv(rlm_bss_table_ESTR_FDR, file="results/rlm_bss_table_ESTR_FDR.csv")

# Forest plot
df_plot <- rlm_bss_table_ESTR_arranged %>%
  filter(!is.na(Beta_BSS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_BSS) %>%
  mutate(
    Estrogens = factor(Estrogens, levels = Estrogens),
    
    label = paste0(
      sprintf("%.6f", Beta_BSS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_bss_ESTR <- ggplot(df_plot, aes(x = Beta_BSS, y = Estrogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Estrogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  
  expand_limits(x = x_text * 1.2) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_BSS-R (95% CI)",
    y = "Estrogens",
    title = "Association between BSS-R and estrogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_bss_ESTR.png", forest_plot_bss_ESTR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_bss_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Estrogens = e,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_bss_table_ESTR <- do.call(rbind, vif_list_bss_ESTR)

vif_bss_table_ESTR_arranged <- vif_bss_table_ESTR %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## Brazelton Scale + maternal age + type of delivery ----
results_rlm_braz_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Estrogens = e,
    N = nrow(db_join_s_na),
    Beta_braz = s$coefficients["braz_total", "Estimate"],
    CI_low = ic["braz_total", 1],
    CI_upp = ic["braz_total", 2],
    P_value = s$coefficients["braz_total", "Pr(>|t|)"]
  )
})

rlm_braz_table_ESTR <- do.call(rbind, results_rlm_braz_ESTR)

rlm_braz_table_ESTR_arranged <- rlm_braz_table_ESTR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(rlm_braz_table_ESTR_arranged, file="results/rlm_braz_table_ESTR_arranged.csv")

rlm_braz_table_ESTR_FDR <- rlm_braz_table_ESTR_arranged
rlm_braz_table_ESTR_FDR$FDR <- p.adjust(rlm_braz_table_ESTR_FDR$P_value, method = "fdr")
# View(rlm_braz_table_ESTR_FDR)
write.csv(rlm_braz_table_ESTR_FDR, file="results/rlm_braz_table_ESTR_FDR.csv")

# Forest plot
df_plot <- rlm_braz_table_ESTR_arranged %>%
  filter(!is.na(Beta_braz), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_braz) %>%
  mutate(
    Estrogens = factor(Estrogens, levels = Estrogens),
    
    label = paste0(
      sprintf("%.6f", Beta_braz), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_braz_ESTR <- ggplot(df_plot, aes(x = Beta_braz, y = Estrogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Estrogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_Brazelton Scale (95% CI)",
    y = "Estrogens",
    title = "Association between Brazelton Scale and estrogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_braz_ESTR.png", forest_plot_braz_ESTR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_braz_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Estrogens = e,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_braz_ESTR_table <- do.call(rbind, vif_list_braz_ESTR)

vif_braz_ESTR_table_arranged <- vif_braz_ESTR_table %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## PBQ-16 + maternal age + breastfeeding ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_lact = factor(lact_0,
                        levels = c("Breastfeeding", "Formula"))
  )

results_rlm_pbq_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ pbq_total + edat_1avis + grupo_lact")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Estrogens = e,
    N = nrow(db_join_s_na),
    Beta_PBQ = s$coefficients["pbq_total", "Estimate"],
    CI_low = ic["pbq_total", 1],
    CI_upp = ic["pbq_total", 2],
    P_value = s$coefficients["pbq_total", "Pr(>|t|)"]
  )
})

table_rlm_pbq_ESTR <- do.call(rbind, results_rlm_pbq_ESTR)

table_rlm_pbq_ESTR_arranged <- table_rlm_pbq_ESTR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_pbq_ESTR_arranged, file="results/table_rlm_pbq_ESTR_arranged.csv")

table_rlm_pbq_ESTR_FDR <- table_rlm_pbq_ESTR_arranged
table_rlm_pbq_ESTR_FDR$FDR <- p.adjust(table_rlm_pbq_ESTR_FDR$P_value, method = "fdr")
# View(table_rlm_pbq_ESTR_FDR)
write.csv(table_rlm_pbq_ESTR_FDR, file="results/table_rlm_pbq_ESTR_FDR.csv")

# Forest plot
df_plot <- table_rlm_pbq_ESTR_arranged %>%
  filter(!is.na(Beta_PBQ), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_PBQ) %>%
  mutate(
    Estrogens = factor(Estrogens, levels = Estrogens),
    
    label = paste0(
      sprintf("%.6f", Beta_PBQ), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_pbq_ESTR <- ggplot(df_plot, aes(x = Beta_PBQ, y = Estrogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Estrogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_PBQ-16 (95% CI)",
    y = "Estrogens",
    title = "Association between PBQ-16 and estrogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_pbq_ESTR.png", forest_plot_pbq_ESTR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_pbq_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(e, "~ pbq_total + edat_1avis + grupo_lact")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  # If model fails, return NAs
  if(is.null(modelo)) {
    return(data.frame(
      Estrogens = e,
      Variable = c("pbq_total", "edat_1avis", "pbq_total"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("pbq_total", "edat_1avis", "pbq_total")
  
  
  data.frame(
    estrogens = e,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_pbq_ESTR_table <- do.call(rbind, vif_list_pbq_ESTR)

vif_pbq_ESTR_table_arranged <- vif_pbq_ESTR_table %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## EPDS + maternal age + ATD ----
ctrl <- lmrob.control(
  maxit.scale = 500,
  max.it = 500,
  k.max = 2000
)

results_rlm_epds_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ epds_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na,
    control = ctrl
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Estrogens = e,
    N = nrow(db_join_s_na),
    Beta_EPDS = s$coefficients["epds_total", "Estimate"],
    CI_low = ic["epds_total", 1],
    CI_upp = ic["epds_total", 2],
    P_value = s$coefficients["epds_total", "Pr(>|t|)"]
  )
})

table_rlm_epds_ESTR <- do.call(rbind, results_rlm_epds_ESTR)

table_rlm_epds_ESTR_arranged <- table_rlm_epds_ESTR %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_epds_ESTR_arranged, file="results/table_rlm_epds_ESTR_arranged.csv")

table_rlm_epds_ESTR_FDR <- table_rlm_epds_ESTR_arranged
table_rlm_epds_ESTR_FDR$FDR <- p.adjust(table_rlm_epds_ESTR_FDR$P_value, method = "fdr")
# View(table_rlm_epds_ESTR_FDR)
write.csv(table_rlm_epds_ESTR_FDR, file="results/table_rlm_epds_ESTR_FDR.csv")

# Forest plot
df_plot <- table_rlm_epds_ESTR_arranged %>%
  filter(!is.na(Beta_EPDS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_EPDS) %>%
  mutate(
    Estrogens = factor(Estrogens, levels = Estrogens),
    
    label = paste0(
      sprintf("%.6f", Beta_EPDS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_epds_ESTR <- ggplot(df_plot, aes(x = Beta_EPDS, y = Estrogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Estrogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_EPDS (95% CI)",
    y = "Estrogens",
    title = "Association between EPDS and estrogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_epds_ESTR.png", forest_plot_epds_ESTR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_epds_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(e, "~ epds_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Estrogens = e,
      Variable = c("epds_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("epds_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    estrogens = e,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_epds_ESTR_table <- do.call(rbind, vif_list_epds_ESTR)

vif_epds_ESTR_table_arranged <- vif_epds_ESTR_table %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)

## HDRS + maternal age + ATD ----
results_rlm_hamd_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- lmrob(
    as.formula(paste(e, "~ hamd_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Estrogens = e,
    N = nrow(db_join_s_na),
    Beta_HDRS = s$coefficients["hamd_total", "Estimate"],
    CI_low = ic["hamd_total", 1],
    CI_upp = ic["hamd_total", 2],
    P_value = s$coefficients["hamd_total", "Pr(>|t|)"]
  )
})

table_rlm_hamd_ESTR <- do.call(rbind, results_rlm_hamd_ESTR)

table_rlm_hamd_ESTR_arranged <- table_rlm_hamd_ESTR %>%
  filter(!is.na(P_value)) %>%        
  arrange(P_value)

write.csv(table_rlm_hamd_ESTR_arranged, file="results/table_rlm_hamd_ESTR_arranged.csv")

table_rlm_hamd_ESTR_FDR <- table_rlm_hamd_ESTR_arranged
table_rlm_hamd_ESTR_FDR$FDR <- p.adjust(table_rlm_hamd_ESTR_FDR$P_value, method = "fdr")
# View(table_rlm_hamd_ESTR_FDR)
write.csv(table_rlm_hamd_ESTR_FDR, file="results/table_rlm_hamd_ESTR_FDR.csv")

# Forest plot
df_plot <- table_rlm_hamd_ESTR_arranged %>%
  filter(!is.na(Beta_HDRS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_HDRS) %>%
  mutate(
    Estrogens = factor(Estrogens, levels = Estrogens),
    
    label = paste0(
      sprintf("%.6f", Beta_HDRS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_hamd_ESTR <- ggplot(df_plot, aes(x = Beta_HDRS, y = Estrogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Estrogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_HDRS (95% CI)",
    y = "Estrogens",
    title = "Association between HDRS and estrogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  )

ggsave("results/forest_plot_hamd_ESTR.png", forest_plot_hamd_ESTR, width = 18, height = 7, dpi = 300)

# VIF
vif_list_hamd_ESTR <- lapply(ESTR, function(e) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(e, "~ hamd_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(z) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Estrogens = e,
      Variable = c("hamd_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(z) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("hamd_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    estrogens = e,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_hamd_ESTR_table <- do.call(rbind, vif_list_hamd_ESTR)

vif_hamd_ESTR_table_arranged <- vif_hamd_ESTR_table %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)





# Robust lineal regression PROG ----

## BSS-R + maternal age + delivery ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_part = factor(tipus_part_via,
                        levels = c("Vaginal", "C-section"))
  )

results_rlm_bss_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Progestogens = p,
    N = nrow(db_join_s_na),
    Beta_BSS = s$coefficients["bss_total", "Estimate"],
    CI_low = ic["bss_total", 1],
    CI_upp = ic["bss_total", 2],
    P_value = s$coefficients["bss_total", "Pr(>|t|)"]
  )
})

rlm_bss_table_PROG <- do.call(rbind, results_rlm_bss_PROG)

rlm_bss_table_PROG_arranged <- rlm_bss_table_PROG %>%
  filter(!is.na(P_value)) %>%  
  arrange(P_value)

write.csv(rlm_bss_table_PROG_arranged, file="results/rlm_bss_table_PROG_arranged.csv")

rlm_bss_table_PROG_FDR <- rlm_bss_table_PROG_arranged
rlm_bss_table_PROG_FDR$FDR <- p.adjust(rlm_bss_table_PROG_FDR$P_value, method = "fdr")
# View(rlm_bss_table_PROG_FDR)
write.csv(rlm_bss_table_PROG_FDR, file="results/rlm_bss_table_PROG_FDR.csv")

# Forest plot
df_plot <- rlm_bss_table_PROG_arranged %>%
  filter(!is.na(Beta_BSS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_BSS) %>%
  mutate(
    Progestogens = factor(Progestogens, levels = Progestogens),
    
    label = paste0(
      sprintf("%.6f", Beta_BSS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05
x_limite <- x_text * 1.6

forest_plot_bss_PROG <- ggplot(df_plot, aes(x = Beta_BSS, y = Progestogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Progestogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_BSS-R (95% CI)",
    y = "Progestogens",
    title = "Association between BSS-R and progestogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    # plot.margin = margin(5.5, 100, 5.5, 5.5)
    plot.margin = margin(5.5, 140, 5.5, 5.5)
  )

ggsave("results/forest_plot_bss_PROG.png", forest_plot_bss_PROG, width = 18, height = 7, dpi = 300)

# VIF
vif_list_bss_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ bss_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Progestogens = p,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_bss_table_PROG <- do.call(rbind, vif_list_bss_PROG)

vif_bss_table_PROG_arranged <- vif_bss_table_PROG %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## Brazelton Scale + maternal age + type of delivery ----
results_rlm_braz_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Progestogens = p,
    N = nrow(db_join_s_na),
    Beta_braz = s$coefficients["braz_total", "Estimate"],
    CI_low = ic["braz_total", 1],
    CI_upp = ic["braz_total", 2],
    P_value = s$coefficients["braz_total", "Pr(>|t|)"]
  )
})

rlm_braz_table_PROG <- do.call(rbind, results_rlm_braz_PROG)

rlm_braz_table_PROG_arranged <- rlm_braz_table_PROG %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(rlm_braz_table_PROG_arranged, file="results/rlm_braz_table_PROG_arranged.csv")

rlm_braz_table_PROG_FDR <- rlm_braz_table_PROG_arranged
rlm_braz_table_PROG_FDR$FDR <- p.adjust(rlm_braz_table_PROG_FDR$P_value, method = "fdr")
# View(rlm_braz_table_PROG_FDR)
write.csv(rlm_braz_table_PROG_FDR, file="results/rlm_braz_table_PROG_FDR.csv")

# Forest plot
df_plot <- rlm_braz_table_PROG_arranged %>%
  filter(!is.na(Beta_braz), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_braz) %>%
  mutate(
    Progestogens = factor(Progestogens, levels = Progestogens),
    
    label = paste0(
      sprintf("%.6f", Beta_braz), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_braz_PROG <- ggplot(df_plot, aes(x = Beta_braz, y = Progestogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Progestogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_Brazelton Scale (95% CI)",
    y = "Progestogens",
    title = "Association between Brazelton Scale and progestogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
    # plot.margin = margin(5.5, 140, 5.5, 5.5)
  )

ggsave("results/forest_plot_braz_PROG.png", forest_plot_braz_PROG, width = 18, height = 7, dpi = 300)

# VIF
vif_list_braz_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ braz_total + edat_1avis + grupo_part")),
    data = db_join_s_na
  )
  
  vif_vals <- vif(modelo)
  
  data.frame(
    Progestogens = p,
    Variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_braz_PROG_table <- do.call(rbind, vif_list_braz_PROG)

vif_braz_PROG_table_arranged <- vif_braz_PROG_table %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## PBQ-16 + maternal age + breastfeeding ----
db_join_s_na <- db_join_s_na %>%
  mutate(
    grupo_lact = factor(lact_0,
                        levels = c("Breastfeeding", "Formula"))
  )

results_rlm_pbq_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ pbq_total + edat_1avis + grupo_lact")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(e) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Progestogens = p,
    N = nrow(db_join_s_na),
    Beta_PBQ = s$coefficients["pbq_total", "Estimate"],
    CI_low = ic["pbq_total", 1],
    CI_upp = ic["pbq_total", 2],
    P_value = s$coefficients["pbq_total", "Pr(>|t|)"]
  )
})

table_rlm_pbq_PROG <- do.call(rbind, results_rlm_pbq_PROG)

table_rlm_pbq_PROG_arranged <- table_rlm_pbq_PROG %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_pbq_PROG_arranged, file="results/table_rlm_pbq_PROG_arranged.csv")

table_rlm_pbq_PROG_FDR <- table_rlm_pbq_PROG_arranged
table_rlm_pbq_PROG_FDR$FDR <- p.adjust(table_rlm_pbq_PROG_FDR$P_value, method = "fdr")
# View(table_rlm_pbq_PROG_FDR)
write.csv(table_rlm_pbq_PROG_FDR, file="results/table_rlm_pbq_PROG_FDR.csv")

# Forest plot
df_plot <- table_rlm_pbq_PROG_arranged %>%
  filter(!is.na(Beta_PBQ), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_PBQ) %>%
  mutate(
    Progestogens = factor(Progestogens, levels = Progestogens),
    
    label = paste0(
      sprintf("%.6f", Beta_PBQ), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_pbq_PROG <- ggplot(df_plot, aes(x = Beta_PBQ, y = Progestogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Progestogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +

  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_PBQ-16 (95% CI)",
    y = "Progestogens",
    title = "Association between PBQ-16 and progestogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
    # plot.margin = margin(5.5, 140, 5.5, 5.5)
  )

ggsave("results/forest_plot_pbq_PROG.png", forest_plot_pbq_PROG, width = 18, height = 7, dpi = 300)

# VIF
vif_list_pbq_PROG <- lapply(PROG, function(p) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(p, "~ pbq_total + edat_1avis + grupo_lact")),
          data = db_join_s_na),
    error = function(e) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Progestogens = p,
      Variable = c("pbq_total", "edat_1avis", "pbq_total"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(e) rep(NA_real_, 3) # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("pbq_total", "edat_1avis", "pbq_total")
  
  
  data.frame(
    progestogens = p,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_pbq_PROG_table <- do.call(rbind, vif_list_pbq_PROG)

vif_pbq_PROG_table_arranged <- vif_pbq_PROG_table %>%
  filter(!is.na(VIF)) %>%
  arrange(VIF)

# Pause
Sys.sleep(1)

## EPDS + maternal age + ATD ----
ctrl <- lmrob.control(
  maxit.scale = 500,
  max.it = 500,
  k.max = 2000
)

results_rlm_epds_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ epds_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na,
    control = ctrl
  )
  
  s <- summary(modelo)
  ic <- tryCatch(
    confint(modelo),
    error = function(z) matrix(NA, nrow = length(coef(modelo)), ncol = 2,
                               dimnames = list(names(coef(modelo)), c("2.5 %", "97.5 %")))
  )
  
  data.frame(
    Progestogens = p,
    N = nrow(db_join_s_na),
    Beta_EPDS = s$coefficients["epds_total", "Estimate"],
    CI_low = ic["epds_total", 1],
    CI_upp = ic["epds_total", 2],
    P_value = s$coefficients["epds_total", "Pr(>|t|)"]
  )
})

table_rlm_epds_PROG <- do.call(rbind, results_rlm_epds_PROG)

table_rlm_epds_PROG_arranged <- table_rlm_epds_PROG %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(table_rlm_epds_PROG_arranged, file="results/table_rlm_epds_PROG_arranged.csv")

table_rlm_epds_PROG_FDR <- table_rlm_epds_PROG_arranged
table_rlm_epds_PROG_FDR$FDR <- p.adjust(table_rlm_epds_PROG_FDR$P_value, method = "fdr")
# View(table_rlm_epds_PROG_FDR)
write.csv(table_rlm_epds_PROG_FDR, file="results/table_rlm_epds_PROG_FDR.csv")

# Forest plot
df_plot <- table_rlm_epds_PROG_arranged %>%
  filter(!is.na(Beta_EPDS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_EPDS) %>%
  mutate(
    Progestogens = factor(Progestogens, levels = Progestogens),
    
    label = paste0(
      sprintf("%.6f", Beta_EPDS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_epds_PROG <- ggplot(df_plot, aes(x = Beta_EPDS, y = Progestogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Progestogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_EPDS (95% CI)",
    y = "Progestogens",
    title = "Association between EPDS and progestogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
    # plot.margin = margin(5.5, 140, 5.5, 5.5)
  )

ggsave("results/forest_plot_epds_PROG.png", forest_plot_epds_PROG, width = 18, height = 7, dpi = 300)

# VIF
vif_list_epds_PROG <- lapply(PROG, function(p) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(p, "~ epds_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(e) NULL
  )
  
  if(is.null(modelo)) {
    return(data.frame(
      Progestogens = p,
      Variable = c("epds_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(e) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("epds_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    progestogens = p,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

vif_epds_PROG_table <- do.call(rbind, vif_list_epds_PROG)

vif_epds_PROG_table_arranged <- vif_epds_PROG_table %>% arrange(VIF) %>% filter(!is.na(VIF))

# Pause
Sys.sleep(1)

## HDRS + maternal age + ATD ----
results_rlm_hamd_PROG <- lapply(PROG, function(p) {
  
  modelo <- lmrob(
    as.formula(paste(p, "~ hamd_total + edat_1avis + fcodic.factor")),
    data = db_join_s_na
  )
  
  s <- summary(modelo)
  ic <- confint(modelo)
  
  data.frame(
    Progestogens = p,
    N = nrow(db_join_s_na),
    Beta_HDRS = s$coefficients["hamd_total", "Estimate"],
    CI_low = ic["hamd_total", 1],
    CI_upp = ic["hamd_total", 2],
    P_value = s$coefficients["hamd_total", "Pr(>|t|)"]
  )
})

table_rlm_hamd_PROG <- do.call(rbind, results_rlm_hamd_PROG)

table_rlm_hamd_PROG_arranged <- table_rlm_hamd_PROG %>%
  filter(!is.na(P_value)) %>%        
  arrange(P_value)

write.csv(table_rlm_hamd_PROG_arranged, file="results/table_rlm_hamd_PROG_arranged.csv")

table_rlm_hamd_PROG_FDR <- table_rlm_hamd_PROG_arranged
table_rlm_hamd_PROG_FDR$FDR <- p.adjust(table_rlm_hamd_PROG_FDR$P_value, method = "fdr")
# View(table_rlm_hamd_PROG_FDR)
write.csv(table_rlm_hamd_PROG_FDR, file="results/table_rlm_hamd_PROG_FDR.csv")

# Forest plot
df_plot <- table_rlm_hamd_PROG_arranged %>%
  filter(!is.na(Beta_HDRS), !is.na(CI_low), !is.na(CI_upp)) %>%
  arrange(Beta_HDRS) %>%
  mutate(
    Progestogens = factor(Progestogens, levels = Progestogens),
    
    label = paste0(
      sprintf("%.6f", Beta_HDRS), " (",
      sprintf("%.6f", CI_low), " – ",
      sprintf("%.6f", CI_upp), ")   p=",
      sprintf("%.6f", P_value)
    ),
    
    Significance = ifelse(P_value < 0.05, "< 0.05", "> 0.05")
  )

x_text <- max(df_plot$CI_upp, na.rm = TRUE) * 1.05

forest_plot_hamd_PROG <- ggplot(df_plot, aes(x = Beta_HDRS, y = Progestogens)) +
  
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_upp, color = Significance),
    width  = 0.2,
    orientation = "y"
  ) +
  
  geom_point(aes(color = Significance), size = 2.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_text(
    aes(label = label),
    x = x_text,
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = length(df_plot$Progestogens) + 1,
    label = "Beta (95% CI)                        P-value",
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "grey50")) +
  
  expand_limits(x = x_text * 1.4) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Beta_HDRS (95% CI)",
    y = "Progestogens",
    title = "Association between HDRS and progestogens"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, hjust = 0.41,
                                margin = margin(r = 100, l = 0)),
    legend.position = "bottom",
    legend.justification = c(0.5, 0),
    legend.margin = margin(0, 100, 0, 0),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "grey50"),
    plot.margin = margin(5.5, 100, 5.5, 5.5)
    # plot.margin = margin(5.5, 140, 5.5, 5.5)
  )

ggsave("results/forest_plot_hamd_PROG.png", forest_plot_hamd_PROG, width = 18, height = 7, dpi = 300)

# VIF
vif_list_hamd_PROG <- lapply(PROG, function(p) {
  
  modelo <- tryCatch(
    lmrob(as.formula(paste(p, "~ hamd_total + edat_1avis + fcodic.factor")),
          data = db_join_s_na),
    error = function(e) NULL
  )
  
  
  if(is.null(modelo)) {
    return(data.frame(
      Progestogens = p,
      Variable = c("hamd_total", "edat_1avis", "fcodic.factor"),
      VIF = NA_real_
    ))
  }
  
  vif_vals <- tryCatch(
    vif(modelo),
    error = function(e) rep(NA_real_, 3)  # 3 default variables
  )
  
  # Ensure column names are correct if model fails
  if(is.null(names(vif_vals))) names(vif_vals) <- c("hamd_total", "edat_1avis", "fcodic.factor")
  
  
  data.frame(
    progestogens = p,
    variable = names(vif_vals),
    VIF = as.numeric(vif_vals)
  )
})

table_vif_hamd_PROG <- do.call(rbind, vif_list_hamd_PROG)

table_vif_hamd_PROG_arranged <- table_vif_hamd_PROG %>% arrange(VIF) %>% filter(!is.na(VIF))






# Table 1. Sample characteristics of the included participants ----

db_join_s_na$sex_baby <- dplyr::recode(
  db_join_s_na$sex_baby,
  "0" = "Female",
  "1" = "Male"
)

db_join_s_na$ethnicity <- dplyr::recode(
  db_join_s_na$ethnicity,
  "1" = "Caucasian",
  "7" = "Hispanic",
  "6" = "Caribbean",
  "2" = "Gypsy"
)

db_join_s_na$educ_madre <- dplyr::recode(
  db_join_s_na$educ_madre,
  "0" = "Low",
  "1" = "High"
)

db_join_s_na$lact_0 <- factor(
  db_join_s_na$lact_0,
  levels = c("Breastfeeding", "Formula")
)

db_join_s_na$tipus_part_via <- factor(
  db_join_s_na$tipus_part_via,
  levels = c("Vaginal", "C-section")
)

db_join_s_na$fcodic.factor <- factor(
  db_join_s_na$fcodic.factor,
  levels = c("No ATD use", "ATD use at some point")
)

# Mean (SD)
mean_sd_sep <- (function(x) {
  c(
    main = sprintf("%.2f", mean(x, na.rm = TRUE)),
    spread = sprintf("%.2f", sd(x, na.rm = TRUE))
  )
})

# Median (IQR)
median_iqr_sep <- (function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  c(
    main = sprintf("%.2f", q[2]),
    spread = paste0("[", sprintf("%.2f", q[1]), " – ", sprintf("%.2f", q[3]), "]")
  )
})

# n (%)
n_percent_sep <- (function(x) {
  tbl <- table(x)
  prop <- prop.table(tbl) * 100
  
  c(
    main = paste0(
      names(tbl), ": ", tbl, " (", sprintf("%.1f", prop), "%)",
      collapse = "; "
    ),
    spread = "-"
  )
})

db_join_s_na$tpal_4_v2 <- as.integer(db_join_s_na$tpal_4_v2)
db_join_s_na$gest_setm <- as.numeric(db_join_s_na$gest_setm)
table1 <- data.frame(
  
  Variable = c(
    "N",
    "Age, years",
    "Gestational age, weeks",
    "Infant sex",
    "Number of children",
    "Ethnicity",
    "Education level",
    "Type of breastfeeding",
    "Type of delivery",
    "Antidepressant use",
    "BSS-R",
    "PBQ-16",
    "HDRS",
    "EPDS",
    "Brazelton Scale"
  ),
  
  Value = c(
    nrow(db_join_s_na),
    mean_sd_sep(db_join_s_na$edat_1avis)["main"],
    mean_sd_sep(db_join_s_na$gest_setm)["main"],
    n_percent_sep(db_join_s_na$sex_baby)["main"],
    median_iqr_sep(db_join_s_na$tpal_4_v2)["main"],
    n_percent_sep(db_join_s_na$ethnicity)["main"],
    n_percent_sep(db_join_s_na$educ_madre)["main"],
    n_percent_sep(db_join_s_na$lact_0)["main"],
    n_percent_sep(db_join_s_na$tipus_part_via)["main"],
    n_percent_sep(db_join_s_na$fcodic.factor)["main"],
    mean_sd_sep(db_join_s_na$bss_total)["main"],
    median_iqr_sep(db_join_s_na$pbq_total)["main"],
    median_iqr_sep(db_join_s_na$hamd_total)["main"],
    median_iqr_sep(db_join_s_na$epds_total)["main"],
    median_iqr_sep(db_join_s_na$braz_total)["main"]
  ),
  
  Spread = c(
    "-",
    mean_sd_sep(db_join_s_na$edat_1avis)["spread"],
    mean_sd_sep(db_join_s_na$gest_setm)["spread"],
    "-",
    median_iqr_sep(db_join_s_na$tpal_4_v2)["spread"],
    "-",
    "-",
    "-",
    "-",
    "-",
    mean_sd_sep(db_join_s_na$bss_total)["spread"],
    median_iqr_sep(db_join_s_na$pbq_total)["spread"],
    median_iqr_sep(db_join_s_na$hamd_total)["spread"],
    median_iqr_sep(db_join_s_na$epds_total)["spread"],
    median_iqr_sep(db_join_s_na$braz_total)["spread"]
  )
)


table1 %>%
  kable(
    caption = "Table 1. Sample characteristics of the included participants.",
    align = "lccc"
  ) %>%
  kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE
  )





# Table 2. Associations between subjective childbirth experience (BSS-R) and steroid hormone profiles across major endocrine pathways: Robust linear regression models. ----
a <- rlm_bss_table_CORT_arranged %>% rename(Hormone = Corticosteroids)
b <- rlm_bss_table_ANDR_arranged %>% rename(Hormone = Androgens)
c <- rlm_bss_table_ESTR_arranged %>% rename(Hormone = Estrogens)
d <- rlm_bss_table_PROG_arranged %>% rename(Hormone = Progestogens)
df_total <- dplyr::bind_rows(a,b,c,d)

df_total %>%
  kable(
    caption = "Table 2. Associations between subjective childbirth experience (BSS-R) and steroid hormone profiles across major endocrine pathways: Robust linear regression models.",
    align = "lccc"
  ) %>%
  kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE
  )

df_total_arranged <- df_total %>% arrange(P_value)

df_total_arranged %>%
  kable(
    caption = "Table 2. Associations between subjective childbirth experience (BSS-R) and steroid hormone profiles across major endocrine pathways: Robust linear regression models.",
    align = "lccc"
  ) %>%
  kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE
  )
