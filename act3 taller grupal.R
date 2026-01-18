set.seed(202506)

library(tidyverse)
library(factoextra)
library(pheatmap)
library(gtsummary)
library(broom)
library(performance)
library(pROC)

#1.	Paso 1 - Carga la base de datos y limpieza

df <- read_csv("Base de datos Los Simpson completa.csv", show_col_types = FALSE)

glimpse(df)
summary(df)

genes <- df %>% select(ADCY3:UHMK1)

colSums(is.na(genes))

df <- df[complete.cases(genes), ]

any(is.na(df %>% select(ADCY3:UHMK1)))

#PASO 2. Evaluación de la normalidad (Shapiro–Wilk)

shapiro_genes <- df %>% select(ADCY3:UHMK1) %>% pivot_longer(cols = everything(),
    names_to = "gen",
    values_to = "expresion") %>%
  group_by(gen) %>%
  summarise(
    p_value = shapiro.test(expresion)$p.value,
    normalidad = ifelse(p_value > 0.05, "Compatible con normalidad", "No normal"), .groups = "drop")

shapiro_genes

#Se evaluó la normalidad de cada variable génica mediante el test de Shapiro–Wilk. Los valores p obtenidos se utilizaron para determinar si las distribuciones podían considerarse normales o no. Con este test se pudo mirar que todas las vaiables génicas tienen un p valor inferior 0.05 entoces no seguen una distibuícion normal

#PASO 3. Análisis de Componentes Principales (PCA)

X <- df %>% select(ADCY3:UHMK1)

pca <- prcomp(X, center = TRUE, scale. = TRUE)
summary(pca)

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 15))

scores_pca <- as.data.frame(pca$x[, 1:6])
colnames(scores_pca) <- paste0("PC", 1:6)
df_pca <- bind_cols(df, scores_pca)

#Se realizó un Análisis de Componentes Principales (PCA) utilizando exclusivamente las variables de expresión génica. Se extrajeron las seis primeras componentes principales, que representan la mayor parte de la variabilidad del conjunto de datos.

#Variables del PCA coloreadas por cos2
fviz_pca_var(pca, col.var = "cos2",
  gradient.cols = c("#00AFBB", "#9acd32", "#B452CD"),
  repel = TRUE)

#Clustering de variables 
fviz_pca_var(pca, kmeans = TRUE, k = 3, repel = TRUE)

#Contribución de variables a las componentes
fviz_contrib(pca, choice = "var", axes = 1, top = 15)
fviz_contrib(pca, choice = "var", axes = 2, top = 15)

#Individuos (PC1–PC2) coloreados por IMC
df_pca <- df_pca %>% mutate(IMC_cat = case_when(
  imc_kg_m2 < 25 ~ "Normal",
  imc_kg_m2 >= 25 & imc_kg_m2 < 30 ~ "Sobrepeso",
  imc_kg_m2 >= 30 ~ "Obesidad"))

grafico_IMC <- fviz_pca_ind(pca, geom.ind = "point",
  col.ind = df_pca$IMC_cat,
  palette = c("Normal" = "#00AFBB", "Sobrepeso" = "#9acd32","Obesidad" = "#B452CD"),
  addEllipses = TRUE,
  legend.title = "IMC")

# Clustering de individuos 
Ind_pca <- pca$x[, 1:2]
set.seed(3)

km_ind <- kmeans(Ind_pca, centers = 3, nstart = 25)

fviz_cluster(km_ind, data = Ind_pca, geom = "point")

# Paso5 Heatmap
# Matriz de correlaciones Spearman genes vs componentes principales
cor_matrix <- cor(
  df_pca %>% select(ADCY3:UHMK1),
  df_pca %>% select(PC1:PC6),
  method = "spearman")

cor_matrix <- cor_matrix[ apply(cor_matrix, 1, function(x) all(is.finite(x))), apply(cor_matrix, 2, function(x) all(is.finite(x)))]

pheatmap(cor_matrix,
  color = colorRampPalette(c("#00AFBB", "white", "#B452CD"))(60),
  clustering_method = "complete",
  fontsize_row = 9,
  fontsize_col= 11,
  main = "Correlación Spearman entre genes y componentes principales")

#Paso 6 – Tabla descriptiva por terciles de los componentes principales
make_terciles <- function(x) {
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)
  cut(x,
    breaks = c(-Inf, q, Inf),
    labels = c("t1", "t2", "t3"),
    include.lowest = TRUE)}

df_pca <- df_pca %>%
  mutate(
    PC1_t = make_terciles(PC1),
    PC2_t = make_terciles(PC2),
    PC3_t = make_terciles(PC3))

tabla_desc_PC1 <- df_pca %>%
  select(PC1_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC1_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tabla_desc_PC2 <- df_pca %>%
  select(PC2_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC2_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tabla_desc_PC3 <- df_pca %>%
  select(PC3_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC3_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tabla_desc_PCs <- tbl_merge(
  tbls = list(tabla_desc_PC1, tabla_desc_PC2, tabla_desc_PC3),
  tab_spanner = c(
    "Terciles PC1",
    "Terciles PC2",
    "Terciles PC3"))

tabla_desc_PCs

#7.	Paso 7 - Regresión logística
#7.1 Justificación:
##Para poder usar regresión logística, la variable debe ser binaria
####IL-6 es una citoquina inflamatoria de interés clínico. Su producción puede provenir tanto de las células tumorales como de células del sistema inmune, niveles elevados se han asociado con progresión tumoral y peor pronóstico clínico.

###Las variables predictoras serán:
#Terciles de las componentes principales (PC1_t, PC2_t, PC3_t)
#Variables de ajuste clínicas/sociodemográficas relevantes.

p_il6 <- median(df_pca$il6_pg_ml, na.rm = TRUE)

df_pca <- df_pca %>%
  mutate(
    il6_cat = ifelse(il6_pg_ml < p_il6, "Baja", "Alta"),
    il6_cat = factor(il6_cat, levels = c("Baja", "Alta")))

table(df_pca$il6_cat)

df_pca <- df_pca %>%
  mutate(
    PC1_t = relevel(factor(PC1_t), ref = "t1"),
    PC2_t = relevel(factor(PC2_t), ref = "t1"),
    PC3_t = relevel(factor(PC3_t), ref = "t1"),
    sexo = factor(sexo))

modelo_log <- glm(
  il6_cat ~ PC1_t + PC2_t + PC3_t + edad_anios + sexo + imc_kg_m2,
  data = df_pca,
  family = binomial)

summary(modelo_log)

install.packages("broom.helpers")
library(broom.helpers)

tabla_logistica <- tbl_regression(
  modelo_log,
  exponentiate = TRUE) %>%
  bold_p(t = 0.05) %>%
  bold_labels()

tabla_logistica

#Se ajustó un modelo de regresión logística binaria para evaluar la asociación entre los terciles de las componentes principales derivadas del PCA y los niveles de IL-6 (alta vs baja), ajustando por edad, sexo e IMC. 
###Los resultados muestran que la edad se asocia de forma significativa con una mayor probabilidad de presentar niveles elevados de IL-6, incrementando el riesgo con el aumento de la edad
###El IMC y el t3 de la PC2_t mostraron una tendencia, sugiriendo un posible efecto sobre el estado inflamatorio. No se observaron asociaciones estadísticamente significativas para las restantes componentes principales ni para el sexo.



