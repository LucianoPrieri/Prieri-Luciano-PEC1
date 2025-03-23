### 1. Carga de paquetes ----------------------------------------------------
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if (!require("readxl")) install.packages("readxl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")

install.packages("tinytex")
tinytex::install_tinytex()  # Instala TinyTeX

library(SummarizedExperiment)
library(readxl)
library(ggplot2)
library(pheatmap)

### 2. Carga y preparación de datos -----------------------------------------
# Leemos el  archivo de  Excel
file_path <- "Gastric_NMR.xlsx"
raw_data <- read_excel(file_path, sheet = "data")
raw_data
head(raw_data)
# Seleccionar columnas 1-7 (metadatos)
sample_metadata <- raw_data[, 1:7]
sample_metadata

# Renombrar columnas 
colnames(sample_metadata) <- c("Idx", "Date", "Sample_Type", "QC", "Batch", "Order", "Sample_ID")


# Separamos los  componentes
# Columnas  para metabolitos (M1-M129)
expression_matrix <- as.matrix(raw_data[, 8:ncol(raw_data)])  # Columnas 8-136
rownames(expression_matrix) <- raw_data$Sample_id

# Estructura de la matriz
dim(expression_matrix)  
rownames(expression_matrix) 
colnames(expression_matrix)  

# La matriz de expresión tiene 140 filas y 129 columnas, lo que significa que:
# 140 muestras están incluidas en el análisis.
# 129 metabolitos están cuantificados en cada muestra.
# Esto confirma que la matriz tiene la estructura esperada, con muestras en las filas y metabolitos en las columnas.

# Se muestra una lista de 140 muestras, etiquetadas como "sample_1" a "sample_140".

# Se listan 129 metabolitos, etiquetados como "M1" a "M129". Cada columna representa un metabolito cuantificado en todas las muestras
# Si hubiera nombres duplicados o inconsistencias, podría afectar el análisis

# La matriz de expresión contiene 140 muestras y 129 metabolitos, asegurando una estructura adecuada para el análisis. Se verificó que las etiquetas de muestra y metabolitos están correctamente asignadas, lo que permite continuar con los análisis estadísticos sin problemas estructurales




# Se observa que cada metabolito (M1, M2, M3, …, M129) está correctamente listado. La columna "Metabolite_ID" contiene los nombres de los metabolitos, que son las columnas de la matriz expression_matrix.
# Se asignan estos mismos nombres como índices del DataFrame (row.names), permitiendo que cada fila se identifique con su propio metabolito

# Metadatos de rasgos (nombres de metabolitos)
metabolite_metadata <- data.frame(
  Metabolite_ID = colnames(expression_matrix),
  row.names = colnames(expression_matrix)
)
metabolite_metadata

# Se observa que cada metabolito (M1, M2, M3, …, M129) está correctamente listado
# Se generó un DataFrame de metadatos que contiene los identificadores de los metabolitos presentes en la matriz de expresión. Esta estructura permite vincular información adicional y facilita la interpretación de los resultados en análisis posteriores


###  Manejo de Valores Faltantes -------------------------------

#  Calculamos el porcentaje de valores faltantes en la matriz de expresión
na_percentage <- mean(is.na(expression_matrix)) * 100  
cat("Porcentaje de valores NA:", na_percentage, "%\n")

#  Exploramos las primeras filas de la matriz de datos para identificar problemas
head(raw_data[, 7:ncol(raw_data)])

# Revisamos la estructura de los datos
str(raw_data[, 7:ncol(raw_data)])

#  Buscamos valores no numéricos en cada columna (pueden ser caracteres extraños)
sapply(raw_data[, 7:ncol(raw_data)], function(x) sum(!grepl("^-?[0-9.]+$", x)))



# Se detectó que aproximadamente el 5% de los valores están ausentes, lo que puede afectar los análisis posteriores si no se maneja adecuadamente
# Se identificaron valores no numéricos en algunas columnas, lo que indica que puede haber datos mal formateados que deben ser corregidos.




### . Limpieza de Datos --------------------------------------------

# Copiamos el dataset original para trabajar con datos limpios
raw_data_clean <- raw_data

#  Eliminamos espacios en blanco dentro de los valores
raw_data_clean[, 7:ncol(raw_data_clean)] <- apply(raw_data_clean[, 7:ncol(raw_data_clean)], 2, function(x) gsub(" ", "", x))  

#  Reemplazamos comas por puntos (para asegurar formato numérico correcto)
raw_data_clean[, 7:ncol(raw_data_clean)] <- apply(raw_data_clean[, 7:ncol(raw_data_clean)], 2, function(x) gsub(",", ".", x))

#  Reemplazamos valores no numéricos con NA (para evitar errores en la conversión)
raw_data_clean[, 7:ncol(raw_data_clean)] <- apply(raw_data_clean[, 7:ncol(raw_data_clean)], 2, function(x) ifelse(grepl("^[0-9.]+$", x), x, NA))

#  Convertimos todas las columnas a tipo numérico
expression_matrix <- apply(raw_data_clean[, 7:ncol(raw_data_clean)], 2, as.numeric)

#  Eliminamos la columna "Sample_id" de la matriz de expresión, ya que no es numérica
expression_matrix <- expression_matrix[, -1]  

#  Contamos los valores NA restantes
cat("Cantidad de valores NA después de la limpieza:", sum(is.na(expression_matrix)), "\n")


# Se llevó a cabo un proceso de limpieza de datos para garantizar la calidad del análisis. Se eliminaron espacios en blanco, se ajustó el formato numérico y se identificaron valores atípicos, dejando una matriz de expresión puramente numérica. Tras la limpieza, se detectaron 915 valores faltantes, los cuales requieren estrategias adicionales para su manejo antes del análisis exploratorio

###  Imputación de Valores Faltantes ----------------------------------

#  Imputamos los valores faltantes usando la mediana de cada columna (metabolito)
for(i in 1:ncol(expression_matrix)) {
  expression_matrix[is.na(expression_matrix[,i]), i] <- median(expression_matrix[,i], na.rm = TRUE)
}

#  Asignamos los nombres de las muestras a las filas de la matriz de expresión
rownames(expression_matrix) <- raw_data$Sample_id

# Verificamos que los nombres de las muestras se han asignado correctamente
rownames(expression_matrix)

###  Creación del Objeto SummarizedExperiment -------------------------

library(SummarizedExperiment)

# Se creó un objeto SummarizedExperiment utilizando el paquete SummarizedExperiment de R para almacenar y analizar los datos de expresión de metabolitos.
# El objeto contiene dos matrices de expresión:

# raw: Matriz original de expresión de metabolitos (filas = metabolitos, columnas = muestras).
# log2: Matriz transformada en log2 para normalizar los datos y reducir el efecto de valores extremos, sumando 1 antes de realizar la transformación logarítmica.
# Metadatos donde: 
  
#colData: Contiene información sobre las muestras (e.g., identificadores, fechas, condiciones experimentales).
# rowData: Contiene información sobre los metabolitos (e.g., identificador del metabolito).


# Creamos el objeto SummarizedExperiment con los datos procesados
se <- SummarizedExperiment(
  assays = list(
    raw = t(expression_matrix),  # Trasponemos: filas = metabolitos, columnas = muestras
    log2 = log2(t(expression_matrix) + 1)  # Transformación log2 para normalizar los datos
  ),
  colData = sample_metadata,     # Metadatos de las muestras
  rowData = metabolite_metadata  # Metadatos de los metabolitos
)

#  Mostramos la estructura del objeto SummarizedExperiment
se


summary(assay(se, "raw"))  # Resumen de los valores de la matriz de expresión
boxplot(assay(se, "raw"), main = "Distribución de Intensidades por Metabolito", las = 2)


# Este análisis muestra que las distribuciones de los valores de expresión tienen una gran variabilidad, con valores mínimos que pueden ser muy bajos (incluso cercanos a 0), pero con valores máximos que superan los 10,000 en algunas muestras.
# La media de los valores de expresión es considerablemente más alta que la mediana, lo que sugiere la presencia de algunos valores extremos en los datos.


# Guardar objeto
save(se, file = "se_metabolomics.Rda")


###  Análisis de Calidad ---------------------------------------------------

#  Función para calcular el Coeficiente de Variación (CV)
calculate_cv <- function(x) sd(x) / mean(x) * 100  # CV = (Desviación estándar / Media) * 100

# Identificamos las muestras de control de calidad (QC)
qc_samples <- which(colData(se)$QC == 1)

#  Calculamos el CV solo para los controles de calidad
cv_values <- apply(assay(se, "raw")[, qc_samples], 1, calculate_cv)

#  Verificamos si hay valores NA en los CV calculados
cat("Número de valores NA en CV:", sum(is.na(cv_values)), "\n") 

# Número de valores NA en CV: 0.
# Esto indica que todos los metabolitos tienen un coeficiente de variación bien definido.


#  Filtramos los metabolitos con CV ≤ 20% en las muestras QC
se_filtered <- se[!is.na(cv_values) & cv_values <= 20, ]

# El coeficiente de variación (CV) se calculó para cada metabolito en las muestras QC. 
# Se eliminaron aquellos con CV > 20% para asegurar estabilidad

#  Eliminamos muestras con valores NA en los datos transformados (log2)
se_filtered <- se_filtered[, colSums(is.na(assay(se_filtered, "log2"))) == 0]

#  Mostramos la nueva dimensión del objeto filtrado
dim(se_filtered)

# Donde X es la cantidad de metabolitos filtrados y Y la cantidad de muestras sin valores NA.





###   Análisis Exploratorio ------------------------------------------------

# Cargamos las librerías necesarias para el análisis exploratorio
library(ggplot2)
library(FactoMineR)
library(factoextra)

#  Aplicamos un Análisis de Componentes Principales (PCA) sobre los datos log-transformados
pca_data <- prcomp(t(assay(se_filtered, "log2")), scale. = TRUE)

#  Mostramos los resultados del PCA
summary(pca_data)

# Aquí podemos observar qué porcentaje de la variabilidad total explican las primeras componentes.
# Generalmente, las dos primeras explican la mayor parte


# Se realizó un Análisis de Componentes Principales (PCA) para visualizar la variabilidad de los datos.
# La Figura X muestra la distribución de las muestras en los dos primeros componentes principales

# Convertimos datos a data frame
pca_df <- as.data.frame(pca_data$x)
pca_df$Batch <- se_filtered$Batch 

# Graficamos PCA 
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.9) +  # Puntos más visibles
  scale_color_gradient(low = "darkblue", high = "lightblue") +  
  labs(title = "PCA por Lote Experimental", x = "PC1", y = "PC2", color = "Lote") +
  theme_minimal(base_size = 14)  


# El gráfico representa un Análisis de Componentes Principales (PCA) de los datos experimentales, diferenciando los lotes mediante una escala de color. 
# El eje X (PC1) y el eje Y (PC2) corresponden a las dos primeros componentes principales, que explican la mayor variabilidad en los datos.

# El análisis PCA muestra que los lotes experimentales no parecen influir significativamente en la variabilidad explicada por las dos primeras


#  Imputamos NA en log2_data 
log2_data <- assay(se_filtered, "log2")
log2_data[is.na(log2_data)] <- median(log2_data, na.rm = TRUE)
assay(se_filtered, "log2") <- log2_data

#  Recalculamos la matriz de correlación 
sample_cor <- cor(assay(se_filtered, "log2"))  
# Calculamos la  varianza por muestra
variances <- apply(assay(se_filtered, "log2"), 2, var)

# Filtramos muestras con varianza > 0
se_filtered <- se_filtered[, variances > 0]

#  Heatmap de correlación entre muestras
library(pheatmap)


pheatmap(
  sample_cor,
  annotation_col = as.data.frame(colData(se_filtered)[, "Batch", drop = FALSE]),  # <-- ¡Cierra aquí el as.data.frame!
  main = "Matriz de Correlación entre Muestras",
  show_rownames = FALSE,   # Opcional para mejorar visualización
  show_colnames = FALSE
)





# Este heatmap muestra la correlación entre todas las muestras
# Colores más cercanos al azul indican alta correlación entre muestras
# Colores más cercanos al rojo indican baja correlación o posibles errores en la adquisición de datos.

# El análisis de evaluación sugiere que las muestras tienen una alta similitud en general, con algunas diferencias entre ciertos subconjuntos
# No se observa un fuerte efecto de lote en la estructura de evaluación, lo que indica que otros factores pueden estar influyendo más en la variabilidad de los datos





# Boxplot de intensidades por muestra
boxplot(assay(se_filtered, "log2"), 
        main = "Distribución de Intensidades (log2)",
        xlab = "Muestras", 
        ylab = "Intensidad",
        col = ifelse(colData(se_filtered)$QC == 1, "red", "blue"))

# Este boxplot muestra la distribución de las intensidades de cada muestra.
# Las muestras QC están en rojo
# Las muestras experimentales están en azul


# El gráfico presentado muestra un boxplot de la distribución de intensidades en escala log2 para diferentes muestras, etiquetadas como sample_1, sample_33, sample_67, etc.
# Este tipo de visualización es útil para evaluar la calidad de los datos en estudios ómicos, como los de metabolómica o transcriptómica
# Las muestras QC están dentro del rango esperado, lo que sugiere que los datos son técnicamente consistentes




###  7. Normalización --------------------------------------------------------

#  Normalización por suma total para corregir diferencias en la intensidad total
total_intensity <- colSums(assay(se_filtered, "raw"))

# Evitamos problemas de división por cero
total_intensity[total_intensity == 0] <- 1e-6  

# Aplicamos la normalización
assay(se_filtered, "normalized") <- assay(se_filtered, "raw") / total_intensity * 1e6

# Verificamos la nueva matriz normalizada
head(assay(se_filtered, "normalized"))

# Este comando mostrará las primeras filas de la matriz normalizada.
# Los valores ahora están escalados a 1,000,000 (1e6) para hacer comparaciones más fáciles entre muestras.

# Si la normalización es correcta, los valores deberían ser similares entre muestras.

# Tras aplicar normalización por suma total, los valores de intensidad están más homogéneos entre muestras, reduciendo efectos de batch


### 8. Exportación de resultados --------------------------------------------
# Guardar gráficos
ggsave("PCA_plot.pdf", width = 8, height = 6)
pdf("Heatmap.pdf", width = 10, height = 8)
pheatmap(sample_cor)
dev.off()

# Generar reporte de metadatos
write.table(colData(se_filtered), file = "metadata_report.txt", sep = "\t")
