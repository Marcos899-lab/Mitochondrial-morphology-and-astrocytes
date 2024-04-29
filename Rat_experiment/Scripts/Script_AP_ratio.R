#############################################################################
#############SCRIPT RSTUDIO PARA ORGANIZACION DE DATASET#####################
#############################################################################


#################################
###########Librerias#############
#################################

library(ggplot2)
library(tibble)
library(dplyr)


###############Directorios##############
getwd()
setwd("C:/Users/usuario/Desktop/PP/Master_Biomedicina/Practicas/Imaging/Experimentoss/Rata_tratamientos/")

directory <- "C:/Users/usuario/Desktop/PP/Master_Biomedicina/Practicas/Imaging/Experimentoss/Rata_tratamientos"

# Lista de archivos .csv con los que vamos a trabajar
file_names <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Hacemos una lista vac?a para introducir los resutlados
results_list<- list()

# Bucle para automatizar el proceso
for (file_name in file_names) {
  # Leemos el dataset
  mito_data <- read.csv(file_name)
  
  # Creamos un clasificador
  mito_data$Classification <- NA
  
  #Filtrado
  mito_data <- mito_data[mito_data$Area >= 0.22 & mito_data$Perim. >= 1.5, ]
  
  #Creamos una division de Area entre Perimetro para hacer un ratio y meterlo en una columna vacia
  mito_data$AP <- NA
  
  mito_data$AP <- mito_data$Area/mito_data$Perim.
  
  # thresholds
  puncta_threshold_AP <- 0.130
  rods_threshold_AP<- 0.210
  
  # Clasificacion segun los tresholds en area y perimetro
  mito_data$Classification[ mito_data$AP < puncta_threshold_AP] <- "Puncta"
  mito_data$Classification[ mito_data$AP >= puncta_threshold_AP & mito_data$AP < rods_threshold_AP] <- "Rods"
  mito_data$Classification[ mito_data$AP >= rods_threshold_AP] <- "Network"
  
  # Contamos cada nivel de las Clasificaciones
  level_counts <- table(mito_data$Classification)
  
  # Porcentaje de los conteos
  total_mitos <- sum(level_counts)
  percentages <- round(level_counts / total_mitos * 100, 2)
  
  # Creamos un dataset con los conteos
  counts_df<- data.frame(
    puncta = ifelse("Puncta" %in% names(level_counts), percentages["Puncta"], 0),
    rods = ifelse("Rods" %in% names(level_counts), percentages["Rods"], 0),
    network = ifelse("Network" %in% names(level_counts), percentages["Network"], 0),
    se_puncta = sqrt((percentages["Puncta"] * (100 - percentages["Puncta"])) / total_mitos),
    se_rods = sqrt((percentages["Rods"] * (100 - percentages["Rods"])) / total_mitos),
    se_network = sqrt((percentages["Network"] * (100 - percentages["Network"])) / total_mitos)
  )
    
  
  # guardamos los resultados en la lista vacia
  results_list[[basename(file_name)]] <- cbind(counts_df, total_mitos)
}

# combinamos en un dataframe
combined_results_AP <- do.call(rbind, results_list)# Añadimos la columna llamada Mitos
combined_results_AP <- rownames_to_column(combined_results_AP, var = "Mitos")

#Agrupamos por experimento y sacamos la media de porcentajes

grouped_resultsAP <- combined_results_AP %>%
  mutate(Group = rep(1:ceiling(n()/15), each = 15)[1:n()]) %>%
  group_by(Group) %>%
  summarize(
    mean_puncta = mean(puncta),
    mean_rods = mean(rods),
    mean_network = mean(network),
    se_puncta = sd(puncta) / sqrt(n()),
    se_rods = sd(rods) / sqrt(n()),
    se_network = sd(network) / sqrt(n()),
    T_mitos=sum(total_mitos)
  )

grouped_resultsAP$Group <- c( "DMEM","KH_2H", "KH_4H", "KH+carnitina", "KH+etomoxir", "KH+glucosa", "KH+glutamato", "KH+oligomicina",  "KH+TBOA+glutamato", "KH+TBOA", rep("", nrow(grouped_resultsAP) - 10))

# Guardamos en un .csv
write.csv(combined_results_AP, "combined_mitochondria_counts_AP_ratio.csv", row.names = FALSE)
write.csv(grouped_resultsAP, "grouped_results_AP_ratio.csv", row.names = FALSE)

#######################################################################################
########################ANALISIS ESTADISTICO EXPERIMENTAL##############################
#######################################################################################


###############################ANOVA##############################

#Construimos un dataset para hacer el analisis
combined_results_AP$Mitos[1:15] <- "KH"
combined_results_AP$Mitos[16:30] <- "KH+glucose"
combined_results_AP$Mitos[31:45] <- "KH+oligomicine"

# Dataset para elanalisis Puncta
anova_dataset_Puncta <- data.frame(
  Mitos = combined_results_AP$Mitos,  # Factor
  Puncta = combined_results_AP$puncta  # Response variable (puncta counts)
)
anova_dataset_Puncta$Mitos <- as.factor(anova_dataset_Puncta$Mitos)
boxplot(Puncta ~ Mitos, data = anova_dataset_Puncta)

# Perform ANOVA para Puncta
anova_result_Puncta <- aov(Puncta~ Mitos, data = anova_dataset_Puncta)

# Summary del ANOVA
summary(anova_result_Puncta)

#Validación
tuckey_test_puncta <- TukeyHSD(anova_result_Puncta)
tuckey_test_puncta

# Dataset para elanalisis Rods
anova_dataset_Rods <- data.frame(
  Mitos = combined_results_AP$Mitos,  # Factor
  Rods = combined_results_AP$rods  # Response variable (puncta counts)
)
anova_dataset_Rods$Mitos <- as.factor(anova_dataset_Rods$Mitos)
boxplot(Rods ~ Mitos, data = anova_dataset_Rods)

# ANOVA para Rods
anova_result_Rods <- aov(Rods~ Mitos, data = anova_dataset_Rods)

# Summary del ANOVA
summary(anova_result_Rods)

#Validación
tuckey_test_rods <- TukeyHSD(anova_result_Rods)
tuckey_test_rods

# Dataset para el analisis de Network
anova_dataset_Network <- data.frame(
  Mitos = combined_results_AP$Mitos,  # Factor
  Network = combined_results_AP$network  # Response variable (puncta counts)
)
anova_dataset_Network$Mitos <- as.factor(anova_dataset_Network$Mitos)
boxplot(Network ~ Mitos, data = anova_dataset_Network)

# ANOVA para Puncta
anova_result_Network<- aov(Network ~ Mitos, data = anova_dataset_Network)

# Summary del ANOVA
summary(anova_result_Network)

#Validación
tuckey_test_network <- TukeyHSD(anova_result_Network)
tuckey_test_network

###########################################################################
####################Graficas analisis inicial##############################
###########################################################################


##Puncta
ggplot(combined_results_AP, aes(x = Mitos, y = puncta)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(title = "Counts of Mitochondria Morphologies",
       x = "Mitos",
       y = "Counts (puncta)") +
  theme_minimal()

##Rods
ggplot(combined_results_AP, aes(x = Mitos, y = rods)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(title = "Counts of Mitochondria Morphologies",
       x = "Mitos",
       y = "Counts (rods)") +
  theme_minimal()

##Network
ggplot(combined_results_AP, aes(x = Mitos, y = network)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(title = "Counts of Mitochondria Morphologies",
       x = "Mitos",
       y = "Counts (network)") +
  theme_minimal()

######Plots de porcentajes

###Puncta
ggplot(grouped_resultsAP, aes(x = factor(Group), y = mean_puncta)) +
  geom_bar(stat = "identity", fill = "red") +
  geom_errorbar(aes(ymin= mean_puncta - se_puncta, ymax = mean_puncta + se_puncta),
                width= 0.2) + 
  labs(title = "Mean Percentage of Puncta",
       x = "Group",
       y = "Mean Percentage") +
  theme_minimal()

###Rods
ggplot(grouped_resultsAP, aes(x = factor(Group), y = mean_rods)) +
  geom_bar(stat = "identity", fill = "green") +
  geom_errorbar(aes(ymin= mean_rods - se_rods, ymax = mean_rods + se_rods),
                width= 0.2) + 
  labs(title = "Mean Percentage of Rods",
       x = "Group",
       y = "Mean Percentage") +
  theme_minimal()

###Network
ggplot(grouped_resultsAP, aes(x = factor(Group), y = mean_network)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_errorbar(aes(ymin= mean_network - se_network, ymax = mean_network + se_network),
                width= 0.2) + 
  labs(title = "Mean Percentage of Network",
       x = "Group",
       y = "Mean Percentage") +
  theme_minimal()

