# Examples to NMetC early career night


# Load necessary libraries
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(42)

# Parameters
n_samples <- 10
n_features <- 30000
n_significant <- 200

# Create column names for samples
sample_names <- c(paste0("heart_", 1:n_samples), 
                  paste0("liver_", 1:n_samples), 
                  paste0("blood_", 1:n_samples))

# Initialize the data matrix with random values (normal distribution)
data <- matrix(rnorm(n_features * length(sample_names), mean = 0, sd = 1), 
               nrow = n_features, 
               ncol = length(sample_names))

# Assign row names (feature names)
rownames(data) <- paste0("feature_", 1:n_features)

# Add signal to the first 200 features to make them statistically significant
# Heart will have a higher mean in the first 100 features
data[1:100, 1:n_samples] <- data[1:100, 1:n_samples] + rnorm(100, mean = 5, sd = 1)

# Liver will have a higher mean in the next 100 features
data[101:200, (n_samples + 1):(2 * n_samples)] <- data[101:200, (n_samples + 1):(2 * n_samples)] + rnorm(100, mean = 5, sd = 1)

# Convert the data matrix to a data frame
df <- as.data.frame(data)
colnames(df) <- sample_names


# Verify with ANOVA to confirm the number of statistically significant features
anova_results <- apply(df, 1, function(x) {
  fit <- aov(x ~ factor(rep(c("heart", "liver", "blood"), each = n_samples)))
  summary(fit)[[1]][["Pr(>F)"]][1]
})



# Count how many features have p-value < 0.05
significant_features <- sum(anova_results < 0.05)
cat("Number of significant features (p < 0.05):", significant_features, "\n")

# Output the first few rows of the data frame
head(df)

################################################################################
################################################################################
################################################################################

# Pca and heatmap

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)

# Perform PCA on the data
pca_res <- prcomp(t(df), scale. = TRUE)

# Plot the PCA
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], 
                     Sample = factor(rep(c("heart", "liver", "blood"), each = n_samples)))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Metabolomics Data", x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")


################
################
#Suffle the samples in illustration purpose

# Create a vector to store the original sample types
sample_types <- c(rep("heart", n_samples), rep("liver", n_samples), rep("blood", n_samples))

# Shuffle the columns while keeping track of the order
shuffled_indices <- sample(ncol(df))
df <- df[, shuffled_indices]
shuffled_sample_names <- sample_names[shuffled_indices]
shuffled_sample_types <- sample_types[shuffled_indices]

# Print the shuffled column names and their corresponding sample types
print(data.frame(Shuffled_Column = shuffled_sample_names, Sample_Type = shuffled_sample_types))



# Take 50 features using ANOVA
top_features <- names(sort(anova_results)[1:50])

# Subset the data frame to include only the top 50 significant features
top_df <- df[top_features, ]

# Z-score normalization for better visualization
z_score_df <- t(scale(t(top_df)))

# Set colors for each sample type
sample_colors <- c("heart" = "#64a85c", "liver" = "#b75bf3", "blood" = "#b88f3e")

# Map the shuffled sample types to their respective colors
shuffled_colors <- sample_colors[shuffled_sample_types]




# Create a heatmap using ComplexHeatmap
Heatmap(z_score_df,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("#0099c8", "#ffffff", "#ff6400")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(Sample = factor(shuffled_sample_types, 
                                                           levels = c("heart", "liver", "blood")),
                                           col = list(Sample = shuffled_colors)))
        # top_annotation = HeatmapAnnotation(Sample = factor(rep(c("heart", "liver", "blood"), each = n_samples)),
        #                                    col = list(Sample = sample_colors)))

# No clustering
Heatmap(z_score_df,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("#0099c8", "#ffffff", "#ff6400")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(Sample = factor(shuffled_sample_types, 
                                                           levels = c("heart", "liver", "blood")),
                                           col = list(Sample = shuffled_colors)))
        # top_annotation = HeatmapAnnotation(Sample = factor(rep(c("heart", "liver", "blood"), each = n_samples)),
        #                                    col = list(Sample = sample_colors)))



##############################################################################
##############################################################################
##############################################################################
# With notame data

# Install the notame package (if not already installed)
# install.packages("remotes")
# remotes::install_github("boennemannlab/notame")

# Install ComplexHeatmap if not installed
# install.packages("ComplexHeatmap")

# Load the necessary libraries
library(notame)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


# Load the example dataset from notame
data("example_set")

# Extract the data and metadata
data_matrix <- t(exprs(example_set))  # Transpose to have samples as columns
metadata <- pData(example_set)

# Perform PCA
pca_res <- prcomp(data_matrix)

# Prepare PCA plot data
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], 
                     Time = metadata$Time)

# Plot the PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Time)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Example Set Data", x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")

##################
##################
##################
# Or easier way
plot_pca(drop_qcs(example_set))
plot_pca(drop_qcs(example_set), color = "Time")

# ANOVA to find significant features
anova_results <- apply(data_matrix, 2, function(x) {
  fit <- aov(x ~ metadata$Time)
  summary(fit)[[1]][["Pr(>F)"]][1]
})

# Find top 50 significant features
top_features <- names(sort(anova_results))

# Subset the data for top 50 features
top_data <- data_matrix[, top_features]

# Z-score normalization for better visualization
z_score_data <- t(scale(top_data))

# Define colors for annotations
group_colors <- c("1" = "red", "2" = "blue", "QC" = "green")

# Create a heatmap using ComplexHeatmap
Heatmap(z_score_data,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(Group = metadata$Time,
                                           col = list(Group = group_colors)))


####################################
####################################
# Heatmap

# Find top significant features
top_features <- names(sort(anova_results, decreasing = T))
top_features <- head(top_features, n= 15)

# Subset the data for top 50 features
top_data <- data_matrix[, top_features]

# Z-score normalization for better visualization
z_score_data <- t(scale(top_data))

group_colors <- c("A" = "#64a85c", "B" = "#b75bf3", "QC" = "#b88f3e")
Heatmap(z_score_data,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("#0099c8", "#ffffff", "#ff6400")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(Group = metadata$Group,
                                           col = list(Group = group_colors)
                                           )
        )

# Some boxplots from the data
library(ggpubr)

# Select the data that you want to boxplot
feat_selected <- head(top_features, n = 4)
selected_metaboset <- example_set[example_set@featureData@data$Feature_ID %in% feat_selected,]

# Selfmade ggpubr extensions

# plot_list <- save_group_boxplots(object = drop_qcs(selected_metaboset),
#                              save = F, all_features = F, color = "Time")

# notame_boxplot(example_set, separate = "Time",feat_name = "HILIC_pos_108_1065a2_6121", free_y = T)
# example_set$Time <- factor(example_set$Time, levels = c("2", "QC", "1"))
# notame_boxplot(example_set, separate = "Time",feat_name = "HILIC_pos_446_9413a2_315", free_y = T)

ggarrange(plotlist = plot_list, labels = "auto")
# make long data to be boxplotted

## One way to make boxplot

metabolite <- as.data.frame(t(exprs(drop_qcs(example_set)[fData(drop_qcs(example_set))$Feature_ID == "HILIC_pos_108_1065a2_6121"])))
metabolite <- cbind(metabolite, pData(drop_qcs(example_set))[,"Time"])

colnames(metabolite) <- c("metabol", "sep")
metabolite$sep <- as.factor(metabolite$sep)

# Start from zero!
ggpubr::ggboxplot(metabolite, x = "sep", y = "metabol", fill = "sep", palette = c("#00AFBB", "#E7B800"),
                 add = c("mean_sd"), ylim = c(10000, 55000))

