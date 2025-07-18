library(dplyr)
library(data.table)
library(plotly)

x <- list.files("/blue/guralnick/millerjared/SuperBlooms/data/processed/", pattern = "singularAdditive", full.names = TRUE)

store <- list()
for(i in 1:length(x)){
  store[[i]] <- readRDS(x[i])
}

results_all <- bind_rows(store)
fwrite(results_all, "/blue/guralnick/millerjared/SuperBlooms/data/processed/simulation-predictions-additive.csv")
results_all <- fread("/blue/guralnick/millerjared/SuperBlooms/data/processed/simulation-predictions-additive.csv")
# Summarize accuracy per parameter combination (too complex for ease of visuals)
accuracy_summary_full <- results_all %>%
  mutate(correct_prediction = ifelse(S_y_site == S_y_site_predict, TRUE, FALSE)) %>%
  group_by(param_set) %>%
  summarize(
    correct_pred = sum(correct_prediction == TRUE),
    incorrect_pred = sum(correct_prediction == FALSE),
    accuracy = (correct_pred/(correct_pred + incorrect_pred)),
    pdetwday = unique(pdetwday),
    pdetwend = unique(pdetwend),
    pdetminyear = unique(min_det),
    pdetmaxyear = unique(max_det),
    pdetlowspatial = unique(low_sdet),
    pdethighspatial = unique(high_sdet),
    .groups = "drop"
  )

# summarize by means (for simpler views)
accuracy_summary <- results_all %>% 
  mutate(correct_prediction = ifelse(S_y_site == S_y_site_predict, TRUE, FALSE)) %>% 
  group_by(param_set) %>% 
  summarize(
    correct_pred = sum(correct_prediction == TRUE), 
    incorrect_pred = sum(correct_prediction == FALSE), 
    accuracy = (correct_pred/(correct_pred + incorrect_pred)),
    pdet_day_avg = (unique(pdetwday) + unique(pdetwend))/2 ,
    pdet_year_avg = (unique(min_det) + unique(max_det))/2 ,
    pdet_spatial_avg = (unique(low_sdet) + unique(high_sdet))/2 ,
    .groups = "drop"
  )


fig <- plot_ly(
  data = accuracy_summary,
  x = ~pdet_day_avg,
  y = ~pdet_year_avg,
  z = ~pdet_spatial_avg,
  color = ~accuracy,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5),
  text = ~paste(
    "Accuracy:", round(accuracy, 2),
    "<br>pdet_day:", pdet_day_avg,
    "<br>pdet_year:", pdet_year_avg,
    "<br>pdet_spatial:", pdet_spatial_avg
  )
) %>%
  layout(
    scene = list(
      xaxis = list(title = "pdet_day"),
      yaxis = list(title = "pdet_year"),
      zaxis = list(title = "pdet_spatial")
    ),
    title = "3D Detection Probability Grid vs Accuracy"
  )
fig

library(tidyr)
library(scatterplot3d)

# Bin accuracy
accuracy_summary$acc_bin <- cut(accuracy_summary$accuracy, breaks = c(0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1),
                                labels = c('60-70%', '70-80%', '80-90%', '90-95%', '95-99%', '99-100%'), include.lowest=TRUE)

palette <- c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850")
point_colors <- palette[accuracy_summary$acc_bin]

s3d <- scatterplot3d(accuracy_summary$pdet_day_avg, accuracy_summary$pdet_year_avg, accuracy_summary$pdet_spatial_avg,
                     pch = 19, color = point_colors,
                     xlab = "pdet_day_avg", ylab = "pdet_year_avg", zlab = "pdet_spatial")
legend("right", legend = levels(accuracy_summary$acc_bin), col = palette, pch = 19, title = "Accuracy Bin")


point_colors <- setNames(palette, levels(df$acc_bin))
df$color <- point_colors[as.character(df$acc_bin)]

library(rgl)
# Open 3D plot
plot3d(
  df$pdet_day, df$pdet_year, df$pdet_spatial,
  col = df$color,
  size = 8,
  type = "s",
  xlab = "pdet_day",
  ylab = "pdet_year",
  zlab = "pdet_spatial"
)

# Add a legend
legend3d("topright", legend = levels(df$acc_bin), pch = 16, col = palette, cex=1.2, title = "Accuracy Bin")


## try a thin plate interp 
library(predicts)
# Prediction as before
X <- as.matrix(accuracy_summary[, c("pdet_day_avg", "pdet_year_avg", "pdet_spatial_avg")])
y <- accuracy_summary$accuracy

tps_model <- Tps(X, y)

# Define grid
gx <- seq(min(X[,1]), max(X[,1]), length = 30)
gy <- seq(min(X[,2]), max(X[,2]), length = 30)
gz <- seq(min(X[,3]), max(X[,3]), length = 30)

# Ensure cube ordering
grid_df <- expand.grid(
  pdet_day_avg = gx,
  pdet_year_avg = gy,
  pdet_spatial_avg = gz
)

pred_vals <- predictSurface(tps_model)
surface(pred_vals)
# Check if any values fall in your isomin:isomax range
range(pred_vals, na.rm = TRUE)

library(fields)
X <- as.matrix(accuracy_summary[, c("pdet_day_avg", "pdet_year_avg", "pdet_spatial_avg")])
y <- accuracy_summary$accuracy
tps_model <- Tps(X, y)
gx <- seq(min(X[,1]), max(X[,1]), length = 50)
gy <- seq(min(X[,2]), max(X[,2]), length = 50)
grid2d <- expand.grid(
  pdet_day_avg = gx,
  pdet_year_avg = gy
)

# Step 1: Find the global min/max for all slices you intend to plot
slice_seq <- seq(0, 1, by = 0.2)
pred_matrix <- sapply(slice_seq, function(s) {
  grid2d$pdet_spatial_avg <- s
  predict(tps_model, as.matrix(grid2d))
})
global_min <- min(pred_matrix, na.rm=TRUE)
global_max <- max(pred_matrix, na.rm=TRUE)

# Step 2: When plotting, use `zlim = c(global_min, global_max)`
for (i in seq_along(slice_seq)) {
  s <- slice_seq[i]
  z_matrix <- matrix(pred_matrix[,i], nrow = length(gx), ncol = length(gy))
  filled.contour(
    gx, gy, z_matrix,
    zlim = c(global_min, global_max),
    plot.title = title(
      main = paste("Accuracy Surface at pdet_spatial_avg =", round(s,2)),
      xlab = "pdet_day_avg",
      ylab = "pdet_year_avg"
    ),
    key.title = title(main = "accuracy")
  )
}
