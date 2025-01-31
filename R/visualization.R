#' Generate Pie Plot Data Frame
#'
#' This function generates a dsata frame for creating pie plots, where each row corresponds to a spatial cell with information about cell types, weights, and coordinates. The function computes weights for each cell type based on the first and second cell type annotations from RCTD results and supports different modes of analysis.
#'
#' @param rctd An object containing RCTD results, including spatial and cell type data.
#' @param mode A character string specifying the mode of computation. Currently, only `"doublet"` mode is implemented. Available options include:
#'   \itemize{
#'     \item `"doublet"`: Computes weights for first and second cell types for doublet cells.
#'     \item `"all"`: Not yet implemented.
#'     \item `"by_spot_class"`: Not yet implemented.
#'   }
#'
#' @return A data frame containing the following columns:
#'   \itemize{
#'     \item `cell_id`: Unique identifier for each spatial cell.
#'     \item `x`, `y`: Spatial coordinates of the cell.
#'     \item `weight_first_type`, `weight_second_type`: Weights for the first and second cell types.
#'     \item `first_type`, `second_type`: Annotations of the first and second cell types.
#'     \item `w1_larger_w2`: Logical indicator of whether the first type weight is larger than the second type weight.
#'     \item One column for each cell type, containing the computed weight of that cell type for each spatial cell.
#'   }
#'
#' @note Only the `"doublet"` mode is currently implemented. Other modes will result in an error.
#'
#' @importFrom dplyr %>% select
#'
#' @export

get_pieplot_df <- function(
    rctd,
    mode = c("doublet", "all", "by_spot_class")[1]
){

  cell_types <- rctd@cell_type_info[[1]][[2]]
  cell_types_N <- length(cell_types)

  pie_df <- rctd@results$results_df %>%
    select(weight_first_type, weight_second_type, first_type, second_type, w1_larger_w2, x, y)
  pie_df$cell_id <- rownames(pie_df)

  if(mode == "doublet"){
    pie_celltype_weigts <- apply(pie_df, 1, FUN = function(row){
      res <- rep(0, cell_types_N)
      names(res) <- cell_types
      ft <- row["first_type"] %>% unname()
      st <- row["second_type"] %>% unname()
      w1 <- row["weight_first_type"] %>% as.numeric()
      if(is.na(st)){
        res[ft] <- 1
      } else {
        res[ft] <- w1
        res[st] <- 1 - w1
      }

      return(res)

    }) %>% t()
  } else {
    stop("Only doublet mode is implemented so far!")
  }

  pie_celltype_weigts <- as.data.frame(pie_celltype_weigts)
  pie_celltype_weigts$cell_id <- rownames(pie_celltype_weigts)


  # Merge the two data frames by `cell_id`
  merged_pie_df <- merge(
    pie_df,
    pie_celltype_weigts,
    by = "cell_id",
    all.x = TRUE  # Keep all rows from pie_df
  )
  rownames(merged_pie_df) <- merged_pie_df$cell_id
  return(merged_pie_df)
}

#' Crop pie_df Data to a Specific Cell
#'
#' This function filters a data frame to retain only the rows within a specified radius around a given cell.
#'
#' @param df A data frame containing spatial cell information.
#' @param cell_id The unique identifier of the central cell.
#' @param radius A numeric value specifying the radius for cropping (default: 100).
#'
#' @importFrom dplyr %>% filter
#' @return A filtered data frame containing only the cells within the specified radius.
#' @export
#'
crop_pie_df_to_cell <- function(
    df,
    cell_id,
    radius = 100
){
  x_center <- df[cell_id, "x"]
  y_center <- df[cell_id, "y"]
  print(x_center)
  print(y_center)
  result <- df %>% filter(
    (x < (x_center + radius) & x > (x_center - radius) &
       y < (y_center + radius) & y > (y_center - radius))
  )
  return(result)
}

#' Generate a Pie Plot
#'
#' This function creates a pie plot from a given data frame containing spatial cell type proportions.
#'
#' @param pie_df A data frame containing spatial cell type proportions and coordinates.
#' @param cols An optional named vector of colors for cell types.
#' @param pie_scale A numeric value to scale pie sizes (default: 1).
#'
#' @return A ggplot object representing the spatial pie chart.
#'
#' @importFrom scatterpie geom_scatterpie
#' @importFrom ggplot2 ggplot coord_fixed scale_color_manual scale_fill_manual aes theme_void guides
#' @importFrom dplyr %>%

#' @export

plot_pie <- function(pie_df, cols = NULL, pie_scale = 1){

  cell_types <- unique(c(pie_df$first_type,pie_df$second_type))
  cell_types <- cell_types[!is.na(cell_types)]

  p <- ggplot() +
    geom_scatterpie(
      aes(x = x, y = -y, group = cell_id, color = first_type),  # Define the position and grouping
      data = pie_df,
      cols = cell_types ,               # Specify columns for the pie slices
      pie_scale = pie_scale
    )
  if(!is.null(cols)){
    p <- p +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols)
  }

  p <- p +
    coord_fixed() +
    theme_void() +
    guides(color = "none")  # Hide the color legend but keep the fill legend

  return(p)
}


#' Plot Pie Chart by Coordinate Range
#'
#' This function filters the pie data frame by x and y limits and generates a pie plot.
#'
#' @param pie_df A data frame containing spatial cell type proportions and coordinates.
#' @param x_lims A numeric vector specifying x-axis limits.
#' @param y_lims A numeric vector specifying y-axis limits.
#' @param cols An optional named vector of colors for cell types.
#' @param pie_scale A numeric value to scale pie sizes (default: 1).
#'
#' @importFrom dplyr %>% filter
#' @return A ggplot object representing the cropped spatial pie chart.
#' @export

plot_pie_by_coordinates <- function(
    pie_df,
    x_lims, y_lims,
    cols = NULL,
    pie_scale = 1
){

  pie_df_crop <- pie_df %>% filter(x > min(x_lims) & x < max(x_lims) & y > min(y_lims) & y < max(y_lims))
  p <- plot_pie(pie_df_crop, cols = cols, pie_scale = pie_scale)

  return(p)
}

#' Plot Pie Chart Around a Specific Cell
#'
#' This function crops the data around a specified cell within a given radius and generates a pie plot.
#'
#' @param pie_df A data frame containing spatial cell type proportions and coordinates.
#' @param cell_id The unique identifier of the central cell.
#' @param radius A numeric value specifying the radius for cropping (default: 100).
#' @param pie_scale A numeric value to scale pie sizes (default: 1).
#' @param cols An optional named vector of colors for cell types.
#' @param DO_highlight_cell A logical indicating whether to highlight the central cell (default: TRUE).
#'
#' @importFrom ggplot2 annotate
#' @return A ggplot object representing the spatial pie chart around the specified cell.
#' @export

plot_pie_around_cell <- function(
    pie_df,
    cell_id,
    radius = 100,
    pie_scale = 1,
    cols = NULL,
    DO_highlight_cell = TRUE
){

  pie_df_crop <- crop_pie_df_to_cell(pie_df, cell_id = cell_id, radius = radius)
  p <- plot_pie(pie_df_crop, cols = cols, pie_scale = pie_scale)

  if(DO_highlight_cell){
    x_center <- pie_df_crop[cell_id, "x"]
    y_center <- pie_df_crop[cell_id, "y"]
    p <- p +
        annotate("point", x = x_center, y = -y_center, color = "black", size = 1, shape = 8)
  }

  return(p)
}

#' Save Pie Plot to File
#'
#' This function saves a given ggplot object as an image file.
#'
#' @param p A ggplot object representing the pie chart.
#' @param filename A character string specifying the output file name.
#' @param ... Additional arguments passed to `ggsave()`, such as width, height, and dpi.
#'
#' @return Saves the plot as a file with a transparent background.
#' @export

save_pieplot <- function(
    p,
    filename,
    ...
){
  ggsave(filename, plot = p, bg = "transparent", ... )
}

