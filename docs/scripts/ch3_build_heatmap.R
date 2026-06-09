#!/usr/bin/env Rscript

pkgs <- c("readr", "dplyr", "plotly", "htmlwidgets", "seqinr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

data_dir <- "data"
output_html <- file.path("docs", "index.html")
dir.create("docs", showWarnings = FALSE, recursive = TRUE)

site_start <- 216
fasta_file <- file.path(data_dir, "Fc_prot.fasta")

fc_seq <- seqinr::read.fasta(file = fasta_file, as.string = TRUE, seqtype = "AA")[[1]]
fc_len <- nchar(fc_seq)
sites <- seq(site_start, length.out = fc_len)

files <- c(
  FCGR2B = file.path(data_dir, "FcgR2b_fucosminusafucos_fucos_minus_afucos_scores.csv"),
  FCGR3AF = file.path(data_dir, "FcgR3a158F_fucosminusafucos_fucos_minus_afucos_scores.csv")
)

read_dataset <- function(path, receptor_name) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      receptor = receptor_name,
      site = as.integer(site),
      mutation = mutation,
      score = as.numeric(scaledscore)
    )
}

aa_of <- function(mut) substr(mut, nchar(mut), nchar(mut))

group_map <- function(aa) {
  if (aa %in% c("D", "E")) return("Negative")
  if (aa %in% c("K", "R")) return("Positive")
  if (aa %in% c("H", "C", "S", "T", "N", "Q")) return("Polar")
  if (aa %in% c("G", "A", "V", "L", "I", "M", "P")) return("Non-Polar")
  if (aa %in% c("F", "Y", "W")) return("Aromatic")
  "Other"
}

group_colors <- c(
  Negative = "#d73027",
  Positive = "#2c7bb6",
  Polar = "#7b3294",
  `Non-Polar` = "#fdae61",
  Aromatic = "#1a9850"
)

row_order <- c("D","E","K","R","H","C","S","T","N","Q","G","A","V","L","I","M","P","F","Y","W")

group_spans <- list(
  Negative = c(1, 2),
  Positive = c(3, 4),
  Polar = c(5, 10),
  `Non-Polar` = c(11, 17),
  Aromatic = c(18, 20)
)

df <- bind_rows(Map(read_dataset, files, names(files))) %>%
  mutate(
    aa = aa_of(mutation),
    group = vapply(aa, group_map, character(1))
  ) %>%
  filter(aa %in% row_order) %>%
  mutate(aa = factor(aa, levels = row_order))

real_vals <- df$score[!is.na(df$score)]
data_min <- min(real_vals, na.rm = TRUE)
data_max <- max(real_vals, na.rm = TRUE)
buffer <- (data_max - data_min) / 9
sentinel <- data_min - buffer

make_matrix <- function(dat, rows, sites, sentinel) {
  mat <- matrix(
    sentinel,
    nrow = length(rows),
    ncol = length(sites),
    dimnames = list(rows, as.character(sites))
  )
  for (i in seq_len(nrow(dat))) {
    r <- match(as.character(dat$aa[i]), rows)
    c <- match(as.character(dat$site[i]), colnames(mat))
    if (!is.na(r) && !is.na(c)) mat[r, c] <- ifelse(is.na(dat$score[i]), sentinel, dat$score[i])
  }
  mat
}

matrices <- lapply(names(files), function(rec) {
  dsub <- filter(df, receptor == rec)
  make_matrix(dsub, row_order, sites, sentinel)
})
names(matrices) <- names(files)

zmin <- sentinel
zmax <- data_max

colorscale <- list(
  list(0.000, "#000000"),
  list(0.001, "#d73027"),
  list(0.500, "#f7f7f7"),
  list(1.000, "#2c7bb6")
)

heatmap_base <- plot_ly()
for (i in seq_along(matrices)) {
  rec <- names(matrices)[i]
  mat <- matrices[[i]]
  heatmap_base <- heatmap_base %>%
    add_trace(
      x = colnames(mat),
      y = rownames(mat),
      z = mat,
      type = "heatmap",
      colorscale = colorscale,
      zmin = zmin,
      zmax = zmax,
      visible = i == 1,
      customdata = ifelse(mat == sentinel, "Missing", sprintf("%.3f", mat)),
      hovertemplate = paste0(
        "AA: %{y}<br>",
        "Site: %{x}<br>",
        "Score: %{customdata}<br>",
        "Receptor: ", rec,
        "<extra></extra>"
      ),
      showscale = TRUE
    )
}

buttons <- lapply(seq_along(names(matrices)), function(i) {
  vis <- rep(FALSE, length(matrices))
  vis[i] <- TRUE
  list(method = "restyle", args = list("visible", vis), label = names(matrices)[i])
})

block_shapes <- lapply(names(group_spans), function(g) {
  span <- group_spans[[g]]
  list(
    type = "rect",
    xref = "paper",
    yref = "y",
    x0 = -0.20,
    x1 = -0.02,
    y0 = span[1] - 0.5,
    y1 = span[2] + 0.5,
    line = list(width = 0),
    fillcolor = group_colors[[g]],
    layer = "above"
  )
})

block_text <- lapply(names(group_spans), function(g) {
  span <- group_spans[[g]]
  list(
    x = -0.11,
    y = mean(span),
    xref = "paper",
    yref = "y",
    text = g,
    showarrow = FALSE,
    font = list(size = 11, color = "black"),
    xanchor = "center",
    yanchor = "middle"
  )
})

fig <- heatmap_base %>%
  layout(
    title = list(text = "Fc variant binding heatmap"),
    xaxis = list(
      title = "Site",
      tickmode = "array",
      tickvals = sites,
      ticktext = as.character(sites),
      showticklabels = FALSE
    ),
    yaxis = list(
      title = "AA",
      autorange = "reversed",
      tickmode = "array",
      tickvals = row_order,
      ticktext = row_order,
      tickfont = list(size = 9),
      automargin = TRUE
    ),
    shapes = block_shapes,
    annotations = block_text,
    updatemenus = list(list(
      type = "dropdown",
      x = 0,
      y = 1.15,
      buttons = buttons
    )),
    margin = list(t = 100, l = 220, r = 40, b = 80)
  )

htmlwidgets::saveWidget(fig, file = output_html, selfcontained = FALSE)
message("Wrote interactive heatmap to: ", normalizePath(output_html))
