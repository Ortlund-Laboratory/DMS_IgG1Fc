#!/usr/bin/env Rscript

pkgs <- c("readr", "dplyr", "plotly", "htmlwidgets")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

data_dir <- "data"
output_html <- file.path("docs", "index.html")
dir.create("docs", showWarnings = FALSE, recursive = TRUE)

files <- c(
#  FCGR1 = file.path(data_dir, "fcgr1.csv"),
#  FCGR2A = file.path(data_dir, "fcgr2a.csv"),
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
  if (aa %in% c("G", "A", "V", "L", "I", "M", "P")) return("Non-Polart")
  if (aa %in% c("F", "Y", "W")) return("Aromatic")
  "Other"
}

group_colors <- c(
  Negative = "#d73027",
  Positive = "#2c7bb6",
  Polar = "#7b3294",
  `Non-Polart` = "#fdae61",
  Aromatic = "#1a9850"
)

row_order <- c("D","E","K","R","H","C","S","T","N","Q","G","A","V","L","I","M","P","F","Y","W")
group_spans <- list(
  Negative = c(1, 2),
  Positive = c(3, 4),
  Polar = c(5, 10),
  `Non-Polart` = c(11, 17),
  Aromatic = c(18, 20)
)

df <- bind_rows(Map(read_dataset, files, names(files))) %>%
  mutate(
    aa = aa_of(mutation),
    group = vapply(aa, group_map, character(1))
  ) %>%
  filter(aa %in% row_order) %>%
  mutate(aa = factor(aa, levels = row_order))

sites <- sort(unique(df$site))

make_matrix <- function(dat, rows, sites) {
  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(sites),
                dimnames = list(rows, as.character(sites)))
  for (i in seq_len(nrow(dat))) {
    r <- match(as.character(dat$aa[i]), rows)
    c <- match(as.character(dat$site[i]), colnames(mat))
    if (!is.na(r) && !is.na(c)) mat[r, c] <- dat$score[i]
  }
  mat
}

matrices <- lapply(names(files), function(rec) {
  make_matrix(filter(df, receptor == rec), row_order, sites)
})
names(matrices) <- names(files)

zmin <- min(df$score, na.rm = TRUE)
zmax <- max(df$score, na.rm = TRUE)

colorscale <- list(
  list(0, "#d73027"),
  list(0.5, "#f7f7f7"),
  list(1, "#2c7bb6")
)

wedge_panel <- plot_ly() %>%
  add_trace(
    x = 1, y = 1,
    type = "scatter",
    mode = "markers",
    marker = list(size = 0),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  layout(
    xaxis = list(visible = FALSE, range = c(0, 1), fixedrange = TRUE),
    yaxis = list(visible = FALSE, range = c(20.5, 0.5), fixedrange = TRUE),
    shapes = lapply(names(group_spans), function(g) {
      span <- group_spans[[g]]
      list(
        type = "rect",
        xref = "paper",
        yref = "y",
        x0 = 0.05,
        x1 = 0.95,
        y0 = span[1] - 0.5,
        y1 = span[2] + 0.5,
        line = list(width = 0),
        fillcolor = group_colors[[g]],
        layer = "above"
      )
    }),
    annotations = lapply(names(group_spans), function(g) {
      span <- group_spans[[g]]
      list(
        x = 0.50,
        y = mean(span),
        xref = "paper",
        yref = "y",
        text = g,
        showarrow = FALSE,
        font = list(size = 10, color = "black"),
        xanchor = "center",
        yanchor = "middle"
      )
    }),
    margin = list(l = 0, r = 0, t = 0, b = 0)
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
      hovertemplate = paste0(
        "AA: %{y}<br>",
        "Site: %{x}<br>",
        "Score: %{z:.3f}",
        "<extra>", rec, "</extra>"
      ),
      showscale = TRUE
    )
}

buttons <- lapply(seq_along(names(matrices)), function(i) {
  vis <- rep(FALSE, length(matrices))
  vis[i] <- TRUE
  list(method = "restyle", args = list("visible", vis), label = names(matrices)[i])
})

heatmap_panel <- heatmap_base %>%
  layout(
    title = list(text = "Fc variant binding heatmap"),
    xaxis = list(title = "Site"),
    yaxis = list(
      title = "AA",
      autorange = "reversed",
      tickfont = list(size = 9),
      automargin = TRUE
    ),
    updatemenus = list(list(
      type = "dropdown",
      x = 0,
      y = 1.15,
      buttons = buttons
    )),
    margin = list(t = 100, l = 20, r = 40, b = 80)
  )

fig <- subplot(
  wedge_panel,
  heatmap_panel,
  widths = c(0.06, 0.94),
  margin = 0,
  shareY = FALSE,
  nrows = 1
) %>%
  layout(
    title = list(text = "Fc variant binding heatmap"),
    margin = list(t = 100, l = 20, r = 40, b = 80)
  )

saveWidget(fig, file = output_html, selfcontained = FALSE)
message("Wrote interactive heatmap to: ", normalizePath(output_html))
