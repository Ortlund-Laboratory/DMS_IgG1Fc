library(plotly)
library(dplyr)
library(tidyr)

# -----------------------------
# 1. Load input files
# -----------------------------
files <- list(
  FcgR2b = "data/FcgR2b_fucosminusafucos_fucos_minus_afucos_scores.csv",
  FcgR3a158F = "data/FcgR3a158F_fucosminusafucos_fucos_minus_afucos_scores.csv"
)

df_list <- lapply(names(files), function(rec) {
  d <- read.csv(files[[rec]], stringsAsFactors = FALSE)

  d <- d %>% rename(
    site  = site,
    aa    = mutation,
    score = scaledscore
  )

  d$site <- as.numeric(d$site)
  d$receptor <- rec

  d %>% select(site, aa, score, receptor)
})

df <- bind_rows(df_list)

# -----------------------------
# 2. FASTA
# -----------------------------
read_fasta_simple <- function(file) {
  lines <- readLines(file)
  paste(lines[!grepl("^>", lines)], collapse = "")
}

wt_seq <- read_fasta_simple("data/Fc_prot.fasta")
wt_vec <- strsplit(wt_seq, "")[[1]]

offset <- 215

# -----------------------------
# 3. AA order (YOUR ORDER)
# -----------------------------
aa_levels <- c(
  "D","E",
  "K","R",
  "H","C","S","T","N","Q",
  "G","V","L","I","M","P",
  "F","Y","W"
)

group_colors <- c(
  Negative = "#d73027",
  Positive = "#2166ac",
  Polar    = "#984ea3",
  Nonpolar = "#f0e442",
  Aromatic = "#1b9e77"
)

aa_index <- setNames(seq_along(aa_levels), aa_levels)

group_bounds <- list(
  Negative = c("D","E"),
  Positive = c("K","R"),
  Polar    = c("H","C","S","T","N","Q"),
  Nonpolar = c("G","V","L","I","M","P"),
  Aromatic = c("F","Y","W")
)

sites <- sort(unique(df$site))

df2 <- df %>%
  group_by(receptor) %>%
  complete(site = sites, aa = aa_levels) %>%
  ungroup() %>%
  mutate(
    wt_aa = wt_vec[site - offset],
    status = case_when(
      aa == wt_aa ~ "wt",
      is.na(score) ~ "missing",
      TRUE ~ "score"
    )
  )

bar_data_list <- list()

for (rec in unique(df2$receptor)) {

  d <- df2 %>%
    filter(receptor == rec, status == "score")  # ✅ exclude WT + missing

  site_means <- d %>%
    group_by(site) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

  mu <- mean(site_means$mean_score, na.rm = TRUE)
  sdv <- sd(site_means$mean_score, na.rm = TRUE)

  site_means <- site_means %>%
    mutate(
      color = case_when(
        mean_score > (mu + sdv) ~ "#2166ac",   # blue
        mean_score < (mu - sdv) ~ "#b2182b",   # red
        TRUE ~ "#bdbdbd"                       # grey
      )
    )

  bar_data_list[[rec]] <- site_means
}

# -----------------------------
# 4. Matrix builder
# -----------------------------
make_matrix <- function(dat, rows, sites) {
  z  <- matrix(NA_real_, length(rows), length(sites),
               dimnames = list(rows, sites))
  st <- matrix(NA_character_, length(rows), length(sites),
               dimnames = list(rows, sites))

  for (i in seq_len(nrow(dat))) {
    r <- match(dat$aa[i], rows)
    c <- match(dat$site[i], sites)
    z[r, c]  <- dat$score[i]
    st[r, c] <- dat$status[i]
  }

  list(score = z, status = st)
}

make_hover <- function(score, status) {

  out <- matrix("", nrow(score), ncol(score))

  out[status == "wt"]      <- "WT"
  out[status == "missing"] <- "Missing Datapoint"
  out[status == "score"]   <- sprintf("%.3f", score[status == "score"])

  out
}

# -----------------------------
# 5. Split by receptor
# -----------------------------
receptors <- unique(df2$receptor)

score_mats  <- list()
status_mats <- list()
hover_mats  <- list()

for (rec in receptors) {
  d <- filter(df2, receptor == rec)
  m <- make_matrix(d, aa_levels, as.character(sites))

  score_mats[[rec]]  <- m$score
  status_mats[[rec]] <- m$status
  hover_mats[[rec]]  <- make_hover(m$score, m$status)
}

score_vals <- df$score[!is.na(df$score)]
ext <- max(abs(range(score_vals)))

# -----------------------------
# 6. Plot
# -----------------------------
fig_heat <- plot_ly()
fig_bar  <- plot_ly()

for (i in seq_along(score_mats)) {

  rec <- names(score_mats)[i]
  score  <- score_mats[[rec]]
  status <- status_mats[[rec]]
  hover  <- hover_mats[[rec]]

  vis <- (i == 1)
  xvals <- as.numeric(colnames(score))

  # =========================
  # 1. Build strip FIRST ✅
  # =========================

  strip_z <- matrix(NA, nrow = length(aa_levels), ncol = 1)

  for (g in names(group_bounds)) {
    aas <- group_bounds[[g]]
    idx <- match(aas, aa_levels)
    strip_z[idx, 1] <- g
  }

  group_levels <- names(group_colors)

  strip_numeric <- matrix(
    match(strip_z, group_levels),
    nrow = length(aa_levels),
    ncol = 1
  )

  # =========================
  # 2. SCORE layer
  # =========================

  z <- score
  z[status != "score"] <- NA

  fig_heat <- fig_heat %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = z,
    type = "heatmap",
    zmin = -ext,
    zmax = ext,
    colorscale = list(
      list(0, "#b2182b"),
      list(0.5, "#ffffff"),
      list(1, "#2166ac")
    ),
    showscale = TRUE,
    colorbar = list(title = "Score", len = 0.7),
    text = hover,
    hoverinfo = "text",
    visible = vis
  )

  # =========================
  # 3. STRIP layer ✅
  # =========================

  fig_heat <- fig_heat %>% add_trace(
    x = min(xvals) - 2,
    y = aa_levels,
    z = strip_numeric,
    type = "heatmap",
    showscale = FALSE,
    colorscale = lapply(seq_along(group_levels), function(i) {
      list((i-1)/(length(group_levels)-1), group_colors[[group_levels[i]]])
    }),
    zmin = 1,
    zmax = length(group_levels),
    hoverinfo = "skip",
    visible = vis
  )

  # =========================
  # 4. WT layer
  # =========================

  wt_mask <- ifelse(status == "wt", 1, NA)

  fig_heat <- fig_heat %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = wt_mask,
    type = "heatmap",
    colorscale = list(list(0,"#bdbdbd"), list(1,"#bdbdbd")),
    text = hover,
    hoverinfo = "text",
    showscale = FALSE,
    visible = vis
  )

  # =========================
  # 5. MISSING layer (LAST ✅)
  # =========================

  miss_mask <- ifelse(status == "missing", 1, NA)

  fig_heat <- fig_heat %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = miss_mask,
    type = "heatmap",
    colorscale = list(list(0,"#f0e442"), list(1,"#f0e442")),
    text = hover,
    hoverinfo = "text",
    showscale = FALSE,
    visible = vis
  )

  # =========================
  # 6. BAR trace
  # =========================

  bars <- bar_data_list[[rec]]

  fig_bar <- fig_bar %>% add_trace(
    x = bars$site,
    y = bars$mean_score,
    type = "bar",
    marker = list(color = bars$color),
    visible = vis,
    showlegend = FALSE
  )

  print(paste("Added bar trace for:", rec))
}
# -----------------------------
# 7. X-axis ticks + labels
# -----------------------------
tick_sites <- seq(
  from = ceiling(min(sites)/20)*20,
  to   = floor(max(sites)/20)*20,
  by   = 20
)

tick_annotations <- lapply(tick_sites, function(x) {
  list(
    x = x,
    y = 0,
    yshift = -10,
    text = as.character(x),
    xref = "x",
    yref = "paper",
    showarrow = FALSE,
    xanchor = "center",
    yanchor = "top",
    font = list(size = 12)
  )
})

# -----------------------------
# 8. AA GROUP STRIP
# -----------------------------
strip_shapes <- lapply(names(group_bounds), function(g) {

  aas <- group_bounds[[g]]
  idx <- aa_index[aas]

  list(
    type = "rect",

    # ✅ attach to the HEATMAP subplot explicitly
    xref = "x2 domain",
    yref = "y2",

    # ✅ fixed strip at left edge of heatmap panel
    x0 = 0,
    x1 = 0.025,

    # ✅ categorical AA coordinates
    y0 = aa_levels[max(idx)],
    y1 = aa_levels[min(idx)],

    fillcolor = group_colors[[g]],
    line = list(width = 0),

    layer = "above"
  )
})
# -----------------------------
# 9. GROUP LABELS (pre-midpoint version)
# -----------------------------
group_labels <- list(
  list(x=-0.02, y="E", text="Negative", showarrow=FALSE, xref="paper", yref="y2", xanchor="right"),
  list(x=-0.02, y="R", text="Positive", showarrow=FALSE, xref="paper", yref="y2", xanchor="right"),
  list(x=-0.02, y="T", text="Polar", showarrow=FALSE, xref="paper", yref="y2", xanchor="right"),
  list(x=-0.02, y="I", text="Nonpolar", showarrow=FALSE, xref="paper", yref="y2", xanchor="right"),
  list(x=-0.02, y="Y", text="Aromatic", showarrow=FALSE, xref="paper", yref="y2", xanchor="right")
)

fig <- subplot(
  fig_bar,
  fig_heat,
  nrows = 2,
  heights = c(0.25, 0.75),
  shareX = TRUE
)


# -----------------------------
# 10. Dropdown
# -----------------------------
n <- length(score_mats)

buttons <- list()

for (i in seq_len(n)) {

  vis <- rep(FALSE, length(fig$x$data))

  # ---- BAR traces (first n)
  vis[i] <- TRUE

  # ---- HEATMAP traces (grouped after)
  heat_start <- n + (i - 1) * 4 + 1
  heat_end   <- n + i * 4

  vis[heat_start:heat_end] <- TRUE

  buttons[[i]] <- list(
    method = "restyle",
    args = list("visible", vis),
    label = names(score_mats)[i]
  )
}
# -----------------------------
# 11. FINAL layout
# -----------------------------

fig <- fig %>% layout(

  title = "Fc variant heatmap",

  updatemenus = list(list(
    type = "dropdown",
    buttons = buttons
  )),


  # ✅ annotations still work
  annotations = c(
    tick_annotations,
    group_labels
  ),

  # -----------------------------
  # TOP panel (bar chart)
  # -----------------------------
  yaxis = list(
    title = "Mean Score",
    zeroline = TRUE
  ),

  # -----------------------------
  # BOTTOM panel (heatmap)
  # -----------------------------
  yaxis2 = list(
    autorange = "reversed"
  ),

  # -----------------------------
  # X-axes (shared)
  # -----------------------------
  xaxis = list(
    showticklabels = FALSE
  ),

  xaxis2 = list(
    range = c(min(sites)-2, max(sites)),
    showticklabels = FALSE,
    tickmode = "array",
    tickvals = tick_sites,
    ticks = "outside",
    ticklen = 6
  ),

  # -----------------------------
  # margins
  # -----------------------------
  margin = list(
    l = 140,
    b = 100
  ),

  height = 500
)
dir.create("output", showWarnings = FALSE)

htmlwidgets::saveWidget(fig, "output/index.html", selfcontained = FALSE)
