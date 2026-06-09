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

sites <- sort(unique(df$site))

df2 <- df %>%
  complete(receptor, site = sites, aa = aa_levels) %>%
  mutate(
    wt_aa = wt_vec[site - offset],
    status = case_when(
      aa == wt_aa ~ "wt",
      !is.na(score) ~ "score",
      TRUE ~ "missing"
    )
  )

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
  out[status == "wt"]      <- "Wildtype"
  out[status == "missing"] <- "No score"
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
fig <- plot_ly()

for (i in seq_along(score_mats)) {

  rec <- names(score_mats)[i]
  score  <- score_mats[[rec]]
  status <- status_mats[[rec]]
  hover  <- hover_mats[[rec]]

  vis <- (i == 1)
  xvals <- as.numeric(colnames(score))

  z <- score
  z[status != "score"] <- NA

  # score heatmap
  fig <- fig %>% add_trace(
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
    hoverinfo = "skip",
    visible = vis
  )

  # WT (grey)
  wt_mask <- ifelse(status == "wt", 1, NA)

  fig <- fig %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = wt_mask,
    type = "heatmap",
    colorscale = list(list(0,"#bdbdbd"), list(1,"#bdbdbd")),
    showscale = FALSE,
    hoverinfo = "skip",
    visible = vis
  )

  # Missing (yellow)
  miss_mask <- ifelse(status == "missing", 1, NA)

  fig <- fig %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = miss_mask,
    type = "heatmap",
    colorscale = list(list(0,"#f0e442"), list(1,"#f0e442")),
    showscale = FALSE,
    hoverinfo = "skip",
    visible = vis
  )

  # Hover
  fig <- fig %>% add_trace(
    x = xvals,
    y = rownames(score),
    z = matrix(0, nrow(score), ncol(score)),
    type = "heatmap",
    colorscale = list(list(0,"rgba(0,0,0,0)"), list(1,"rgba(0,0,0,0)")),
    showscale = FALSE,
    customdata = hover,
    hovertemplate =
      paste0(
        "AA: %{y}<br>",
        "Site: %{x}<br>",
        "Value: %{customdata}<br>",
        "Receptor: ", rec,
        "<extra></extra>"
      ),
    visible = vis
  )
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

strip_shapes <- lapply(names(group_bounds), function(g) {

  aas <- group_bounds[[g]]
  idx <- aa_index[aas]

  y_top    <- min(idx) - 1.5
  y_bottom <- max(idx) - 0.5

  list(
    type = "rect",
    xref = "paper",
    yref = "y",
    x0 = -0.02,
    x1 = -0.01,
    y0 = y_bottom,
    y1 = y_top,
    fillcolor = group_colors[[g]],
    line = list(width = 0)
  )
})

# -----------------------------
# 9. GROUP LABELS (pre-midpoint version)
# -----------------------------
group_labels <- list(
  list(x=-0.02, y="E", text="Negative", showarrow=FALSE, xref="paper", yref="y", xanchor="right"),
  list(x=-0.02, y="R", text="Positive", showarrow=FALSE, xref="paper", yref="y", xanchor="right"),
  list(x=-0.02, y="T", text="Polar", showarrow=FALSE, xref="paper", yref="y", xanchor="right"),
  list(x=-0.02, y="I", text="Nonpolar", showarrow=FALSE, xref="paper", yref="y", xanchor="right"),
  list(x=-0.02, y="Y", text="Aromatic", showarrow=FALSE, xref="paper", yref="y", xanchor="right")
)

# -----------------------------
# 10. Dropdown
# -----------------------------
n <- length(score_mats)

buttons <- lapply(seq_len(n), function(i) {
  vis <- rep(FALSE, n*4)
  vis[((i-1)*4+1):(i*4)] <- TRUE
  list(method="restyle", args=list("visible", vis), label=names(score_mats)[i])
})

# -----------------------------
# 11. FINAL layout
# -----------------------------
fig <- fig %>% layout(
  title = "Fc variant heatmap",
  updatemenus = list(list(type="dropdown", buttons=buttons)),

  shapes = strip_shapes,

  annotations = c(
    tick_annotations,
    group_labels
  ),

  height = 500,   # ✅ balanced height

  xaxis = list(
    showticklabels = FALSE,
    tickmode = "array",
    tickvals = tick_sites,
    ticks = "outside",
    ticklen = 6
  ),

  yaxis = list(
    autorange = "reversed"
  ),

  margin = list(
    l = 140,
    b = 100
  )
)
dir.create("output", showWarnings = FALSE)

htmlwidgets::saveWidget(fig, "output/index.html", selfcontained = FALSE)
