library(plotly)
library(dplyr)
library(tidyr)
library(NGLVieweR)
library(htmltools)
library(htmlwidgets)

# -----------------------------
# 1. Load input files
# -----------------------------
files <- list(
  FcgR2b = "data/FcgR2b_fucosminusafucos_fucos_minus_afucos_scores.csv",
  FcgR3a158F = "data/FcgR3a158F_fucosminusafucos_fucos_minus_afucos_scores.csv"
)

df_list <- lapply(names(files), function(rec) {
  d <- read.csv(files[[rec]], stringsAsFactors = FALSE)
  d <- d %>% rename(site = site, aa = mutation, score = scaledscore)
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
# 3. Amino acids and groups
# -----------------------------
aa_levels <- c(
  "D","E","K","R","H","C","S","T","N","Q",
  "G","V","L","I","M","P","F","Y","W"
)

group_bounds <- list(
  Negative = c("D","E"),
  Positive = c("K","R"),
  Polar    = c("H","C","S","T","N","Q"),
  Nonpolar = c("G","V","L","I","M","P"),
  Aromatic = c("F","Y","W")
)

group_colors <- c(
  Negative = "#d73027",
  Positive = "#2166ac",
  Polar    = "#984ea3",
  Nonpolar = "#f0e442",
  Aromatic = "#1b9e77"
)

sites <- sort(unique(df$site))

# -----------------------------
# 4. COMPLETE MATRIX
# -----------------------------
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

# -----------------------------
# 5. BAR DATA + HIGHLIGHT EXTRACTION
# -----------------------------
bar_data_list <- list()
highlight_sites <- list()

for (rec in unique(df2$receptor)) {

  d <- df2 %>% filter(receptor == rec, status == "score")

  site_means <- d %>%
    group_by(site) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

  mu <- mean(site_means$mean_score)
  sdv <- sd(site_means$mean_score)

  site_means <- site_means %>%
    mutate(
      color = case_when(
        mean_score > (mu + sdv) ~ "#2166ac",
        mean_score < (mu - sdv) ~ "#b2182b",
        TRUE ~ "#bdbdbd"
      )
    )

  bar_data_list[[rec]] <- site_means

  highlight_sites[[rec]] <- list(
    blue = site_means$site[site_means$color == "#2166ac"],
    red  = site_means$site[site_means$color == "#b2182b"]
  )
}

# -----------------------------
# 6. MATRIX + HOVER
# -----------------------------
make_matrix <- function(dat, rows, sites) {
  z <- matrix(NA, length(rows), length(sites), dimnames = list(rows, sites))
  st <- z

  for (i in seq_len(nrow(dat))) {
    r <- match(dat$aa[i], rows)
    c <- match(dat$site[i], sites)
    z[r, c] <- dat$score[i]
    st[r, c] <- dat$status[i]
  }

  list(score = z, status = st)
}

make_hover <- function(score, status) {
  out <- matrix("", nrow(score), ncol(score))
  out[status == "wt"] <- "WT"
  out[status == "missing"] <- "Missing Datapoint"
  out[status == "score"] <- sprintf("%.3f", score[status == "score"])
  out
}

score_mats <- list()
status_mats <- list()
hover_mats <- list()

for (rec in unique(df2$receptor)) {
  m <- make_matrix(filter(df2, receptor == rec), aa_levels, as.character(sites))
  score_mats[[rec]] <- m$score
  status_mats[[rec]] <- m$status
  hover_mats[[rec]] <- make_hover(m$score, m$status)
}

ext <- max(abs(df$score), na.rm = TRUE)

# -----------------------------
# 7. PLOTTING
# -----------------------------
fig_heat <- plot_ly()
fig_bar  <- plot_ly()

for (i in seq_along(score_mats)) {

  rec <- names(score_mats)[i]
  score <- score_mats[[rec]]
  status <- status_mats[[rec]]
  hover <- hover_mats[[rec]]

  vis <- (i == 1)
  xvals <- as.numeric(colnames(score))

  z <- score
  z[status != "score"] <- NA

  fig_heat <- fig_heat %>% add_trace(
    x = xvals, y = rownames(score), z = z,
    type = "heatmap",
    zmin = -ext, zmax = ext,
    colorscale = list(
      list(0,"#b2182b"),
      list(0.5,"#ffffff"),
      list(1,"#2166ac")
    ),
    text = hover,
    hoverinfo = "text",
    visible = vis
  )

  bars <- bar_data_list[[rec]]

  fig_bar <- fig_bar %>% add_trace(
    x = bars$site,
    y = bars$mean_score,
    type = "bar",
    marker = list(color = bars$color),
    visible = vis
  )
}
library(plotly)
library(dplyr)
library(tidyr)
library(NGLVieweR)
library(htmltools)
library(htmlwidgets)

# -----------------------------
# 1–7: KEEP ALL YOUR EXISTING CODE ABOVE UNCHANGED
# -----------------------------
# (Your current script up to fig <- subplot(...) stays exactly as-is)

# -----------------------------
# X-axis ticks + labels
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
# GROUP LABELS
# -----------------------------
group_labels <- list(
  list(
    x = -0.02,
    y = "E",
    text = "Negative",
    showarrow = FALSE,
    xref = "paper",
    yref = "y2",
    xanchor = "right"
  ),
  list(
    x = -0.02,
    y = "R",
    text = "Positive",
    showarrow = FALSE,
    xref = "paper",
    yref = "y2",
    xanchor = "right"
  ),
  list(
    x = -0.02,
    y = "T",
    text = "Polar",
    showarrow = FALSE,
    xref = "paper",
    yref = "y2",
    xanchor = "right"
  ),
  list(
    x = -0.02,
    y = "I",
    text = "Nonpolar",
    showarrow = FALSE,
    xref = "paper",
    yref = "y2",
    xanchor = "right"
  ),
  list(
    x = -0.02,
    y = "Y",
    text = "Aromatic",
    showarrow = FALSE,
    xref = "paper",
    yref = "y2",
    xanchor = "right"
  )
)

fig <- subplot(
  fig_bar,
  fig_heat,
  nrows = 2,
  heights = c(0.25, 0.75),
  shareX = TRUE
)

# -----------------------------
# ✅ NEW: Build TWO structures (2b + 3a)
# -----------------------------

make_struct <- function(sites_struct) {

  blue_sel <- if (length(sites_struct$blue) > 0) {
    paste0("resi ", paste(sites_struct$blue, collapse = " or resi "))
  } else "none"

  red_sel <- if (length(sites_struct$red) > 0) {
    paste0("resi ", paste(sites_struct$red, collapse = " or resi "))
  } else "none"

  NGLVieweR("data/Fc_structure.pdb") %>%
    addRepresentation("surface", param = list(opacity = 0.6)) %>%
    addRepresentation(
      "surface",
      param = list(sele = blue_sel, colorValue = "blue")
    ) %>%
    addRepresentation(
      "surface",
      param = list(sele = red_sel, colorValue = "red")
    )
}

struct_2b  <- make_struct(highlight_sites[["FcgR2b"]])
struct_3a  <- make_struct(highlight_sites[["FcgR3a158F"]])

# -----------------------------
# ✅ Wrap structures in divs
# -----------------------------
struct_block <- tagList(
  tags$div(
    id="struct_1",
    style="width:100%; height:600px;",
    htmltools::browsable(struct_2b)
  ),
  tags$div(
    id="struct_2",
    style="width:100%; height:600px
# -----------------------------
# ✅ JavaScript to sync dropdown
# -----------------------------


# -----------------------------
# 10. Dropdown (KEEP YOUR FIXED VERSION)
# -----------------------------
n <- length(score_mats)

buttons <- list()

for (i in seq_len(n)) {

  vis <- rep(FALSE, length(fig$x$data))

  # bar traces
  vis[i] <- TRUE

  # heatmap traces
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
# 11. Layout
# -----------------------------
fig <- fig %>% layout(
  title = "Fc variant heatmap",
  updatemenus = list(list(
    type = "dropdown",
    buttons = buttons
  )),
  annotations = c(tick_annotations, group_labels),
  yaxis = list(title = "Mean Score"),
  yaxis2 = list(autorange = "reversed"),
  xaxis = list(showticklabels = FALSE),
  xaxis2 = list(
    range = c(min(sites)-2, max(sites)),
    showticklabels = FALSE,
    tickmode = "array",
    tickvals = tick_sites
  )
)

# -----------------------------
# ✅ FINAL PAGE WITH STRUCTURE
# -----------------------------

style_fix <- tags$style(HTML("
#struct_1, #struct_2 {
  min-height: 600px;
}
"))


final_page <- htmltools::browsable(
  tagList(
    style_fix,
    h2("Fc Variant Heatmap"),
    fig,
    h2("Structure Mapping"),
    struct_block
  )
)
dir.create("output", showWarnings = FALSE)

htmltools::save_html(
  final_page,
  file = "output/index.html"
)
