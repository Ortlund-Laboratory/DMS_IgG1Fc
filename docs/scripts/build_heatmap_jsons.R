library(plotly)
library(dplyr)
library(tidyr)
library(NGLVieweR)
library(htmltools)
library(htmlwidgets)
library(jsonlite)

# -----------------------------
# 1. Load input files

files <- list(

  # ---- Fucos
  DF001_FcgR1_fucos = "data/FcgR1_fucos_znorm_top_minus_bottom_escape.csv",
  DF002_FcgR2a131H_fucos = "data/FcgR2a131H_fucos_znorm_top_minus_bottom_escape.csv",
  DF003_FcgR2a131R_fucos = "data/FcgR2a131R_fucos_znorm_top_minus_bottom_escape.csv",
  DF004_FcgR2b_fucos = "data/FcgR2b_fucos_znorm_top_minus_bottom_escape.csv",
  DF005_FcgR3a158F_fucos = "data/FcgR3a158F_fucos_znorm_top_minus_bottom_escape.csv",
  DF006_FcgR3a158V_fucos = "data/FcgR3a158V_fucos_znorm_top_minus_bottom_escape.csv",

  # ---- Afucos
  DF007_FcgR2a131H_afucos = "data/FcgR2a131H_afucos_znorm_top_minus_bottom_escape.csv",
  DF008_FcgR2a131R_afucos = "data/FcgR2a131R_afucos_znorm_top_minus_bottom_escape.csv",
  DF009_FcgR2b_afucos = "data/FcgR2b_afucos_znorm_top_minus_bottom_escape.csv",
  DF010_FcgR3a158F_afucos = "data/FcgR3a158F_afucos_znorm_top_minus_bottom_escape.csv",
  DF011_FcgR3a158V_afucos = "data/FcgR3a158V_afucos_znorm_top_minus_bottom_escape.csv",

  # ---- Fucos - Afucos
  DF012_FcgR2a131H__fucos_minus_afucos = "data/FcgR2a131H_fucosminusafucos_znorm_top_minus_bottom_escape.csv",
  DF013_FcgR2a131R_fucos_minus_afucos = "data/FcgR2a131R_fucosminusafucos_znorm_top_minus_bottom_escape.csv",
  DF014_FcgR2b_fucos_minus_afucos = "data/FcgR2b_fucosminusafucos_znorm_top_minus_bottom_escape.csv",
  DF015_FcgR3a158F_fucos_minus_afucos = "data/FcgR3a158F_fucosminusafucos_znorm_top_minus_bottom_escape.csv",
  DF016_FcgR3a158V_fucos_minus_afucos = "data/FcgR3a158V_fucosminusafucos_znorm_top_minus_bottom_escape.csv",

  # ---- Afucos pairwise
  DF017_Afuc_2b_minus_2aH = "data/Afucos_FcgR2b_minus_Afucos_FcgR2a131H_znorm_scores.csv",
  DF018_Afuc_2b_minus_2aR = "data/Afucos_FcgR2b_minus_Afucos_FcgR2a131R_znorm_scores.csv",
  DF019_Afuc_2b_minus_3aF = "data/Afucos_FcgR2b_minus_Afucos_FcgR3a158F_znorm_scores.csv",
  DF020_Afuc_2b_minus_3aV = "data/Afucos_FcgR2b_minus_Afucos_FcgR3a158V_znorm_scores.csv",
  DF021_Afuc_2aH_minus_2aR = "data/Afucos_FcgR2a131H_minus_Afucos_FcgR2a131R_znorm_scores.csv",
  DF022_Afuc_2aH_minus_3aF = "data/Afucos_FcgR2a131H_minus_Afucos_FcgR3a158F_znorm_scores.csv",
  DF023_Afuc_2aH_minus_3aV = "data/Afucos_FcgR2a131H_minus_Afucos_FcgR3a158V_znorm_scores.csv",
  DF024_Afuc_2aR_minus_3aF = "data/Afucos_FcgR2a131R_minus_Afucos_FcgR3a158F_znorm_scores.csv",
  DF025_Afuc_2aR_minus_3aV = "data/Afucos_FcgR2a131R_minus_Afucos_FcgR3a158V_znorm_scores.csv",
  DF026_Afuc_3aF_minus_3aV = "data/Afucos_FcgR3a158F_minus_Afucos_FcgR3a158V_znorm_scores.csv",

  # ---- Fucos pairwise
  DF027_Fuc_2b_minus_2aH = "data/Fucos_FcgR2b_minus_Fucos_FcgR2a131H_znorm_scores.csv",
  DF028_Fuc_2b_minus_2aR = "data/Fucos_FcgR2b_minus_Fucos_FcgR2a131R_znorm_scores.csv",
  DF029_Fuc_2b_minus_3aF = "data/Fucos_FcgR2b_minus_Fucos_FcgR3a158F_znorm_scores.csv",
  DF030_Fuc_2b_minus_3aV = "data/Fucos_FcgR2b_minus_Fucos_FcgR3a158V_znorm_scores.csv",
  DF031_Fuc_2aH_minus_2aR = "data/Fucos_FcgR2a131H_minus_Fucos_FcgR2a131R_znorm_scores.csv",
  DF032_Fuc_2aH_minus_3aF = "data/Fucos_FcgR2a131H_minus_Fucos_FcgR3a158F_znorm_scores.csv",
  DF033_Fuc_2aH_minus_3aV = "data/Fucos_FcgR2a131H_minus_Fucos_FcgR3a158V_znorm_scores.csv",
  DF034_Fuc_2aR_minus_3aF = "data/Fucos_FcgR2a131R_minus_Fucos_FcgR3a158F_znorm_scores.csv",
  DF035_Fuc_2aR_minus_3aV = "data/Fucos_FcgR2a131R_minus_Fucos_FcgR3a158V_znorm_scores.csv",
  DF036_Fuc_3aF_minus_3aV = "data/Fucos_FcgR3a158F_minus_Fucos_FcgR3a158V_znorm_scores.csv",

  # ---- Individual Titrations (all Fucos)
  DF037_FcgR1_fucos_titration = "data/formatted_ordered_FcgR1_mutation_to_Ka.csv",
  DF038_FcgR2a131H_fucos_titration = "data/formatted_ordered_FcgR2a131H_mutation_to_Ka.csv",
  DF039_FcgR2a131R_fucos_titration = "data/formatted_ordered_FcgR2a131R_mutation_to_Ka.csv",
  DF040_FcgR2b_fucos_titration = "data/formatted_ordered_FcgR2b_mutation_to_Ka.csv",
  DF041_FcgR3a158F_fucos_titration = "data/formatted_ordered_FcgR3a158F_mutation_to_Ka.csv",
  DF042_FcgR3a158V_fucos_titration = "data/formatted_ordered_FcgR3a158V_mutation_to_Ka.csv",

  # ---- Differential Titrations (all Fucos)
  DF043_Fuc_titration_2b_minus_2aH = "data/formatted_FcgR2b-FcgR2a-131H_Ka_fractions.csv",
  DF044_Fuc_titration_2b_minus_2aR = "data/formatted_FcgR2b-FcgR2a-131R_Ka_fractions.csv",
  DF045_Fuc_titration_2b_minus_3aF = "data/formatted_FcgR2b-FcgR3a-158F_Ka_fractions.csv",
  DF046_Fuc_titration_2b_minus_3aV = "data/formatted_FcgR2b-FcgR3a-158V_Ka_fractions.csv",
  DF047_Fuc_titration_2aH_minus_2aR = "data/formatted_FcgR2a-131H-FcgR2a-131R_Ka_fractions.csv",
  DF048_Fuc_titration_2aH_minus_3aF = "data/formatted_FcgR2a-131H-FcgR3a-158F_Ka_fractions.csv",
  DF049_Fuc_titration_2aH_minus_3aV = "data/formatted_FcgR2a-131H-FcgR3a-158V_Ka_fractions.csv",
  DF050_Fuc_titration_2aR_minus_3aF = "data/formatted_FcgR2a-131R-FcgR3a-158F_Ka_fractions.csv",
  DF051_Fuc_titration_2aR_minus_3aV = "data/formatted_FcgR2a-131R-FcgR3a-158V_Ka_fractions.csv",
  DF052_Fuc_titration_3aF_minus_3aV = "data/formatted_FcgR3a-158F-FcgR3a-158V_Ka_fractions.csv"
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
  "G","A","V","L","I","M","P","F","Y","W"
)

group_bounds <- list(
  Negative = c("D","E"),
  Positive = c("K","R"),
  Polar    = c("H","C","S","T","N","Q"),
  Nonpolar = c("G","A","V","L","I","M","P"),
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

make_hover <- function(score, status, sites, aa_levels, wt_vec, offset) {

  out <- matrix("", nrow(score), ncol(score))

  for (i in seq_len(nrow(score))) {
    for (j in seq_len(ncol(score))) {

      aa    <- aa_levels[i]
      site  <- as.numeric(colnames(score)[j])
      wt    <- wt_vec[site - offset]

      if (status[i, j] == "wt") {
        value <- "WT"
      } else if (status[i, j] == "missing") {
        value <- "Missing Datapoint"
      } else {
        value <- sprintf("%.3f", score[i, j])
      }

      out[i, j] <- paste0(
        "<b>Site:</b> ", site,
        "<br><b>Mutation:</b> ", wt, " → ", aa,
        "<br><b>Value:</b> ", value
      )
    }
  }

  out
}

score_mats <- list()
status_mats <- list()
hover_mats <- list()

for (rec in unique(df2$receptor)) {
  m <- make_matrix(filter(df2, receptor == rec), aa_levels, as.character(sites))
  score_mats[[rec]] <- m$score
  status_mats[[rec]] <- m$status
  hover_mats[[rec]] <- make_hover(
  m$score,
  m$status,
  sites,
  aa_levels,
  wt_vec,
  offset
)
}

ext <- max(abs(df$score), na.rm = TRUE)

## -----------------------------
## 7B. EXPORT EACH DATASET TO JSON  ✅ ADD HERE
## -----------------------------

dir.create("datasets_json", showWarnings = FALSE)

rec_names <- names(score_mats)

for (i in seq_along(rec_names)) {

  rec <- rec_names[i]

  score  <- score_mats[[rec]]
  status <- status_mats[[rec]]
  hover  <- hover_mats[[rec]]
  bars   <- bar_data_list[[rec]]

  xvals <- as.numeric(colnames(score))

  ## ---- BAR TRACE (matches your plotly call)
  trace_bar <- list(
    x = bars$site,
    y = bars$mean_score,
    type = "bar",
    marker = list(color = bars$color),
    visible = FALSE,
    showlegend = FALSE
  )

  ## ---- HEATMAP (score layer ONLY, same as your plotting logic)
  z <- score
  z[status != "score"] <- NA

  trace_heat <- list(
    x = xvals,
    y = rownames(score),
    z = z,
    type = "heatmap",
    zmin = -ext,
    zmax = ext,
    colorscale = list(
      list(0,"#b2182b"),
      list(0.5,"#ffffff"),
      list(1,"#2166ac")
    ),
    text = hover,
    hovertemplate = "%",
    hoverinfo = "text",
    visible = FALSE
  )

  ## ---- FINAL JSON STRUCTURE (matches your dataset_1.json)
  dataset_json <- list(
    dataset_index = i-1,
    name = rec,
    traces = list(trace_bar, trace_heat)
  )

  ## ---- WRITE FILE
  outfile <- file.path("datasets_json", paste0("dataset_", i-1, ".json"))

  write_json(
    dataset_json,
    outfile,
    auto_unbox = TRUE,
    pretty = FALSE,
    na = "null"
  )
}

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


bars <- bar_data_list[[rec]]   # ✅ MUST be here

  fig_bar <- fig_bar %>% add_trace(
    x = bars$site,
    y = bars$mean_score,
    type = "bar",
    marker = list(color = bars$color),
    visible = vis,
    showlegend = FALSE
  )


# ---- SCORE layer
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
    list(0,"#b2182b"),
    list(0.5,"#ffffff"),
    list(1,"#2166ac")
  ),
  
 text = hover,
  hovertemplate = "%{text}<extra></extra>",				   
  hoverinfo = "text",
  visible = vis
)

# ---- STRIP ✅ MOVE HERE
strip_z <- matrix(NA, length(aa_levels), 1)

for (g in names(group_bounds)) {
  idx <- match(group_bounds[[g]], aa_levels)
  strip_z[idx,1] <- g
}

group_levels <- names(group_colors)

strip_numeric <- matrix(
  match(strip_z, group_levels),
  nrow = length(aa_levels),
  ncol = 1
)

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

 text = hover,
  hovertemplate = "%{text}<extra></extra>",
  visible = vis
)

# ---- WT layer
wt_mask <- ifelse(status == "wt", 1, NA)

fig_heat <- fig_heat %>% add_trace(
  x = xvals,
  y = rownames(score),
  z = wt_mask,
  type = "heatmap",
  
 text = hover,
  hovertemplate = "%{text}<extra></extra>",

  colorscale = list(list(0,"#bdbdbd"), list(1,"#bdbdbd")),
  showscale = FALSE,
  visible = vis
)

# ---- MISSING layer ✅ MUST BE LAST
miss_mask <- ifelse(status == "missing", 1, NA)

fig_heat <- fig_heat %>% add_trace(
  x = xvals,
  y = rownames(score),
  z = miss_mask,
  type = "heatmap",
  
 text = hover,
  hovertemplate = "%{text}<extra></extra>",

  colorscale = list(list(0,"#f0e442"), list(1,"#f0e442")),
  showscale = FALSE,
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
# ✅ NATIVE NGL STRUCTURE
# -----------------------------

# Pass highlight data to JS
highlight_json <- jsonlite::toJSON(highlight_sites, auto_unbox = TRUE)

data_js <- tags$script(HTML(paste0("
  window.highlight_sites = ", highlight_json, ";
")))

# Structure container
struct_block <- tags$div(
  style = "border:2px solid black; width:600px; height:600px;",
  tags$div(id = "ngl_viewer", style = "width:100%; height:100%;")
)

# Load NGL.js
ngl_lib <- tags$script(src = "https://cdn.jsdelivr.net/npm/ngl@latest/dist/ngl.js")
ngl_script <- tags$script(src = "ngl_script_v2.js")

# -----------------------------
# 10. Dropdown (RESTORE THIS)
# -----------------------------

n <- length(score_mats)
traces_per_heat <- 4

buttons <- list()

for (i in seq_len(n)) {

  vis <- rep(FALSE, length(fig$x$data))

  # ---- BAR traces
  vis[i] <- TRUE

  # ---- HEATMAP traces
  start <- n + (i - 1) * traces_per_heat + 1
  end   <- n + i * traces_per_heat

  vis[start:end] <- TRUE

  buttons[[i]] <- list(
    method = "restyle",
    args = list("visible", vis),
    label = names(score_mats)[i]
  )
}

fig <- fig %>% layout(
  title = "Fc variant heatmap",

  width = 2400,
  height = 500,   # ✅ reduced height as requested

  updatemenus = list(list(
    type = "dropdown",
    buttons = buttons,
    x = 0,
    y = 1.1,
    xanchor = "left",
    yanchor = "top"
  )),

  margin = list(
    t = 100,
    l = 170,   # ✅ slightly reduced to match smaller height
    r = 80,
    b = 80
  ),

  annotations = c(tick_annotations, group_labels),

  yaxis = list(
    title = "Mean Score"
  ),

  yaxis2 = list(
    autorange = "reversed",
    type = "category",
    categoryorder = "array",
    categoryarray = aa_levels,
    tickmode = "array",
    tickvals = aa_levels,
    ticktext = aa_levels
  ),

  xaxis = list(
    showticklabels = FALSE
  ),

  xaxis2 = list(
    range = c(min(sites)-2, max(sites)),
    showticklabels = FALSE,
    tickmode = "array",
    tickvals = tick_sites
  )
)

page_spacing <- tags$style(HTML("
  #structure_section {
    margin-top: 60px;   /* ✅ pushes structure section down */
  }
"))

final_page <- htmltools::browsable(
  tagList(

    h2("Fc Variant Heatmap"),
    fig,

    page_spacing,   # ✅ NEW

    tags$div(
      id = "structure_section",

      h2("Structure Mapping"),
      struct_block
    ),

    ngl_lib,
    data_js,
    ngl_script
  )
)


dir.create("output", showWarnings = FALSE)

htmltools::save_html(
  final_page,
  file = "./index.html"
)
