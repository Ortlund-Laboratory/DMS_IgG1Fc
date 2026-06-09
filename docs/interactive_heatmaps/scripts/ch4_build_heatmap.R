library(plotly)
library(dplyr)

# -----------------------------
# 1. Load input files (CORRECTED)
# -----------------------------

files <- list(
  FcgR2b = "data/FcgR2b_fucosminusafucos_fucos_minus_afucos_scores.csv",
  FcgR3a158F = "data/FcgR3a158F_fucosminusafucos_fucos_minus_afucos_scores.csv"
)

df_list <- lapply(names(files), function(rec) {

  d <- read.csv(files[[rec]], stringsAsFactors = FALSE)

  # ✅ Rename to expected columns
  d <- d %>% dplyr::rename(
    site  = site,
    aa    = mutation,
    score = scaledscore
  )

  d$site <- as.numeric(d$site)  # ensure numeric

  d$receptor <- rec

  # keep only what we need
  d <- d %>% dplyr::select(site, aa, score, receptor)

  d
})

df <- dplyr::bind_rows(df_list)

# sanity check
print(head(df))

read_fasta_simple <- function(file) {
  lines <- readLines(file)
  seq <- paste(lines[!grepl("^>", lines)], collapse = "")
  return(seq)
}

wt_seq <- read_fasta_simple("data/Fc_prot.fasta")

# convert to AA vector
wt_vec <- base::strsplit(wt_seq, "")[[1]]

# convert to AA vector
wt_vec <- base::strsplit(wt_seq, "")[[1]]

# ensure site numbering valid
offset <- 215  # FASTA starts at site 216

idx <- df$site - offset

if (any(idx < 1 | idx > length(wt_vec))) {
  stop("Site indexing is outside FASTA range — check offset")
}

# -----------------------------
# 2. Build full grid (ALL AA x sites)
# -----------------------------

library(tidyr)

# define full AA set explicitly (better than unique)
aa_levels <- c("A","C","D","E","F","G","H","I","K","L",
               "M","N","P","Q","R","S","T","V","W","Y")

sites <- sort(unique(df$site))

# expand data safely per receptor
df2 <- df %>%
  tidyr::complete(
    receptor,
    site = sites,
    aa = aa_levels
  )

# -----------------------------
# 4. Assign status automatically
# -----------------------------
offset <- 215

df2 <- df2 %>%
  mutate(
    wt_aa = wt_vec[site - offset],
    status = case_when(
      aa == wt_aa ~ "wt",
      !is.na(score) ~ "score",
      TRUE ~ "missing"
    )
  )

if (any((df2$site - offset) < 1 | (df2$site - offset) > length(wt_vec))) {
  stop("Site indexing is outside FASTA range — check offset!")
}

# -----------------------------
# 5. Matrix builder
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

  out <- matrix("", nrow(score), ncol(score), dimnames = dimnames(score))

  out[status == "wt"]      <- "Wildtype"
  out[status == "missing"] <- "No score"
  out[status == "score"]   <- sprintf("%.3f", score[status == "score"])

  out
}

make_points <- function(status, type) {
  idx <- which(status == type, arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)

  data.frame(
    site = colnames(status)[idx[, 2]],
    aa   = rownames(status)[idx[, 1]]
  )
}

# -----------------------------
# 6. Split by receptor
# -----------------------------
receptors <- unique(df2$receptor)

score_mats  <- list()
status_mats <- list()
hover_mats  <- list()

for (rec in receptors) {

  d <- dplyr::filter(df2, receptor == rec)

  m <- make_matrix(d, aa_levels, as.character(sites))

  score_mats[[rec]]  <- m$score
  status_mats[[rec]] <- m$status
  hover_mats[[rec]]  <- make_hover(m$score, m$status)
}

# -----------------------------
# 7. Color scale limits
# -----------------------------
score_vals <- df$score[!is.na(df$score)]
ext <- max(abs(range(score_vals)))

# -----------------------------
# 8. Plot
# -----------------------------
fig <- plot_ly()

for (i in seq_along(score_mats)) {

  rec <- names(score_mats)[i]
  score  <- score_mats[[rec]]
  status <- status_mats[[rec]]
  hover  <- hover_mats[[rec]]

  z <- score
  z[status != "score"] <- NA

  wt   <- make_points(status, "wt")
  miss <- make_points(status, "missing")

  vis <- (i == 1)

  # heatmap
  fig <- fig %>% add_trace(
    x = as.numeric(colnames(z)), y = rownames(z), z = z,
    type = "heatmap",
    zmin = -ext, zmax = ext,
    colorscale = list(
      list(0, "#b2182b"),
      list(0.5, "#ffffff"),
      list(1, "#2166ac")
    ),
    hoverinfo = "skip",
    visible = vis
  )

wt_mask <- ifelse(status == "wt", 1, NA)

fig <- fig %>% add_trace(
  x = colnames(score),
  y = rownames(score),
  z = wt_mask,
  type = "heatmap",
  colorscale = list(
    list(0, "black"),
    list(1, "black")
  ),
  showscale = FALSE,
  hoverinfo = "skip",
  visible = vis
)


miss_mask <- ifelse(status == "missing", 1, NA)

fig <- fig %>% add_trace(
  x = colnames(score),
  y = rownames(score),
  z = miss_mask,
  type = "heatmap",
  colorscale = list(
    list(0, "#f0e442"),
    list(1, "#f0e442")
  ),
  showscale = FALSE,
  hoverinfo = "skip",
  visible = vis
)

  # Hover layer
  fig <- fig %>% add_trace(
    x = colnames(score), y = rownames(score),
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

tick_sites <- seq(
  from = ceiling(min(sites) / 20) * 20,
  to   = floor(max(sites) / 20) * 20,
  by   = 20
)

tick_annotations <- lapply(tick_sites, function(x) {
  list(
    x = x,
    y = -0.5,                 # slightly below the axis
    text = as.character(x),
    xref = "x",
    yref = "paper",
    showarrow = FALSE,
    xanchor = "center",
    yanchor = "top",
    font = list(size = 12)
  )
})

fig <- fig %>% layout(
  annotations = tick_annotations,
  xaxis = list(
    showticklabels = FALSE   # <-- hide broken default ticks
  )
)
# -----------------------------
# 9. Dropdown
# -----------------------------
n <- length(score_mats)

buttons <- lapply(seq_len(n), function(i) {
  vis <- rep(FALSE, n*4)
  vis[((i-1)*4+1):(i*4)] <- TRUE
  list(method="restyle", args=list("visible", vis), label=names(score_mats)[i])
})

fig <- fig %>% layout(
  title = "Fc variant heatmap",
  updatemenus = list(list(type="dropdown", buttons=buttons)),
  xaxis = list(showticklabels=FALSE),
  yaxis = list(autorange="reversed"),
  margin = list(l=120)
)

htmlwidgets::saveWidget(fig, "output/index.html", selfcontained = FALSE)
