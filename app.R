library(shiny)
library(shinycssloaders)
library(shinyWidgets)
library(colourpicker)
library(shinyhelper)
library(shinyjs)
library(imager)
library(EBImage)
library(Matrix)
library(MASS)
library(pracma)
library(magick)
library(png)
library(dplyr)
library(tidyr)
library(vctrs)
library(sp)
library(grDevices)
library(crayon)
library(ggplot2)
library(scattermore)
library(Seurat)
library(RANN)
library(data.table)
library(jsonlite)
library(zip)
library(base64enc)


options(shiny.maxRequestSize = 2000*1024^2)


`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

save_prior_png <- function(prior, name, outdir, max_px = 1200L) {
  filename <- file.path(outdir, paste0(name, ".png"))
  if (inherits(prior, "Image")) {
    img <- prior
  } else if (is.matrix(prior)) {
    img <- Image(prior, colormode = "Grayscale")
  } else if (is.array(prior) && length(dim(prior)) %in% c(2,3)) {
    img <- Image(prior)
  } else {
    stop("Unknown prior type: ", class(prior))
  }
  
  h <- dim(img)[1]; w <- dim(img)[2]
  if (max(h, w) > max_px) {
    sc <- max_px / max(h, w)
    img <- EBImage::resize(img, w = round(w*sc), h = round(h*sc))
  }
  writeImage(img, filename)
  return(filename)
}

apply_filters <- function(img, v1, v2) {
  thr <- v1
  img[img < thr] <- 0
  
  if (v2 == 1) {
    img <- (img - .5) * 2 + .5
    img <- pmin(pmax(img, 0), 1)
  }
  
  img
}

.EBIDW <- function(data,
                   gene,
                   width = dim(data$image)[2],
                   height = dim(data$image)[1],
                   sigma = 5, 
                   iterations = 3,
                   threshold = .65,
                   shp = 20,
                   normalize = TRUE,
                   verbose = FALSE,
                   passive = FALSE,
                   datatype = "noseu",
                   spot_brush = NULL) {
  
  
  coords <- as.data.frame(data$coords)
  lv2 <- NULL  
  if(gene %in% data$rn){
    v <- as.numeric(data$expression[rownames(data$expression) %in% gene, ])
  }else{
    if(!is.null(data$meta) && is.data.frame(data$meta) &&
       gene %in% colnames(data$meta)){
      
      lv <- data$meta[[gene]]
      if (is.null(lv)) lv <- data$meta[, gene, drop = TRUE]
      
      if (!is.null(dim(lv))) lv <- as.vector(lv[, 1, drop = TRUE])
      
      if (is.numeric(lv)) {
        
        v <- as.numeric(lv)
        
        v[!is.finite(v)] <- 0
        
      } else {
        v   <- as.integer(factor(lv))
        lv2 <- unique(lv)
        v2  <- unique(v)
        lv2 <- lv2[order(v2)]
      }
      
    }else{
      v <- NULL
    }
  }
  
  if (length(v) == 0) stop("Gene not found in object.")
  names(v) <- colnames(data$expression)
  
  if(nrow(coords) != length(v)){
    add0 <- rep(0, (nrow(coords) - length(v)))
    names(add0) <- setdiff(coords$barcode, names(v))
    v <- c(v, add0)
    v <- v[match(coords$barcode, names(v))]
  }
  
  xi <- pmin(pmax(round(coords[, ncol(coords)-1]), 1), width)
  yi <- pmin(pmax(round(coords[, ncol(coords)]), 1), height)
  
  mat <- sparseMatrix(
    i = yi,
    j = xi,
    x = v,
    dims = c(height, width)
  )
  mat <- as.matrix(mat)
  
  if(isFALSE(passive)){
    img <- as.cimg(mat)
    
    combined_sigma <- sqrt(iterations) * sigma
    img <- isoblur(img, combined_sigma)
  }else{
    img <- as.cimg(mat)
    if (!is.null(spot_brush) && is.finite(spot_brush) && spot_brush >= 3) {
      
      b <- as.integer(spot_brush)
      if (b %% 2 == 0) b <- b + 1L
      kern <- makeBrush(b, shape='disc')
    } else if(datatype == "noseu"){
      kern <- makeBrush(11, shape='disc')
    }else{
      kern <- makeBrush(5, shape='disc')
    }
    img <- EBImage::dilate(img,kern)
  }
  
  if (normalize) {
    img <- (img - min(img)) / (max(img) - min(img) + 1e-8)
    p <- EBImage::Image(as.matrix(img), colormode = "Grayscale")
    p <- EBImage::rotate(p,-90)
  }else{
    p <- EBImage::rotate(img,-90)
  }
  
  inv_rectifier <- function(x, thr = threshold, sharpness = shp) {
    1 / (1 + exp(sharpness * (thr - x)))
  }
  
  if(isFALSE(passive)){
    k <- p * inv_rectifier(p)
  }else{
    k <- p
  }
  
  
  if (is.null(lv2)) {
    attr(k, "secatt") <- "skip"
  } else {
    attr(k, "secatt") <- as.character(na.omit(lv2))
  }
  
  return(k)
}

col_gene <- function(img_norm, color = "red") {
  
  img <- NULL
  legend <- NULL
  
  if(!is.null(color)){
    col_rgb <- col2rgb(color) / 255
    R <- img_norm * col_rgb[1]
    G <- img_norm * col_rgb[2]
    B <- img_norm * col_rgb[3]
    
    img <- EBImage::rgbImage(red = R, green = G, blue = B)
    
    nsteps <- 5
    vals <- seq(0, 1, length.out = nsteps)
    cols <- rgb(col_rgb[1] * vals, col_rgb[2] * vals, col_rgb[3] * vals)
    legend <- data.frame(level = vals, color = cols)
    
  }else{
    
    levels <- sort(unique(as.vector(img_norm)))
    levels <- levels[levels != 0]
    
    palette <- rainbow(length(levels))
    
    R <- matrix(0, nrow = nrow(img_norm), ncol = ncol(img_norm))
    G <- R
    B <- R
    
    for (i in seq_along(levels)) {
      mask <- img_norm == levels[i]
      col_rgb <- col2rgb(palette[i]) / 255
      R[mask] <- col_rgb[1]
      G[mask] <- col_rgb[2]
      B[mask] <- col_rgb[3]
    }
    
    img <- EBImage::rgbImage(red = R, green = G, blue = B)
    legend <- data.frame(level = attributes(img_norm)$secatt, color = palette)
  }
  
  return(list(image = img, legend = legend))
  
}

multi_scale_autocorr <- function(img, micron_per_pixel = 1, shifts_um = c(0.1, 1, 10, 50, 100, 500)) {
  if(length(dim(img)) > 2){
    img <- Image(img,colormode = "grayscale")
    img <- img[,,1] + img[,,2] + img[,,3]
    img <- normalize(img)
  }else{
    img <- Image(img,colormode = "grayscale")
  }
  
  img <- as.matrix(img)
  
  nr <- nrow(img)
  nc <- ncol(img)
  
  shift_img <- function(mat, dx = 0, dy = 0) {
    
    out <- matrix(NA, nrow = nr, ncol = nc)
    xs <- max(1, 1 + dx):min(nc, nc + dx)
    ys <- max(1, 1 + dy):min(nr, nr + dy)
    xt <- max(1, 1 - dx):min(nc, nc - dx)
    yt <- max(1, 1 - dy):min(nr, nr - dy)
    out[yt, xt] <- mat[ys, xs]
    out
  }
  
  corr_valid <- function(a, b) {suppressWarnings(cor(as.vector(a), as.vector(b), use = "complete.obs"))}
  
  res <- numeric(length(shifts_um))
  names(res) <- paste0(shifts_um, "um")
  
  for (i in seq_along(shifts_um)) {
    shift_px <- round(shifts_um[i] / micron_per_pixel)
    if (shift_px < 1) shift_px <- 1
    
    img_shift_x <- shift_img(img, dx = shift_px, dy = 0)
    img_shift_y <- shift_img(img, dx = 0, dy = shift_px)
    
    corr_x <- corr_valid(img, img_shift_x)
    corr_y <- corr_valid(img, img_shift_y)
    
    res[i] <- mean(c(corr_x, corr_y), na.rm = TRUE)
    
    
  }
  
  return(data.frame(distance_um = shifts_um, autocorrelation = res))
}

plot_spatial_autocorrelation_profile <- function(df, threshold = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install ggplot2")
  
  ggplot2::ggplot(df, ggplot2::aes(x = distance_um, y = autocorrelation)) +
    ggplot2::geom_line(color = "blue", linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    ggplot2::scale_x_log10(breaks = df$distance_um) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Shift (µm)",
      y = "Autocorrelation",
      title = "Spatial Autocorrelation Decay"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

hist_match <- function(target, reference) {
  
  matched <- target
  for (ch in 1:3) {
    ref_vals <- reference[[ch]]
    tgt_vals <- as.vector(target[,,ch])
    n_samp   <- min(length(tgt_vals), 20000L)
    idx      <- sample.int(length(tgt_vals), n_samp)
    tgt_samp <- tgt_vals[idx]
    tgt_q <- quantile(tgt_samp, probs = seq(0, 1, length.out = 256))
    ref_q <- quantile(ref_vals, probs = seq(0, 1, length.out = 256))
    map <- approx(tgt_q, ref_q, xout = tgt_vals, rule = 2)$y
    matched[,,ch] <- matrix(map, nrow = dim(target)[1])
  }
  return(matched)
}

.emap <- function(img_rgb,
                  e_quant = 0.70, h_quant = 0.40,
                  sigma_grad = 1.5, sigma_tens = 3,
                  gabor_scales = c(2,4,8), gabor_orients = seq(0, pi, length.out=4),
                  w_coh = 0.5, w_gabor = 0.3, w_density = 0.2) {
  
  img_rgb <- Image(img_rgb, colormode = "Color")
  
  colorMode(img_rgb) <- 2
  
  heM <- matrix(c(0.650,0.072,0.268,
                  0.704,0.990,0.570,
                  0.286,0.105,0.776), ncol=3, byrow=TRUE)
  heM <- heM / sqrt(rowSums(heM^2))
  od <- -log((imageData(img_rgb) + 1) / 256)
  stains <- solve(heM, matrix(od, nrow=3, byrow=TRUE))
  H <- normalize(matrix(stains[1, ], nrow=dim(img_rgb)[1]))
  E <- normalize(matrix(stains[2, ], nrow=dim(img_rgb)[1]))
  
  e_thr <- quantile(E, e_quant, na.rm=TRUE)
  h_thr <- quantile(H, h_quant, na.rm=TRUE)
  collagen_gate <- (E >= e_thr) & (H <= h_thr)
  
  kx <- matrix(c(-1,0,1,-2,0,2,-1,0,1), nrow=3)
  ky <- matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3)
  Es <- gblur(E, sigma=sigma_grad)
  gx <- filter2(Es, kx); gy <- filter2(Es, ky)
  Sxx <- gblur(gx*gx, sigma=sigma_tens)
  Syy <- gblur(gy*gy, sigma=sigma_tens)
  Sxy <- gblur(gx*gy, sigma=sigma_tens)
  trace <- Sxx + Syy
  disc  <- sqrt((Sxx - Syy)^2 + 4*Sxy^2)
  lam1 <- 0.5*(trace + disc); lam2 <- 0.5*(trace - disc)
  coherence <- normalize((lam1 - lam2) / (lam1 + lam2 + .Machine$double.eps))
  
  gabor_energy <- Image(0, dim = dim(E))
  for (s in gabor_scales) {
    for (th in gabor_orients) {
      k <- makeBrush(size = 7, shape="gaussian", sigma = s)
      gxth <-  cos(th)*gx + sin(th)*gy
      resp <- abs(gblur(gxth, sigma=s))
      gabor_energy <- pmax(gabor_energy, resp)
    }
  }
  gabor_energy <- normalize(gabor_energy)
  
  density <- normalize( E * (collagen_gate > 0) )
  
  stif_soft <- normalize(w_coh*coherence + w_gabor*gabor_energy + w_density*density)
  
  c_img <- Image(collagen_gate, colormode = "Grayscale")
  c_img <- normalize(c_img)
  c_img <- EBImage::rotate(c_img, 90)
  c_img <- EBImage::flop(c_img)
  
  list(H=H, E=E,
       coherence=coherence,
       gabor=gabor_energy,
       density=density,
       stif_soft=stif_soft,
       collagen_gate=collagen_gate,
       image=c_img)
}

.numap <- function(img_rgb, 
                   brush_size = 3,
                   log_offset = 0.001,
                   threshold_frac = 0.95,
                   min_long_edge = 1500L) {
  library(EBImage)
  
  
  img_input <- Image(img_rgb, colormode = "Color")
  in_dims   <- dim(img_input)
  long_edge <- max(in_dims[1:2])
  upsampled <- FALSE
  up_factor <- 1
  if (long_edge < min_long_edge) {
    up_factor <- min_long_edge / long_edge
    new_w <- as.integer(round(in_dims[1] * up_factor))
    new_h <- as.integer(round(in_dims[2] * up_factor))
    img <- EBImage::resize(img_input, w = new_w, h = new_h, filter = "bilinear")
    upsampled <- TRUE
    
  } else {
    img <- img_input
  }
  
  img_norm <- img
  for (i in 1:3) {
    frame <- img[,,i]
    if (diff(range(frame)) != 0) {
      img_norm[,,i] <- (frame - min(frame)) / (max(frame) - min(frame))
    }
  }
  
  I_od <- -log(img_norm + log_offset)
  w <- dim(img_norm)[1]
  h <- dim(img_norm)[2]
  I_od_mat <- matrix(as.vector(I_od), ncol = 3, byrow = FALSE)
  
  stain_matrix <- matrix(c(0.650, 0.072,
                           0.704, 0.990,
                           0.286, 0.105), nrow=3, byrow=TRUE)[,1:2]
  H_E <- I_od_mat %*% ginv(t(stain_matrix))
  
  H_channel <- H_E[,1]
  
  
  H_channel[!is.finite(H_channel)] <- 0
  H_min <- min(H_channel, na.rm = TRUE)
  H_max <- max(H_channel, na.rm = TRUE)
  if (H_max > H_min) {
    H_channel <- (H_channel - H_min) / (H_max - H_min)
  } else {
    H_channel <- rep(0, length(H_channel))
  }
  H_img <- Image(matrix(H_channel, nrow=w, ncol=h), colormode="Grayscale")
  
  H_data <- imageData(H_img)
  H_data[!is.finite(H_data)] <- 0
  H_data[H_data < 0] <- 0
  if (max(H_data) > 0) H_data <- H_data / max(H_data)
  
  
  H_img <- Image(H_data, colormode = "Grayscale")
  
  th <- otsu(H_img)
  if (upsampled) {
    
    
    h_vals <- as.vector(imageData(H_img))
    h_vals <- h_vals[is.finite(h_vals) & h_vals > 0.05]
    if (length(h_vals) < 100) {
      th_used <- th * 0.95
    } else {
      th_used <- as.numeric(quantile(h_vals, 0.85, na.rm = TRUE))
    }
    th_used <- min(th_used, 0.99)
  } else {
    th_used <- th * threshold_frac
  }
  mask <- H_img > th_used
  
  n_above <- sum(as.vector(imageData(mask)))
  total_px <- length(as.vector(imageData(mask)))
  pct_above <- 100 * n_above / total_px
  method <- if (upsampled) "85th-percentile" else sprintf("Otsu*%.2f", threshold_frac)
  
  
  
  if (upsampled) {
    
    mask <- EBImage::Image(imageData(mask) * 1.0, colormode = "Grayscale")
    brush <- makeBrush(3L, "disc")
    mask <- opening(mask, brush)
    mask <- fillHull(mask)
  } else {
    brush <- makeBrush(brush_size, "disc")
    mask <- opening(mask, brush)
    mask <- fillHull(mask)
  }
  
  nuclei_mask <- bwlabel(mask)
  n_components <- max(as.vector(imageData(nuclei_mask)))
  
  
  
  if (upsampled) {
    nuclei_mask <- EBImage::resize(nuclei_mask,
                                   w = in_dims[1], h = in_dims[2],
                                   filter = "none")
  }
  
  nuclei_mask <- rotate(nuclei_mask, 90)
  nuclei_mask <- flop(nuclei_mask)
  
  return(nuclei_mask)
}

.fmap <- function(img_rgb,
                  win_size = 3) {
  
  img <- Image(img_rgb, colormode = "Color")
  img_array <- imageData(img)
  img_gray <- 0.2989 * img_array[,,1] + 0.5870 * img_array[,,2] + 0.1140 * img_array[,,3]
  img_gray <- Image(img_gray, colormode = "Grayscale")
  
  win <- makeBrush(win_size, shape = "disc")
  
  local_mean <- filter2(img_gray, win / sum(win))
  
  local_var <- filter2(img_gray^2, win / sum(win)) - local_mean^2
  
  eps <- 1e-6
  local_entropy <- - (local_mean * log(local_mean + eps) + 
                        (1-local_mean) * log(1-local_mean + eps))
  
  eim <- EBImage::rotate(local_entropy, 90)
  eim <- EBImage::flop(eim)
  
}

.diffmap <- function(nuclei_mask,
                     ecm_mask) {
  
  nuclei_log <- imageData(nuclei_mask)
  ecm_log    <- imageData(ecm_mask)
  nuclei_log[nuclei_log>0] <- 1
  ecm_log[ecm_log>0] <- 1
  
  other_mask <- nuclei_log + ecm_log
  other_mask[other_mask>0] <- 1
  
  other_mask <- 1-other_mask
  other_mask <- Image(other_mask, colormode = "Grayscale")
  return(other_mask)
}

.hmap <- function(data) {
  df_spots <- data$coords
  if(nrow(df_spots) < 3) stop("Need at least 3 centroids for a hull")
  
  img_dim <- data$img_dims
  nrow_img <- img_dim[1]
  ncol_img <- img_dim[2]
  
  hull_idx <- chull(df_spots$x, df_spots$y)
  hull_poly <- df_spots[hull_idx, ]
  
  grid_x <- rep(1:ncol_img, each = nrow_img)
  grid_y <- rep(1:nrow_img, ncol_img)
  
  inside <- point.in.polygon(grid_x, grid_y,
                             hull_poly$x, hull_poly$y) > 0
  
  hull_mask_mat <- matrix(inside, nrow = nrow_img, ncol = ncol_img, byrow = FALSE)
  
  hull_mask <- Image(hull_mask_mat, colormode = "Grayscale")
  
  hull_mask <- rotate(hull_mask, 90)
  hull_mask <- flop(hull_mask)
  
  
  return(hull_mask)
}

.hmap2 <- function(entropy_mask){
  mask <- Image(entropy_mask, colormode = "Grayscale")
  mask <- mask > otsu(as.array(mask))  
  return(mask)
}

.wmap <- function(img_rgb, 
                  hull_mask, 
                  tol_default = 0.95, 
                  sat_tol_default = 0.05,
                  min_size = 50,
                  morph_radius = 1) {   
  
  img <- Image(img_rgb, colormode = "Color")
  if (colorMode(img) != 2) stop("Image must be color (RGB)")
  
  r <- img[,,1]
  g <- img[,,2]
  b <- img[,,3]
  
  hull_idx <- which(hull_mask > 0, arr.ind = TRUE)
  r_h <- r[hull_mask > 0]
  g_h <- g[hull_mask > 0]
  b_h <- b[hull_mask > 0]
  
  bright_vals <- pmin(r_h, g_h, b_h)
  sat_vals <- pmax(r_h, g_h, b_h) - bright_vals
  
  tol_use <- tol_default
  sat_tol_use <- sat_tol_default
  if (length(bright_vals) > 10) {
    km_b <- try(kmeans(bright_vals, centers = 2), silent = TRUE)
    km_s <- try(kmeans(sat_vals, centers = 2), silent = TRUE)
    if (!inherits(km_b, "try-error") && !inherits(km_s, "try-error")) {
      bright_center <- max(km_b$centers)
      sat_center <- min(km_s$centers)
      tol_use <- max(bright_center * 0.9, tol_default)
      sat_tol_use <- min(sat_center * 1.2, sat_tol_default)
    }
  }
  
  mask_vec <- (r_h > tol_use) & (g_h > tol_use) & (b_h > tol_use) &
    ((pmax(r_h, g_h, b_h) - pmin(r_h, g_h, b_h)) < sat_tol_use)
  
  mask <- matrix(0, nrow = dim(img)[1], ncol = dim(img)[2])
  inside <- hull_mask > 0
  mask[inside] <- mask_vec
  mask <- Image(mask, colormode = "Grayscale")
  
  struct <- makeBrush(2*morph_radius + 1, shape = "disc")
  mask_clean <- opening(mask, struct)
  mask_clean <- closing(mask_clean, struct)
  
  mask_label <- bwlabel(mask_clean)
  obj_sizes <- table(mask_label)
  for (lab in as.numeric(names(obj_sizes))) {
    if (lab == 0) next
    if (obj_sizes[as.character(lab)] < min_size) mask_label[mask_label==lab] <- 0
  }
  mask_final <- mask_label > 0
  
  if (morph_radius > 1) mask_final <- medianFilter(mask_final, size = morph_radius)
  
  mask_final <- rotate(mask_final, 90)
  mask_final <- flop(mask_final)
  
  return(mask_final)
}

.tmap <- function(img_list, 
                  quantile_thresh = 0.95) {
  
  lapply(img_list, function(d) {
    
    img <- if (!inherits(d, "Image")) Image(d, colormode = "Grayscale") else d
    
    
    th <- quantile(d, quantile_thresh, na.rm = TRUE)
    mask <- img > th
    th2 <- quantile(d, 1-quantile_thresh, na.rm = TRUE)
    mask2 <- img < th2
    mask <- mask + mask2
    
    
    bwlabel(mask)
  })
}

.lagg <- function(diff_img, 
                  min_px_per_label = 100) {
  
  labels <- bwlabel(diff_img)
  label_mat <- imageData(labels)
  
  tabs <- table(label_mat)
  small_labels <- as.integer(names(tabs)[tabs < min_px_per_label])
  
  if(length(small_labels) == 0) return(labels)  
  
  mask_small <- matrix(label_mat %in% small_labels, nrow=nrow(label_mat), ncol=ncol(label_mat))
  mask_large <- matrix(!mask_small & label_mat > 0, nrow=nrow(label_mat), ncol=ncol(label_mat))
  
  coords_large <- which(mask_large, arr.ind = TRUE)
  coords_small <- which(mask_small, arr.ind = TRUE)
  
  if(is.null(ncol(coords_large))) coords_large <- matrix(coords_large, ncol=2)
  if(is.null(ncol(coords_small))) coords_small <- matrix(coords_small, ncol=2)
  
  labels_large <- label_mat[mask_large]
  
  nn <- nn2(data = coords_large, query = coords_small, k = 1)
  nearest_idx <- nn$nn.idx[,1]
  
  label_mat[mask_small] <- labels_large[nearest_idx]
  
  bwlabel(Image(label_mat, colormode = "Grayscale"))
}

.sr <- function(data, 
                gene = NULL) {
  
  coords <- data$coords
  if("in_tissue"%in%colnames(coords)){
    coords <- coords[coords$in_tissue == 1, ]
  }
  
  if(is.null(gene)){
    cnd <- intersect(c("ACTB","GAPDH","RPLP0","TUBB","HPRT1","PGK1"), data$rn)
    if (length(cnd) < 1){ 
      cnd <- rownames(data$expression)
    }
    cnd <- cnd[1]
    v <- data$expression[cnd, ]
    if (is.matrix(v)) v <- v[1, ]
    v <- data.frame(barcode = names(v), g = as.numeric(v))
  }else{
    if(gene %in% data$rn){
      cnd <- gene
      v <- data$expression[cnd, ]
      if (is.matrix(v)) v <- v[1, ]
      v <- data.frame(barcode = names(v), g = as.numeric(v))
    }else{
      
      cnd <- data$meta[[gene]]
      if (is.null(cnd)) cnd <- data$meta[, gene, drop = TRUE]
      if (!is.null(dim(cnd))) cnd <- as.vector(cnd[, 1, drop = TRUE])
      ensure_numeric <- function(x) {
        if (is.numeric(x)) {
          x
        } else {
          suppressWarnings(as.numeric(as.factor(x)))
        }
      }
      v <- data.frame(barcode = rownames(data$meta), g = ensure_numeric(cnd))
    }
  }
  
  if("barcode"%in%colnames(coords)){
    v <- v[v$barcode %in% coords$barcode, ]
    v <- v[match(coords$barcode, v$barcode), ]
    stopifnot(all(coords$barcode == v$barcode))
  }
  
  df <- cbind(coords[,colnames(coords)%in%c("x","y")],v$g)
  names(df)[3] <- "expression"
  l <- list()
  l$data <- df
  return(l)
}

.gmap <- function(data,
                  SRdata,
                  radius = 8,
                  fill = 0,
                  hull = TRUE) {
  
  w <- data$img_dims[2]
  h <- data$img_dims[1]
  
  gene_mat <- matrix(fill, nrow = h, ncol = w)
  
  xi <- round(SRdata$data$x)
  yi <- round(SRdata$data$y)
  
  keep <- xi >= 1 & xi <= w & yi >= 1 & yi <= h
  xi <- xi[keep]
  yi <- yi[keep]
  vals <- SRdata$data$expression[keep]
  
  gene_mat[cbind(yi, xi)] <- vals
  
  mask <- gene_mat > 0
  
  se <- makeBrush(size = 2*radius+1, shape = "disc", step = FALSE)
  
  mask_dil <- EBImage::dilate(mask, se)
  
  expanded_mat <- EBImage::rotate(mask_dil * max(gene_mat, na.rm = TRUE), -90)
  
  if(isTRUE(hull)){
    
    expanded_mat[expanded_mat>0] <- 1
    
    library(sp)
    
    compute_hull_from_matrix <- function(bin_mat, return_mask = TRUE) {
      coords <- which(bin_mat == 1, arr.ind = TRUE) 
      if (nrow(coords) < 3) stop("Need at least 3 points for a hull")
      
      colnames(coords) <- c("y","x")
      
      hull_idx <- chull(coords[, "x"], coords[, "y"])
      hull_coords <- coords[c(hull_idx, hull_idx[1]), ] 
      
      mask_hull <- NULL
      if (return_mask) {
        h <- nrow(bin_mat)
        w <- ncol(bin_mat)
        mask_hull <- matrix(0, nrow = h, ncol = w)
        
        grid <- expand.grid(x = 1:w, y = 1:h)
        pip <- point.in.polygon(grid$x, grid$y,
                                hull_coords[, "x"], hull_coords[, "y"])
        mask_hull[cbind(grid$y, grid$x)] <- pip > 0
      }
      
      return(list(
        hull_coords = hull_coords,  
        hull_mask = mask_hull       
      ))
    }
    
    cm <- compute_hull_from_matrix(expanded_mat)
    cm <- cm$hull_mask
    
    return(cm)
  }else{
    return(expanded_mat)
  }
}

BuildGenesis <- function(base_dir, 
                         hires = TRUE, 
                         filter_tissue = FALSE) {
  
  if("filtered_feature_bc_matrix.h5"%in%list.files(base_dir)){
    
    expr_mat <- Read10X_h5(paste0(base_dir,"/filtered_feature_bc_matrix.h5"))
    if(class(expr_mat)=="list"){
      expr_mat <- expr_mat[["Gene Expression"]]
    }
    expr_mat <- expr_mat[rowSums(expr_mat)!=0,]
    
  }else{
    expr_dir <- file.path(base_dir, "filtered_feature_bc_matrix")
    
    expr_mat <- Read10X(expr_dir)
    if(class(expr_mat)=="list"){
      expr_mat <- expr_mat[["Gene Expression"]]
    }
    expr_mat <- expr_mat[rowSums(expr_mat)!=0,]
    
  }
  
  
  spatial_dir <- file.path(base_dir, "spatial")
  
  pos_file <- file.path(spatial_dir, "tissue_positions.csv")
  if (!file.exists(pos_file)) {
    pos_file <- file.path(spatial_dir, "tissue_positions_list.csv")
    if (!file.exists(pos_file)) {
      library(arrow)
      pos_file <- file.path(spatial_dir,"tissue_positions.parquet")
      coords <- read_parquet(pos_file, col_select = NULL, as_tibble = TRUE)
    }else{
      coords <- fread(pos_file, header = FALSE)
      colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col",
                            "pxl_row_in_fullres", "pxl_col_in_fullres")
    }
    
  } else {
    coords <- fread(pos_file)
  }
  
  if ("in_tissue" %in% colnames(coords)) {
    n_before <- nrow(coords)
    coords <- coords[coords$in_tissue == 1, , drop = FALSE]
    
  }
  
  
  scalefactors <- fromJSON(file.path(spatial_dir, "scalefactors_json.json"))
  scale <- if (hires) scalefactors$tissue_hires_scalef else scalefactors$tissue_lowres_scalef
  
  coords$x <- coords$pxl_col_in_fullres * scale
  coords$y <- coords$pxl_row_in_fullres * scale
  
  spot_diameter_fullres <- scalefactors$spot_diameter_fullres
  spot_diameter_scaled <- spot_diameter_fullres * scale
  microns_per_pixel <- 55 / spot_diameter_scaled  
  
  
  img_file <- file.path(spatial_dir, if (hires) "tissue_hires_image.png" else "tissue_lowres_image.png")
  image <- readPNG(img_file)
  
  cyt_file <- file.path(spatial_dir, "cytassist_image.tiff")
  
  to_rgb_array <- function(x) {
    A <- as.array(x)
    d <- dim(A)
    
    
    if (length(d) == 4) {
      if (d[3] == 1 && d[4] >= 3) {
        A <- A[,,1,1:d[4]]
      } else if (d[4] == 1 && d[3] >= 3) {
        A <- A[,,1:d[3],1]
      }
    }
    
    
    if (length(dim(A)) != 3) {
      stop("Cannot coerce image to RGB array")
    }
    
    
    if (dim(A)[3] == 4) {
      A <- A[,,1:3]
    }
    
    
    if (dim(A)[3] == 1) {
      A <- array(rep(A, 3), dim = c(dim(A)[1:2], 3))
    }
    
    stopifnot(dim(A)[3] == 3)
    A
  }
  
  
  
  render_cytassist_to_visium <- function(cyt_tiff, vis_png) {
    
    img <- image_read(cyt_tiff)          
    img <- image_convert(img, "RGB") 
    img <- image_data(img, channels = "rgb")  
    
    
    img <- aperm(img, c(2, 3, 1))
    img <- as.numeric(img) / 255
    I_c <- as.cimg(img)
    
    
    I_v <- vis_png
    
    
    if (spectrum(I_c) == 1) {
      I_c <- as.cimg(rep(I_c, 3), dim = c(dim(I_c)[1:2], 3, 1))
    }
    
    
    sx <- width(I_v)  / width(I_c)
    sy <- height(I_v) / height(I_c)
    scale <- mean(c(sx, sy))
    
    I_ds <- imresize(
      I_c,
      scale = scale,
      interpolation = 3
    )
    
    
    
    
    
    rng <- range(I_ds, finite = TRUE)
    if (diff(rng) > 0) {
      I_ds <- (I_ds - rng[1]) / diff(rng)
    } else {
      I_ds[] <- 0
    }
    
    
    
    A <- to_rgb_array(I_ds)
    
    A <- EBImage::Image(A,colormode = "Color")
    A
  }
  
  
  if (file.exists(cyt_file)) {
    
    cytimage <- render_cytassist_to_visium(
      cyt_tiff = cyt_file,
      vis_png  = image
    )
    
    image <- cytimage
    
    
  }
  
  
  img_dims <- dim(image)
  
  coords$y <- img_dims[1] - coords$y 
  
  
  
  if (filter_tissue) {
    coords <- coords[coords$in_tissue == 1, ]
    
  }else{
    
  }
  
  common_barcodes <- intersect(colnames(expr_mat), coords$barcode)
  if (filter_tissue) {
    expr_mat <- expr_mat[, common_barcodes]
    coords <- coords[match(common_barcodes, coords$barcode), ]
    
  }else{
    coords <- coords[match(common_barcodes, coords$barcode), ]
    
  }
  
  mt <- file.path(spatial_dir,"metadata.txt")
  if(file.exists(mt)){
    mt <- as.data.frame(fread(mt))
    names(mt)[1] <- "sampleIDs"
    mt <- mt[mt$sampleIDs%in%common_barcodes,]
    mt <- mt[match(common_barcodes, mt$sampleIDs), ]
    sl <- mt[vapply(mt, is.numeric, logical(1))]
    sl_mat <- t(as(as.matrix(sl), "dgCMatrix"))
    expr_mat <- rbind(expr_mat, sl_mat)
    
  }else{
    mt <- data.frame("no.metadata")
    
  }
  
  obj <- list(
    expression = expr_mat,
    rn = rownames(expr_mat)[order(rownames(expr_mat))],
    coords = as.data.frame(coords),
    image = image,
    img_dims = img_dims,
    scale = scale,
    spot_diameter = spot_diameter_scaled,
    microns_per_pixel = microns_per_pixel,
    hires = hires,
    meta = mt,
    tissue = coords$in_tissue
  )
  class(obj) <- "genesis"
  
  return(obj)
}

GenesisImageAnalysis <- function(data,
                                 clip_to_hull = TRUE) {
  
  
  
  in_img <- data$image
  if(dim(in_img)[3]>3){
    i1 <- in_img[, , 1:3]   
    i2 <- in_img[,,4]
    if(max(i2)>1){i2 <- i2/max(i2)}
    in_img <- i1*i2+(1-i2)
  }
  
  s1 <- .fmap(in_img)
  
  
  s2 <- .numap(in_img)
  
  
  s3 <- .emap(in_img)
  s4 <- .diffmap(s2,s3$image)
  
  
  s5 <- .hmap(data)
  hs <- .hmap2(s1)
  
  g1 <- .sr(data)
  g2 <- .gmap(data,g1,hull=T)
  
  
  
  s6 <- .wmap(in_img,
              s5)
  
  s7 <- Image(in_img, colormode = "Color")
  s7 <- EBImage::rotate(s7, 90)
  s7 <- EBImage::flop(s7)
  
  s1 <- s1 * !s6
  s2 <- s2 * !s6
  s3$image <- s3$image * !s6
  s3$gabor <- t(s3$gabor) * !s6
  s3$stif_soft <- t(s3$stif_soft) * !s6
  s4 <- s4 * !s6
  
  
  s4 <- s4 * hs
  s1 <- s1 * hs
  s2 <- s2 * hs
  s3$image <- s3$image * hs
  s3$gabor <- s3$gabor * hs
  s3$stif_soft <- s3$stif_soft * hs
  s4 <- s4 * hs
  s3$stif_soft <- Image(s3$stif_soft, colormode = "grayscale")
  
  if(isTRUE(clip_to_hull)){
    s1 <- s1 * g2
    s2 <- s2 * g2
    s3$image <- s3$image * g2
    s3$gabor <- s3$gabor * g2
    s3$stif_soft <- s3$stif_soft * g2
    s4 <- s4 * g2
    
    ll <- list(original_image=s7,
               entropy_mask=s1,
               nuclei_mask=s2,
               total_tissue=s3$stif_soft,
               fibrous_tissue_mask=s3$image,
               soft_tissue_mask=s4,
               hull_mask=g2,
               white_mask=s6,
               .all_tissue_parameters=s3)
    
  }else{
    ll <- list(original_image=s7,
               entropy_mask=s1,
               nuclei_mask=s2,
               total_tissue=s3$stif_soft,
               fibrous_tissue_mask=s3$image,
               soft_tissue_mask=s4,
               hull_mask=hs,
               white_mask=s6,
               .all_tissue_parameters=s3)
  }
  
  attributes(ll)$wf <- "map"
  attributes(ll)$ssz <- data$microns_per_pixel
  
  
  
  
  return(ll)
  
}

.process_HE_image <- function(png_path, brush_size = 7, log_offset = 0.001) {
  img <- readImage(png_path)
  
  img_norm <- img
  for (i in 1:3) {
    frame <- img[,,i]
    if (diff(range(frame)) != 0) {
      img_norm[,,i] <- (frame - min(frame)) / (max(frame) - min(frame))
    } else {
      img_norm[,,i] <- frame
    }
  }
  
  I_od <- -log(img_norm + log_offset)
  
  w <- dim(img_norm)[1]
  h <- dim(img_norm)[2]
  I_od_mat <- matrix(as.vector(I_od), ncol = 3, byrow = FALSE)
  
  stain_matrix <- matrix(c(
    0.650, 0.072,
    0.704, 0.990,
    0.286, 0.105
  ), nrow = 3, byrow = TRUE)[,1:2]
  
  H_E <- I_od_mat %*% ginv(t(stain_matrix)) 
  
  H_channel <- H_E[,1]
  H_channel <- H_channel - min(H_channel)
  H_channel <- H_channel / max(H_channel)
  H_img <- matrix(H_channel, nrow = w, ncol = h)
  H_img <- Image(H_img)
  
  th <- otsu(H_img)
  H_thr <- H_img > th 
  
  brush <- makeBrush(brush_size, "disc")
  H_thr <- opening(H_thr, brush)
  H_thr <- closing(H_thr, brush)
  H_thr <- fillHull(H_thr)
  
  nuclei_mask <- bwlabel(H_thr)
  
  H_thr <- rotate(H_thr, 90)
  H_thr <- flop(H_thr)
  
  return(H_thr)
}

.map_continuous_colors <- function(values, col_low = "white", col_high = "red", na_color = "gray") {
  n <- length(values)
  if (n == 0 || all(is.na(values))) return(rep(na_color, n))
  
  is_na <- is.na(values)
  v <- values[!is_na]
  rng <- range(v, na.rm = TRUE)
  scaled <- (v - rng[1]) / (rng[2] - rng[1])
  
  ramp <- grDevices::colorRampPalette(c(col_low, col_high))
  n_colors <- 256
  palette <- ramp(n_colors)
  indices <- pmin(floor(scaled * (n_colors - 1)) + 1, n_colors)
  
  out <- character(n)
  out[!is_na] <- palette[indices]
  out[is_na] <- na_color
  out
}

seurat_to_genesis <- function(data){
  res <- list()
  
  if(inherits(data, "Seurat")){
    if(inherits(data@images[[1]], "VisiumV2")){
      coord_df <- data@images[[1]]@boundaries$centroids@coords %>%
        as.data.frame()
    } else {
      if(inherits(data@images[[1]], "VisiumV1")){
        coord_df <- data@images[[1]]@coordinates %>%
          as.data.frame()
        coord_df <- coord_df[,c(4,5)]
        names(coord_df) <- c("x","y")
      } else {
        coord_df <- data@images[[1]]@coordinates %>%
          as.data.frame()
      }
    }
  }
  res$coords <- coord_df
  
  res$meta <- as.data.frame(data@meta.data)
  
  
  .extract_expr <- function(seu) {
    assay <- seu@assays$Spatial
    if (is.null(assay)) return(NULL)
    
    .reattach_dimnames <- function(m, assay, seu) {
      if (is.null(m) || length(m) == 0) return(m)
      if (is.null(rownames(m))) {
        rn <- tryCatch(rownames(assay), error = function(e) NULL)
        if (is.null(rn) || length(rn) != nrow(m)) {
          rn <- tryCatch(Features(seu), error = function(e) NULL)
        }
        if (!is.null(rn) && length(rn) == nrow(m)) rownames(m) <- rn
      }
      if (is.null(colnames(m))) {
        cn <- tryCatch(colnames(assay), error = function(e) NULL)
        if (is.null(cn) || length(cn) != ncol(m)) {
          cn <- tryCatch(colnames(seu), error = function(e) NULL)
        }
        if (!is.null(cn) && length(cn) == ncol(m)) colnames(m) <- cn
      }
      m
    }
    
    
    .materialise <- function(m) {
      if (is.null(m)) return(NULL)
      if (inherits(m, "IterableMatrix")) {
        m <- tryCatch(as(m, "dgCMatrix"), error = function(e) {
          tryCatch(as.matrix(m), error = function(e2) NULL)
        })
      }
      m
    }
    
    
    .nonempty <- function(m) {
      !is.null(m) && length(m) > 0 && nrow(m) > 0 && ncol(m) > 0
    }
    
    
    try_via_accessor <- function(layer_name) {
      
      m <- tryCatch(
        SeuratObject::LayerData(seu, assay = "Spatial", layer = layer_name),
        error = function(e) NULL
      )
      m <- .materialise(m)
      if (.nonempty(m)) return(m)
      
      m <- tryCatch(
        SeuratObject::GetAssayData(seu, assay = "Spatial", slot = layer_name),
        error = function(e) NULL
      )
      m <- .materialise(m)
      if (.nonempty(m)) return(m)
      
      NULL
    }
    
    
    try_direct_layer <- function(canonical_name) {
      if (is.null(assay@layers)) return(NULL)
      nms <- names(assay@layers)
      if (is.null(nms) || !length(nms)) return(NULL)
      
      if (canonical_name %in% nms) {
        m <- .materialise(assay@layers[[canonical_name]])
        if (.nonempty(m)) return(m)
      }
      
      pat <- paste0("^", canonical_name, "(\\.|$)")
      hits <- nms[grepl(pat, nms)]
      if (!length(hits)) return(NULL)
      
      
      if (length(hits) == 1) {
        return(.materialise(assay@layers[[hits[1]]]))
      }
      
      
      mats <- lapply(hits, function(h) .materialise(assay@layers[[h]]))
      mats <- mats[vapply(mats, .nonempty, logical(1))]
      if (!length(mats)) return(NULL)
      
      
      tryCatch({
        rn_sets <- lapply(mats, rownames)
        if (all(vapply(rn_sets, function(x) !is.null(x), logical(1)))) {
          common_rn <- Reduce(intersect, rn_sets)
          if (length(common_rn) > 0) {
            mats <- lapply(mats, function(m) m[common_rn, , drop = FALSE])
          }
        }
        do.call(cbind, mats)
      }, error = function(e) {
        
        do.call(cbind, mats)
      })
    }
    
    
    m <- try_via_accessor("data")
    if (.nonempty(m)) {
      m <- .reattach_dimnames(m, assay, seu)
      return(m)
    }
    
    m <- try_direct_layer("data")
    if (.nonempty(m)) {
      m <- .reattach_dimnames(m, assay, seu)
      return(m)
    }
    
    
    counts_mat <- try_via_accessor("counts")
    if (!.nonempty(counts_mat)) {
      counts_mat <- try_direct_layer("counts")
    }
    if (!.nonempty(counts_mat)) return(NULL)
    
    counts_mat <- .reattach_dimnames(counts_mat, assay, seu)
    
    col_sums <- Matrix::colSums(counts_mat)
    col_sums[col_sums == 0] <- 1
    rn_keep <- rownames(counts_mat); cn_keep <- colnames(counts_mat)
    m <- Matrix::t(Matrix::t(counts_mat) / col_sums) * 1e4
    m <- log1p(m)
    if (is.null(rownames(m)) && !is.null(rn_keep)) rownames(m) <- rn_keep
    if (is.null(colnames(m)) && !is.null(cn_keep)) colnames(m) <- cn_keep
    
    m
  }
  
  expr <- .extract_expr(data)
  if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) {
    stop("Could not find expression data in this Seurat object after trying the ",
         "Seurat accessors (LayerData, GetAssayData), direct layer names ",
         "('data', 'counts'), and split-layer patterns ('data.*', 'counts.*'). ",
         "If this is a Seurat v5 BPCells-backed assay, try as(seu[['Spatial']]$counts, 'dgCMatrix') ",
         "before passing the object to NaVis.")
  }
  res$expression <- expr
  
  rn <- rownames(expr)
  if (is.null(rn) || !length(rn)) {
    stop("seurat_to_genesis: the Seurat expression matrix has no gene rownames. ",
         "This usually means the Assay5 layer was accessed without its parent ",
         "assay's rownames — please re-save the Seurat object or run ",
         "NormalizeData() before uploading.")
  }
  res$rn <- rn[order(rn)]
  
  img <- data@images[[1]]@image
  res$image <- img
  
  res$img_dims <- dim(img)
  
  sf <- data@images[[1]]@scale.factors
  
  
  spot_full <- if (!is.null(sf$spot)) sf$spot else
    if (!is.null(sf$fiducial)) sf$fiducial else 1
  long_edge <- max(res$img_dims[1:2])
  .pick <- function(a, b) if (!is.null(a) && !is.na(a)) a else b
  lowres_sf <- .pick(sf$lowres, sf$lowres_scalef)
  hires_sf  <- .pick(sf$hires,  sf$hires_scalef)
  cand <- c(lowres = lowres_sf, hires = hires_sf)
  cand <- cand[!is.na(cand) & !is.null(cand) & cand > 0]
  if (length(cand) == 0) {
    res$scale <- if (!is.null(lowres_sf)) lowres_sf else 1
    res$hires <- FALSE
  } else if (length(cand) == 1) {
    res$scale <- as.numeric(cand[[1]])
    res$hires <- names(cand)[1] == "hires"
  } else {
    pred <- c(lowres = 20000 * cand[["lowres"]],
              hires  = 20000 * cand[["hires"]])
    best <- names(which.min(abs(log(long_edge / pred))))
    res$scale <- as.numeric(cand[[best]])
    res$hires <- (best == "hires")
    if (res$hires) {
      
    }
  }
  
  nx <- res$coords$x * res$scale
  ny <- res$coords$y * res$scale
  res$coords$x <- ny
  res$coords$y <- nx 
  res$coords$y <- res$img_dims[1] - res$coords$y 
  
  res$spot_diameter <- spot_full * res$scale
  
  res$microns_per_pixel <- 55/res$spot_diameter
  
  res$tissue <- 1
  
  class(res) <- "genesis"
  
  
  
  return(res)
}

hv <- readRDS("histvals.RDS")

.step <- function(n, col, title_text, body, extra = NULL) {
  div(style="display:flex;gap:10px;margin-bottom:14px;",
      div(style=sprintf("width:24px;height:24px;border-radius:50%%;background:%s;
                       color:#fff;font-size:11px;font-weight:700;display:flex;
                       align-items:center;justify-content:center;flex-shrink:0;margin-top:1px;", col),
          as.character(n)),
      div(style="flex:1;",
          tags$b(style="font-size:12px;color:#2a3050;", title_text),
          tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.6;", body),
          extra
      )
  )
}

ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    tags$meta(name="viewport", content="width=device-width, initial-scale=1"),
    tags$style(HTML("

* { box-sizing: border-box; }
body, .container-fluid, .shiny-app-content {
  margin: 0; padding: 0;
  background: #e8eaf0 !important;
  font-family: 'Segoe UI', system-ui, sans-serif;
  font-size: 13px;
  color: #1a1c24 !important;
  height: 100vh;
  overflow: hidden;
}

#navis-root {
  width: 100vw;
  height: 100vh;
  position: relative;
  overflow: hidden;
  background: #d4d8e2;
}

#scope-stage {
  position: absolute;
  inset: 0;
  background: #1a1c24;
  display: flex;
  align-items: center;
  justify-content: center;
  overflow: hidden;
  cursor: grab;
}
#zoom-wrap {
  width: 100%; height: 100%;
  display: flex; align-items: center; justify-content: center;
  overflow: hidden; position: relative;
}
#zoom-inner {
  transform-origin: center center;
  transform: scale(1) translate(0px,0px);
  user-select: none;
}
#zoom-inner img {
  display: block;
  max-width: 90vw; max-height: 90vh;
  border-radius: 2px;
  pointer-events: none;
  box-shadow: 0 8px 40px rgba(0,0,0,0.5);
}

.fp {
  position: absolute;
  background: #ffffff;
  border: 1px solid #c8ccd8;
  border-radius: 10px;
  box-shadow: 0 4px 24px rgba(0,0,0,0.18);
  overflow: hidden;
  min-width: 240px;
  max-height: calc(100vh - 80px);
  z-index: 200;
  display: flex;
  flex-direction: column;
}
.fp-header {
  background: #f0f2f7;
  border-bottom: 1px solid #d0d4e0;
  padding: 7px 12px 7px 10px;
  display: flex;
  align-items: center;
  gap: 8px;
  cursor: grab;
  user-select: none;
  border-radius: 10px 10px 0 0;
}
.fp-header:active { cursor: grabbing; }
.fp-title {
  font-size: 11px;
  font-weight: 700;
  letter-spacing: 0.5px;
  text-transform: uppercase;
  color: #3a4060;
  flex: 1;
}
.fp-close {
  background: none; border: none;
  color: #8890aa; font-size: 15px;
  cursor: pointer; padding: 0 2px;
  line-height: 1; border-radius: 3px;
}
.fp-close:hover { background: #e0e4ef; color: #333; }
.fp-min {
  background: none; border: none;
  color: #8890aa; font-size: 15px;
  cursor: pointer; padding: 0 4px;
  line-height: 1; border-radius: 3px;
  margin-left: auto;
  font-weight: 700;
}
.fp-min:hover { background: #e0e4ef; color: #333; }
.fp.minimized .fp-body { display: none; }
.fp.minimized { height: auto !important; resize: none; }
.fp-body {
  padding: 10px 12px;
  overflow-y: auto;
  flex: 1 1 auto;
  min-height: 0;
  max-height: calc(100vh - 140px);
}
.fp-body::-webkit-scrollbar { width: 4px; }
.fp-body::-webkit-scrollbar-thumb { background: #c8ccd8; border-radius: 4px; }

.fp-section {
  font-size: 9px;
  font-weight: 700;
  letter-spacing: 1.5px;
  text-transform: uppercase;
  color: #7880a0;
  margin: 12px 0 5px 0;
  padding-bottom: 3px;
  border-bottom: 1px solid #e8eaf0;
}
.fp-section:first-child { margin-top: 0; }

.fp label, .fp .control-label {
  color: #4a5070 !important;
  font-size: 12px !important;
  font-weight: 600;
  margin-bottom: 3px;
}
.fp .form-control, .fp select, .fp .selectize-input {
  background: #f4f6fb !important;
  color: #1a1c24 !important;
  border: 1px solid #c8ccd8 !important;
  border-radius: 5px !important;
  font-size: 12px !important;
}
.fp .selectize-dropdown {
  background: #fff !important;
  color: #1a1c24 !important;
  border: 1px solid #c8ccd8 !important;
}
.fp .selectize-dropdown .option:hover,
.fp .selectize-dropdown .option.active {
  background: #e8edf8 !important;
}
.fp .radio label, .fp .checkbox label {
  color: #3a4060 !important;
  font-size: 12px !important;
  font-weight: 400;
}

.btn-fp-primary {
  background: #2c5be0 !important;
  color: #fff !important;
  border: none !important;
  border-radius: 6px !important;
  font-size: 12px !important;
  font-weight: 600 !important;
  padding: 8px 12px !important;
  width: 100% !important;
  margin-top: 4px !important;
  letter-spacing: 0.3px;
}
.btn-fp-primary:hover { background: #1e4ac8 !important; }
.btn-fp-secondary {
  background: #f0f2f7 !important;
  color: #3a4060 !important;
  border: 1px solid #c8ccd8 !important;
  border-radius: 6px !important;
  font-size: 11px !important;
  padding: 5px 10px !important;
  width: 100% !important;
  margin-top: 3px !important;
}
.btn-fp-secondary:hover { background: #e0e4ef !important; border-color: #2c5be0 !important; }

.btn-fp-icon {
  background: transparent !important;
  color: #7880a0 !important;
  border: 1px solid #d0d4e0 !important;
  border-radius: 4px !important;
  font-size: 11px !important;
  padding: 2px 7px !important;
  min-width: 28px !important;
  width: auto !important;
  margin: 0 !important;
  line-height: 1.2 !important;
  cursor: pointer;
}
.btn-fp-icon:hover { background: #e8ecf5 !important; color: #2c5be0 !important; border-color: #2c5be0 !important; }

a.btn-fp-icon { display: inline-block; text-decoration: none !important; }

.insight-section-row {
  display: flex;
  align-items: center;
  gap: 5px;
  margin: 12px 0 2px;
  padding-bottom: 4px;
  border-bottom: 1px solid #e0e4ef;
}
.insight-section-row .insight-section-title {
  flex: 1;
  font-size: 10px;
  font-weight: 700;
  letter-spacing: 1.5px;
  color: #3a4060;
  text-transform: uppercase;
}
.btn-fp-green {
  background: #1e7a45 !important;
  color: #fff !important;
  border: none !important;
  border-radius: 6px !important;
  font-size: 12px !important;
  font-weight: 600 !important;
  padding: 8px 12px !important;
  width: 100% !important;
  margin-top: 4px !important;
}
.btn-fp-green:hover { background: #155c34 !important; }

.fp .irs--shiny .irs-bar { background: #2c5be0; border: none; }
.fp .irs--shiny .irs-bar-edge { background: #2c5be0; border: none; }
.fp .irs--shiny .irs-handle { background: #2c5be0; border-color: #2c5be0; }
.fp .irs--shiny .irs-from, .fp .irs--shiny .irs-to,
.fp .irs--shiny .irs-single { background: #2c5be0; }
.fp .irs--shiny .irs-line { background: #d0d4e0; border: none; }
.fp .irs--shiny .irs-grid-text { color: #7880a0; }

.channel-card {
  display: flex; align-items: center; gap: 8px;
  background: #f8f9fc;
  border: 1.5px solid #dde0ec;
  border-radius: 7px;
  padding: 6px 8px;
  cursor: pointer;
  margin-bottom: 5px;
  transition: border-color .12s, background .12s;
}
.channel-card:hover { background: #eef0f8; border-color: #aab0cc; }
.channel-card.active { border-color: #2c5be0; background: #e8edf8; }
.channel-thumb {
  width: 38px; height: 38px;
  border-radius: 4px; object-fit: cover;
  border: 1px solid #dde0ec; flex-shrink: 0;
}
.channel-label {
  font-size: 11px; color: #2a3050;
  overflow: hidden; text-overflow: ellipsis;
  white-space: nowrap; flex: 1;
}
.channel-dot { width: 9px; height: 9px; border-radius: 50%; flex-shrink: 0; }

#navis-toolbar {
  position: absolute;
  top: 12px; left: 50%;
  transform: translateX(-50%);
  z-index: 400;
  display: flex; align-items: center; gap: 6px;
  background: rgba(255,255,255,0.93);
  border: 1px solid #c8ccd8;
  border-radius: 30px;
  padding: 5px 14px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.12);
  white-space: nowrap;
}
.tb-btn {
  background: none; border: none;
  color: #3a4060; font-size: 12px;
  padding: 5px 10px; border-radius: 20px;
  cursor: pointer; display: flex; align-items: center;
  gap: 5px; transition: background .12s;
}
.tb-btn:hover { background: #e8edf8; }
.tb-btn.active { background: #2c5be0; color: #fff; }
.tb-sep { width: 1px; height: 20px; background: #d0d4e0; margin: 0 2px; }
.tb-logo {
  display: flex; align-items: center;
  margin-right: 10px; flex-shrink: 0;
}
.tb-logo img {
  height: 42px; width: auto;
  filter: drop-shadow(0 1px 2px rgba(0,0,0,.15));
}

#zoom-controls {
  position: absolute;
  bottom: 20px; right: 20px;
  z-index: 300;
  display: flex; flex-direction: column; gap: 4px;
}
.zoom-btn-hw {
  width: 32px; height: 32px;
  background: rgba(255,255,255,0.9);
  border: 1px solid #c8ccd8;
  color: #3a4060; font-size: 17px;
  line-height: 30px; text-align: center;
  border-radius: 6px; cursor: pointer;
  box-shadow: 0 1px 4px rgba(0,0,0,0.1);
}
.zoom-btn-hw:hover { background: #fff; border-color: #2c5be0; color: #2c5be0; }

#click-marker {
  position: absolute; width: 20px; height: 20px;
  border: 2px solid #ffaa00; border-radius: 50%;
  pointer-events: none;
  transform: translate(-50%, -50%);
  display: none; z-index: 500;
  box-shadow: 0 0 8px #ffaa00aa;
}
#click-marker::before, #click-marker::after {
  content: ''; position: absolute; background: #ffaa00;
}
#click-marker::before { width: 1px; height: 10px; left:50%; top:-7px; transform:translateX(-50%); }
#click-marker::after  { width: 10px; height:1px; top:50%; left:-7px; transform:translateY(-50%); }

.compare-marker {
  position: absolute;
  width: 26px; height: 26px;
  border-radius: 50%;
  pointer-events: none;
  transform: translate(-50%, -50%);
  display: flex; align-items: center; justify-content: center;
  font-family: -apple-system, system-ui, sans-serif;
  font-weight: 700; font-size: 13px; color: #fff;
  z-index: 500;
  border: 2px solid #fff;
  box-shadow: 0 0 0 2px rgba(0,0,0,0.4), 0 0 8px rgba(0,0,0,0.5);
}
#compare-marker-1 { background: #2c5be0; }
#compare-marker-2 { background: #ff8844; }

.prior-row {
  display: flex; align-items: center; gap: 8px; margin-bottom: 6px;
}
.prior-label { font-size: 11px; color: #4a5070; width: 110px; flex-shrink:0; }
.prior-bar-bg {
  flex: 1; height: 8px; background: #e8eaf0;
  border-radius: 4px; overflow: hidden;
}
.prior-bar-fill { height: 100%; border-radius: 4px; transition: width .3s ease; }
.prior-value { font-size: 10px; color: #7880a0; width: 34px; text-align:right; flex-shrink:0; }

.spot-card {
  background: #f4f6fb; border: 1px solid #dde0ec;
  border-radius: 6px; padding: 8px 10px; margin-bottom: 6px;
}
.spot-card-title { font-size: 9px; font-weight:700; letter-spacing:1px; text-transform:uppercase; color:#7880a0; margin-bottom:3px; }
.spot-card-value { font-size: 22px; font-weight:800; color:#2c5be0; line-height:1.1; }
.spot-card-sub   { font-size: 9px; color:#9098b8; margin-top:2px; }

#privacy-notice {
  background: #fffbea;
  border: 1px solid #e8d080;
  border-left: 3px solid #d4a010;
  border-radius: 5px;
  padding: 7px 10px;
  font-size: 11px; color: #5a4a10;
  margin-bottom: 8px;
}

.modal-content { background: #fff !important; color: #1a1c24 !important; }
.modal-header  { border-bottom: 1px solid #e0e4ef !important; }
.modal-footer  { border-top:    1px solid #e0e4ef !important; }
.modal-title   { color: #2c5be0 !important; font-weight:700; }

.shiny-spinner-output-container .shiny-spinner { border-top-color: #2c5be0 !important; }

.btn-dl {
  background: #f0f2f7 !important; color: #3a4060 !important;
  border: 1px solid #c8ccd8 !important; border-radius:5px !important;
  font-size: 11px !important; padding: 5px 8px !important;
}
.btn-dl:hover { border-color: #2c5be0 !important; color:#2c5be0 !important; }

      #analysis-spinner {
        display: none;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        padding: 28px 0;
        gap: 10px;
        color: #7880a0;
        font-size: 12px;
      }
      .nav-spin {
        width: 28px; height: 28px;
        border: 3px solid #e0e4ef;
        border-top-color: #2c5be0;
        border-radius: 50%;
        animation: navspin 0.75s linear infinite;
      }
      @keyframes navspin { to { transform: rotate(360deg); } }

      #stage-busy-badge {
        position: absolute;
        top: 10px; right: 52px;
        background: rgba(44,91,224,0.9);
        color: #fff;
        font-size: 11px;
        font-weight: 600;
        padding: 4px 10px;
        border-radius: 20px;
        z-index: 400;
        display: none;
        align-items: center;
        gap: 6px;
        pointer-events: none;
      }
      .stage-spin {
        width: 12px; height: 12px;
        border: 2px solid rgba(255,255,255,0.4);
        border-top-color: #fff;
        border-radius: 50%;
        animation: navspin 0.6s linear infinite;
        flex-shrink: 0;
      }

      .enrich-bar-row {
        display:flex; align-items:center; gap:8px; margin-bottom:5px;
      }
      .enrich-label { font-size:10px; color:#4a5070; width:110px; flex-shrink:0; }
      .enrich-bar-bg {
        flex:1; height:8px; background:#e8eaf0; border-radius:4px; overflow:hidden;
        position:relative;
      }
      .enrich-bar-fill { height:100%; border-radius:4px; transition:width .3s; }
      .enrich-ref { position:absolute; top:0; left:50%; width:1px; height:100%;
                    background:#999; opacity:.5; }
      .enrich-val { font-size:10px; color:#7880a0; width:40px; text-align:right; flex-shrink:0; }
      .decouple-score {
        font-size:32px; font-weight:800; color:#2c5be0;
        text-align:center; padding:8px 0 4px;
      }
      .decouple-label { font-size:11px; color:#7880a0; text-align:center; margin-bottom:4px; }

    "))
  ),
  
  div(id = "navis-root",
      
      div(id = "scope-stage",
          div(id = "zoom-wrap",
              div(id = "zoom-inner",
                  withSpinner(
                    imageOutput("main_image", width="100%", height="auto", inline=FALSE),
                    color="#2c5be0"
                  )
              )
          ),
          div(id = "click-marker"),
          div(id = "stage-busy-badge",
              div(class = "stage-spin"),
              "Rendering…"
          ),
          div(id = "zoom-controls",
              div(class="zoom-btn-hw", id="btn-zoom-in",    "+"),
              div(class="zoom-btn-hw", id="btn-zoom-reset", "⊙"),
              div(class="zoom-btn-hw", id="btn-zoom-out",   "−")
          )
      ),
      
      div(id = "navis-toolbar",
          div(class="tb-logo",
              tags$img(src="navis_logo.png", alt="NaVis",
                       onerror="this.style.display='none'")
          ),
          div(class="tb-sep"),
          tags$button(class="tb-btn", id="tb-data",     onclick="togglePanel('fp-data')",     "📂 Data"),
          tags$button(class="tb-btn", id="tb-acquire",  onclick="togglePanel('fp-acquire')",  "⬡ Acquire"),
          tags$button(class="tb-btn", id="tb-channels", onclick="togglePanel('fp-channels')", "📋 Channels"),
          tags$button(class="tb-btn", id="tb-stage",    onclick="togglePanel('fp-stage')",    "🎚 Image Controls"),
          tags$button(class="tb-btn", id="tb-overlay",  onclick="togglePanel('fp-overlay')",  "🔆 H&E Overlay"),
          tags$button(class="tb-btn", id="tb-merge",    onclick="togglePanel('fp-merge')",    "⊕ Merge & Interface"),
          tags$button(class="tb-btn", id="tb-meta",     onclick="togglePanel('fp-meta')",     "🏷 Metadata"),
          tags$button(class="tb-btn", id="tb-enrich",   onclick="togglePanel('fp-enrich')",   "🔭 Insights"),
          div(class="tb-sep"),
          tags$button(class="tb-btn", id="tb-analyse",  onclick="toggleAnalyse()",            "🔬 Analyze"),
          div(class="tb-sep"),
          
          downloadLink("dl_session_all", "💾 Save all",
                       class="tb-btn", style="text-decoration:none;"),
          tags$button(class="tb-btn", onclick="togglePanel('fp-help')",   "? Help"),
          tags$button(class="tb-btn", onclick="togglePanel('fp-privacy')", "🔒 Privacy"),
          tags$a(href="https://github.com/Izzilab", target="_blank",
                 class="tb-btn", style="text-decoration:none;", "GitHub")
      ),
      
      div(id="fp-data", class="fp",
          style="top:80px; left:20px; width:270px; display:none;",
          div(class="fp-header",
              span("📂", style="font-size:13px;"),
              span("Data Source", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-data')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-data')", "✕")
          ),
          div(class="fp-body",
              div(id="privacy-notice",
                  HTML("🔒 <b>No data is stored.</b> All processing is session-only.
                For unpublished data, <a href='https://github.com/Izzilab'
                target='_blank' style='color:#8a6010;'>deploy locally</a>.")
              ),
              div(class="fp-section", "Upload type"),
              radioButtons("upload_type", NULL,
                           choices = c("NaVis / Seurat (.rds)" = "proc", "Raw Visium ZIP" = "raw"),
                           inline = FALSE),
              conditionalPanel("input.upload_type == 'proc'",
                               fileInput("proc_file", "Upload .rds file",
                                         accept=c(".rds",".Rds",".RDS"),
                                         buttonLabel="Browse…")
              ),
              conditionalPanel("input.upload_type == 'raw'",
                               uiOutput("wizard_ui")
              ),
              tags$button(class="btn-fp-secondary", id="load_demo_btn",
                          onclick="Shiny.setInputValue('load_demo', Math.random(), {priority:'event'})",
                          "▶  Load Example Dataset")
          )
      ),
      
      div(id="fp-acquire", class="fp",
          style="top:80px; left:310px; width:300px; display:none;",
          div(class="fp-header",
              span("⬡", style="font-size:13px;"),
              span("Acquire Channel", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-acquire')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-acquire')", "✕")
          ),
          div(class="fp-body",
              
              div(class="fp-section", "Variable"),
              selectizeInput("selected_gene", NULL, choices=NULL, multiple=TRUE,
                             options=list(maxOptions=50,
                                          placeholder="Type one or more genes… (multi = signature)",
                                          plugins=list("remove_button"),
                                          onItemAdd = I("function(){ this.close(); }"))),
              
              div(class="fp-section", "Appearance"),
              fluidRow(
                column(6,
                       colourpicker::colourInput("gene_color", "Channel colour", value="#e63060")
                ),
                column(6,
                       div(style="padding-top:2px;",
                           tags$p("Auto params", style="font-size:10px;color:#7880a0;margin:0 0 2px 0;"),
                           materialSwitch("auto_mode", NULL, value=TRUE, status="primary")
                       )
                )
              ),
              
              tags$div(style="display:none;",
                       materialSwitch("mlbp", NULL, value=TRUE, status="primary")
              ),
              
              conditionalPanel(condition="!input.auto_mode",
                               div(style="background:#f4f6fb;border-radius:6px;padding:8px;margin-top:6px;",
                                   div(class="fp-section", "Manual parameters"),
                                   sliderInput("sigma",      "Smoothing (σ)",  min=0, max=20, value=10, step=1),
                                   sliderInput("iterations", "Iterations",     min=1, max=10, value=3,  step=1),
                                   sliderInput("threshold",  "Threshold",      min=0, max=1,  value=0.5,step=0.05),
                                   sliderInput("sharpness",  "Sharpness",      min=1, max=50, value=30, step=1),
                                   div(style="display:flex;align-items:center;gap:8px;",
                                       materialSwitch("localtex","Custom decay image",value=FALSE,status="danger"),
                                       selectInput("imgcust",NULL,choices=NULL,width="110px")
                                   )
                               )
              ),
              
              actionButton("gen_gene", "⬡  Acquire Channel Image",
                           class="btn-fp-primary"),
              
              tags$details(
                style="margin-top:10px;",
                tags$summary(
                  style="font-size:10px;font-weight:700;letter-spacing:1px;
                   color:#7880a0;text-transform:uppercase;cursor:pointer;
                   padding:4px 0;user-select:none;",
                  "🧹 Leakage correction"
                ),
                div(style="background:#f4f6fb;border-radius:6px;padding:8px;margin-top:4px;",
                    tags$p(style="font-size:10px;color:#7880a0;margin:0 0 6px 0;line-height:1.5;",
                           "Subtracts a fraction of each spot's neighbour mean before interpolation,
               correcting for Visium RNA diffusion leakage (SpotClean method).
               0 = off. Typical range: 0.10 – 0.20. SELECT A GENE FIRST!"),
                    sliderInput("decontam_alpha", "Leakage coefficient (α)",
                                min=0, max=0.3, value=0.10, step=0.01, width="100%"),
                    div(style="display:flex;gap:6px;align-items:center;margin-top:4px;",
                        actionButton("auto_decontam", "⚡ Auto-detect α",
                                     class="btn-fp-secondary",
                                     style="flex:1;margin:0;font-size:10px;padding:4px 6px;"),
                        uiOutput("auto_decontam_status_ui")
                    ),
                    tags$p(style="font-size:9px;color:#aab0c0;margin:6px 0 0;",
                           "Auto-detect sweeps α values and finds the elbow of the architecture-coupling curve.")
                )
              ),
              
              tags$details(
                style="margin-top:10px;",
                tags$summary(
                  style="font-size:10px;font-weight:700;letter-spacing:1px;
                   color:#7880a0;text-transform:uppercase;cursor:pointer;
                   padding:4px 0;user-select:none;",
                  "💧 Dropout sensitivity (signatures)"
                ),
                div(style="background:#f4f6fb;border-radius:6px;padding:8px;margin-top:4px;",
                    tags$p(style="font-size:10px;color:#7880a0;margin:0 0 6px 0;line-height:1.5;",
                           "Visium has substantial dropout: a zero count at a single spot
               whose neighbours are positive is more likely a sampling artefact
               than true biological absence. This slider sets the fraction of
               6 nearest-neighbour spots that must also be zero for a zero to
               count as a ", tags$em("true zero"),
                           ". Below this threshold, the zero is treated as dropout and
               imputed from the mean of non-zero neighbours before the
               signature is gated."),
                    sliderInput("dropout_tau", "Neighbour-zero threshold (τ)",
                                min=0, max=1, value=0.5, step=0.05, width="100%"),
                    tags$p(style="font-size:9px;color:#aab0c0;margin:6px 0 0;line-height:1.5;",
                           tags$b("τ = 0"), " → impute every zero (permissive; signature covers
               maximal area, risks false positives). ",
                           tags$b("τ = 1"), " → never impute (strict AND; robust to halos
               but loses real positives to dropout). ",
                           tags$b("τ = 0.5 (default)"), " → impute only isolated zeros whose
               neighbours clearly disagree. Applied only when 2+ genes are
               selected.")
                )
              )
          )
      ),
      
      div(id="fp-channels", class="fp",
          style="top:80px; right:20px; width:250px; display:none;",
          div(class="fp-header",
              span("📋", style="font-size:13px;"),
              span("Channel List", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-channels')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-channels')", "✕")
          ),
          div(class="fp-body", style="max-height:400px; overflow-y:auto;",
              uiOutput("channel_rail_ui") %>% withSpinner(color="#2c5be0", size=0.5)
          )
      ),
      
      div(id="fp-overlay", class="fp",
          style="top:80px; right:280px; width:340px; display:none;",
          div(class="fp-header",
              span("🔆", style="font-size:13px;"),
              span("H&E Overlay", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-overlay')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-overlay')", "✕")
          ),
          div(class="fp-body",
              div(class="fp-section", "Overlay channel on H&E"),
              tags$p(style="font-size:11px;color:#7880a0;margin:0 0 8px 0;",
                     "Move the sliders to preview the blend live. When satisfied, press ",
                     tags$b("Acquire"), " to save it as a new channel."),
              selectInput("he_overlay_gene", "Gene channel", choices=NULL),
              
              div(style="height:180px; overflow:hidden; margin-bottom:10px;
                   display:flex; align-items:center; justify-content:center;
                   background:#000; border-radius:5px; padding:3px;",
                  withSpinner(imageOutput("he_overlay_preview",
                                          width="auto", height="170px", inline=TRUE),
                              color="#c46010", size=0.5)
              ),
              sliderInput("he_opacity", "H&E opacity (0 = gene only, 1 = H&E only)",
                          min=0, max=1, value=0.5, step=0.05, width="100%"),
              sliderInput("gene_brightness", "Gene brightness",
                          min=0.5, max=4, value=1.5, step=0.1, width="100%"),
              actionButton("do_he_overlay", "🔆  Acquire this Overlay",
                           class="btn-fp-primary",
                           style="background:linear-gradient(135deg,#7a3a00,#c46010) !important;
                            border-color:#c46010 !important;")
          )
      ),
      
      div(id="fp-merge", class="fp",
          style="top:80px; left:630px; width:320px; display:none;",
          div(class="fp-header",
              span("⊕", style="font-size:13px;"),
              span("Merge & Interface", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-merge')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-merge')", "✕")
          ),
          div(class="fp-body",
              div(class="fp-section", "Channels"),
              selectInput("img1", "Channel A", choices=NULL),
              selectInput("img2", "Channel B", choices=NULL),
              
              div(class="fp-section", "Interface parameters"),
              tags$div(
                style="font-size:11px; color:#7880a0; margin-bottom:4px;",
                "Spatial scale at which channels A and B are considered ",
                tags$b("adjacent"), " without overlapping."
              ),
              sliderInput("interface_reach_um",
                          "Interface reach (\u00b5m)",
                          min = 5, max = 100, value = 15, step = 1),
              
              tags$div(
                style="font-size:11px; color:#7880a0; margin-bottom:4px;",
                "Defines what intensity counts as ",
                tags$b("present"), " in each channel."
              ),
              selectInput("interface_thr_mode",
                          "Threshold method",
                          choices = c("Adaptive (Otsu, clamped)" = "adaptive",
                                      "Fixed value"             = "fixed"),
                          selected = "adaptive"),
              conditionalPanel(
                condition = "input.interface_thr_mode == 'fixed'",
                sliderInput("interface_thr_fixed",
                            "Fixed threshold",
                            min = 0.02, max = 0.50, value = 0.10, step = 0.01)
              ),
              
              actionButton("do_blend", "\u2295  Compute composite + agreement + interface",
                           class="btn-fp-green"),
              
              div(class="fp-section", "Last computation"),
              uiOutput("interface_readout_ui")
          )
      ),
      
      div(id="fp-stage", class="fp",
          style="bottom:70px; left:20px; width:290px; display:none;",
          div(class="fp-header",
              span("🎚", style="font-size:13px;"),
              span("Image Controls", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-stage')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-stage')", "✕")
          ),
          div(class="fp-body",
              uiOutput("layer_controls_ui"),
              div(class="fp-section", "Tools"),
              actionButton("calcs",     "📈 Spatial Autocorrelation",  class="btn-fp-secondary"),
              actionButton("lgnd",      "🎨 Show Legend",              class="btn-fp-secondary"),
              actionButton("allvals",   "📋 Image Parameters",         class="btn-fp-secondary"),
              actionButton("crisp",     "✨ Apply crispness pass",       class="btn-fp-secondary"),
              actionButton("freeze_img", "📸 Save current view",         class="btn-fp-secondary"),
              div(class="fp-section", "Export"),
              fluidRow(
                column(4, downloadButton("dn1","Image",  class="btn-dl")),
                column(4, downloadButton("dn2","All",    class="btn-dl")),
                column(4, downloadButton("dn3","Object", class="btn-dl"))
              )
          )
      ),
      
      div(id="fp-analysis", class="fp",
          style="top:80px; right:20px; width:290px; display:none;",
          div(class="fp-header",
              span("🔬", style="font-size:13px;"),
              span("Spot Analysis", class="fp-title"),
              tags$button(class="fp-min",
                          onclick="togglePanelMinimize('fp-analysis')",
                          "−"),
              tags$button(class="fp-close",
                          onclick="document.getElementById('fp-analysis').style.display='none';
                             document.getElementById('tb-analyse').classList.remove('active');
                             var s=document.getElementById('scope-stage');
                             var m=document.getElementById('click-marker');
                             if(s) s.style.cursor='grab';
                             if(m) m.style.display='none';
                             analyseMode=false;",
                          "✕")
          ),
          div(class="fp-body",
              div(id="analysis-hint",
                  style="color:#5a6080; font-size:11px; padding:14px 12px;",
                  div(style="font-weight:700; color:#3a4060; margin-bottom:8px;
                     font-size:12px;",
                      "Click on tissue to analyse:"),
                  div(style="display:grid; grid-template-columns:auto 1fr;
                     gap:6px 10px; align-items:start; line-height:1.5;",
                      tags$span(style="color:#2c5be0; font-weight:600;", "Left-click"),
                      tags$span("Inspect a single spot"),
                      tags$span(style="color:#ff8844; font-weight:600;", "Right-click"),
                      tags$span("Compare two regions (right-click twice)"),
                      tags$span(style="color:#1e7a45; font-weight:600;", "Shift+drag"),
                      tags$span("Draw a region of interest (ROI)")
                  )
              ),
              
              div(id="analysis-spinner",
                  div(class="nav-spin"),
                  "Analysing…"
              ),
              div(id="analysis-content", style="display:none;",
                  div(class="fp-section", "📍 Location"),
                  uiOutput("spot_coords_ui"),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🔬 H&E Crop"),
                      actionButton("zoom_spot_crop", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadButton("dl_spot_crop", "💾", class="btn-fp-icon",
                                     style="padding:2px 6px; margin-left:4px;",
                                     title="Download PNG")
                  ),
                  withSpinner(imageOutput("spot_crop", width="100%", height="120px"),
                              color="#2c5be0", size=0.5),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🎨 H&E Colour Histogram"),
                      actionButton("zoom_spot_hist", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadButton("dl_spot_hist", "💾", class="btn-fp-icon",
                                     style="padding:2px 6px; margin-left:4px;",
                                     title="Download CSV")
                  ),
                  withSpinner(plotOutput("spot_histogram", height="110px"),
                              color="#2c5be0", size=0.5),
                  
                  div(class="fp-section", "🧬 Expression at Spot"),
                  uiOutput("spot_expression_ui"),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🧱 Tissue Context"),
                      downloadButton("dl_spot_context", "💾", class="btn-fp-icon",
                                     style="padding:2px 6px; margin-left:4px;",
                                     title="Download CSV")
                  ),
                  uiOutput("spot_priors_ui"),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "📊 Channel vs Tissue Correlations"),
                      actionButton("zoom_spot_corr", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadButton("dl_spot_corr", "💾", class="btn-fp-icon",
                                     style="padding:2px 6px; margin-left:4px;",
                                     title="Download CSV")
                  ),
                  withSpinner(plotOutput("spot_correlations", height="130px"),
                              color="#2c5be0", size=0.5),
                  
                  div(style="margin-top:14px; text-align:center;",
                      downloadButton("dl_spot_all", "💾 Save all as ZIP",
                                     class="btn-fp",
                                     style="width:auto;")
                  )
              ),
              
              div(id="compare-content", style="display:none;",
                  uiOutput("compare_drawer_ui")
              ),
              
              div(id="roi-content", style="display:none;",
                  uiOutput("roi_drawer_ui")
              )
          )
      ),
      
      div(id="fp-meta", class="fp",
          style="top:80px; left:340px; width:300px; display:none;",
          div(class="fp-header",
              span("🏷", style="font-size:13px;"),
              span("Metadata Channel", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-meta')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-meta')", "✕")
          ),
          div(class="fp-body",
              div(class="fp-section", "Generate image from metadata"),
              tags$p(style="font-size:11px;color:#7880a0;margin:0 0 8px 0;",
                     "Select a metadata column to render it as a NaVis channel.
                     Categorical columns are rainbow-coded; numeric columns can be
                     super-resolved via AD-EBIDW using the current Acquire parameters."),
              selectInput("meta_col", "Metadata column", choices=NULL, width="100%"),
              uiOutput("meta_col_typehint"),
              colourpicker::colourInput("meta_color", "Colour (numeric only)", value="#44bb88"),
              actionButton("gen_meta", "🏷  Generate Metadata Image",
                           class="btn-fp-primary",
                           style="background:linear-gradient(135deg,#3a006b,#7a20c0) !important;
                            border-color:#7a20c0 !important;"),
              
              div(class="fp-section", style="margin-top:14px;", "Super-resolve numeric column"),
              tags$p(style="font-size:11px;color:#7880a0;margin:0 0 6px 0;",
                     "Re-render the column with AD-EBIDW using the current Acquire parameters.
                     Active only for numeric columns."),
              uiOutput("meta_param_preview"),
              actionButton("super_meta", "🔬  Super-resolve column",
                           class="btn-fp-primary",
                           style="background:linear-gradient(135deg,#0a3c5e,#1976d2) !important;
                            border-color:#1976d2 !important; margin-top:6px;")
          )
      ),
      
      div(id="fp-enrich", class="fp",
          style="top:80px; left:640px; width:320px; display:none;",
          div(class="fp-header",
              span("🔭", style="font-size:13px;"),
              span("Spatial Insights", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-enrich')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-enrich')", "✕")
          ),
          div(class="fp-body",
              div(id="enrich-hint",
                  style="color:#7880a0;font-size:11px;text-align:center;padding:16px 0 10px;",
                  "Select a gene or signature channel, then press ",
                  tags$b("Run Analysis"),
                  " to compute spatial insights."
              ),
              
              uiOutput("enrich_target_ui"),
              actionButton("run_insights", "▶  Run Spatial Analysis",
                           class="btn-fp-primary",
                           style="background:linear-gradient(135deg,#2a4a7a,#2c5be0) !important;
                            border-color:#2c5be0 !important; margin-bottom:4px;"),
              actionButton("run_simgenes", "🔎  Find Similar Genes",
                           class="btn-fp-primary",
                           style="background:linear-gradient(135deg,#0f4d2e,#1e7a45) !important;
                            border-color:#1e7a45 !important; margin-bottom:6px;"),
              tags$p(style="font-size:10px;color:#7880a0;margin:2px 0 10px;text-align:center;",
                     "Spatial analysis runs on a 400-pixel downsample. Similar genes use spot-level expression."),
              
              div(id="enrich-content", style="display:none;",
                  
                  div(style="display:flex; gap:6px; margin-bottom:12px;",
                      downloadButton("dl_all_insights", "💾 Save all as ZIP",
                                     class="btn-fp-secondary",
                                     style="flex:1; margin:0; text-decoration:none;")
                  ),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🧭 Compartment Enrichment"),
                      actionButton("zoom_enrich", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadLink("dl_enrich",   "💾", class="btn-fp-icon", title="Save CSV")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Observed / expected signal per tissue compartment. ★ p<0.05, ★★ p<0.01 (permutation)."),
                  div(style="margin-bottom:14px;",
                      withSpinner(plotOutput("enrich_plot", height="140px"), color="#2c5be0", size=0.5)
                  ),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "✏️ Expression Boundaries"),
                      actionButton("zoom_boundary", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadLink("dl_boundary",   "💾", class="btn-fp-icon", title="Save PNG")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Edge gradient map. Bright = sharp expression boundary."),
                  div(style="height:180px; overflow:hidden; margin-bottom:10px;
                     display:flex; align-items:center; justify-content:center;
                     background:#000; border-radius:5px; padding:4px;",
                      withSpinner(imageOutput("boundary_img", width="auto", height="170px",
                                              inline=TRUE),
                                  color="#2c5be0", size=0.5)
                  ),
                  div(style="margin-bottom:14px;",
                      uiOutput("boundary_align_ui")
                  ),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "📡 Spatial Cross-correlation"),
                      actionButton("zoom_xcorr", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadLink("dl_xcorr",   "💾", class="btn-fp-icon", title="Save CSV")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Cross-correlation with each tissue prior at multiple lag distances (µm)."),
                  div(style="margin-bottom:14px;",
                      withSpinner(plotOutput("xcorr_plot", height="130px"), color="#2c5be0", size=0.5)
                  ),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🏝 Expression Islands"),
                      actionButton("zoom_island", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadLink("dl_island",   "💾", class="btn-fp-icon", title="Save CSV")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Connected expression domains above threshold: size, circularity, count."),
                  div(style="margin-bottom:6px;",
                      withSpinner(plotOutput("island_plot", height="110px"), color="#2c5be0", size=0.5)
                  ),
                  div(style="margin-bottom:14px;",
                      uiOutput("island_summary_ui")
                  ),
                  
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🔗 Histology Decoupling"),
                      actionButton("zoom_decouple", "🔍", class="btn-fp-icon", title="Zoom")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Fraction of spatial variance NOT explained by tissue architecture (1−R²).
             High = novel spatial pattern independent of H&E."),
                  div(style="margin-bottom:10px;",
                      uiOutput("decouple_ui")
                  )
              ),  
              
              div(id="simgenes-content", style="display:none;",
                  div(class="insight-section-row",
                      div(class="insight-section-title", "🔎 Similar Spatial Profiles"),
                      actionButton("zoom_simgenes", "🔍", class="btn-fp-icon", title="Zoom"),
                      downloadLink("dl_simgenes",   "💾", class="btn-fp-icon", title="Save CSV")
                  ),
                  tags$p(style="font-size:10px;color:#7880a0;margin:0 0 4px 0;",
                         "Top 10 genes with the most similar spatial expression profile (Pearson r, spot-level)."),
                  div(style="margin-bottom:10px;",
                      withSpinner(uiOutput("simgenes_ui"), color="#1e7a45", size=0.5)
                  )
              )
          )
      ),
      
      div(id="fp-help", class="fp",
          style="top:80px; left:50%; transform:translateX(-50%); width:520px; display:none;",
          div(class="fp-header",
              span("?", style="font-size:13px;font-weight:800;"),
              span("About & Tutorial", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-help')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-help')", "✕")
          ),
          div(class="fp-body",
              
              div(class="fp-section", "About NaVis"),
              HTML("
          <p style='color:#2a3050;font-size:12px;line-height:1.6;margin:0 0 8px;'>
          <b>NaVis</b> is an image-based framework for Visium spatial transcriptomics.
          Rather than treating expression as a table with coordinates as metadata, NaVis
          reconstructs each gene as a continuous spatial image and exposes five
          quantitative image-analytic operations unavailable in matrix-based tools:
          tissue-compartment enrichment, expression boundary detection, spatial
          cross-correlation at lag distances, expression island morphology, and
          histology decoupling scoring.
          </p>
          <div style='background:#fff3cd;border-left:3px solid #e8a020;border-radius:4px;
               padding:8px 10px;font-size:11px;color:#5a4010;margin-bottom:10px;'>
          <b>⚠ Spatial resolution.</b>
          Maps are <em>interpolated</em> via the AD-EBIDW algorithm.
          Visual sharpness is high, but precision is bounded by the technology
          (Visium: 55 µm spot diameter, 100 µm centre-to-centre pitch).
          These are spatially continuous inference maps, not single-cell measurements.
          </div>
        "),
              
              div(class="fp-section", "Tutorial — DCIS Breast Cancer Demo"),
              HTML("<p style='font-size:11px;color:#7880a0;margin:0 0 10px;'>
          Click <b>📂 Data → Load Demo</b> to start. Results described below
          were obtained on the bundled FFPE DCIS dataset.</p>"),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "1"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Load data and inspect the tissue"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Open ", tags$b("📂 Data"), " → ", tags$b("Load Demo"), ".
               A brief modal describes the dataset; loading runs in the background.
               Close the modal when you are ready. Five tissue-architecture channels appear in ",
                           tags$b("📋 Channels"), " — H&E Image, Nuclei, Fibrillar ECM, Soft Tissue,
               and Entropy. Click each to inspect them. Use mouse-wheel to zoom,
               click-drag to pan, ⊙ to reset.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "2"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Generate a single-gene map"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Open ", tags$b("⬡ Acquire"), ". Type ", tags$code("ERBB2"),
                           " and select it from the dropdown — it closes automatically on selection.
               Choose red as the colour. Leave ", tags$b("Auto"), " on.
               Expand ", tags$b("🧹 Leakage correction"), " — the slider is at α = 0.10 by
               default. Click ", tags$b("⚡ Auto-detect α"), " to let NaVis sweep α values
               and find the architecture-coupling elbow for this gene
               (typically 0.05–0.15 on FFPE Visium). Then click ",
                           tags$b("⬡ Acquire Channel Image"), ".
               The map appears in Channels and is selected automatically.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "3"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Tune the channel with Image Controls"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "The newly acquired ERBB2 channel still carries some halo from low-confidence
               spots — the AD-EBIDW interpolation extends each spot's signal a few hundred
               microns beyond its physical footprint, and the lowest-intensity tail of that
               extension is where most spurious-looking signal lives. Open ",
                           tags$b("🎚 Image Controls"), " to clean it up live."),
                    tags$ul(style="font-size:11px;color:#3a4060;margin:6px 0 0;
                           padding-left:16px;line-height:1.7;",
                            tags$li(tags$b("Noise gate"), " — pixels below this intensity are forced to
                       zero. Default 0.30; raise to 0.40–0.50 to drop diffuse background
                       and keep only the regions where the gene was confidently expressed.
                       The channel image redraws live as you drag the slider."),
                            tags$li(tags$b("Contrast boost"), " — applies an unsharp-mask-style boost
                       to the surviving signal. Default 0 (off); raise to 0.3–0.6 after
                       gating to recover crispness of the regions you kept.")
                    ),
                    tags$p(style="font-size:11px;color:#3a4060;margin:8px 0 0;line-height:1.55;",
                           "Three buttons under ", tags$b("Tools"), " expose secondary readouts:"),
                    tags$ul(style="font-size:11px;color:#3a4060;margin:4px 0 0;
                           padding-left:16px;line-height:1.7;",
                            tags$li(tags$b("📈 Spatial Autocorrelation"), " — computes Moran's-style
                       autocorrelation of the active channel at multiple lag distances and
                       plots the decay profile. The shape of the curve characterises the
                       channel's spatial scale: a sharp drop signals fine-grained patterning,
                       a slow decay signals broad gradients. Useful for comparing how
                       'clumpy' vs 'smooth' two channels are."),
                            tags$li(tags$b("🎨 Show Legend"), " — displays the colour-mapping legend for
                       the active channel (the gene name, the colour assigned, the dynamic
                       range). The legend is also embedded in any saved channel ZIP."),
                            tags$li(tags$b("📋 Image Parameters"), " — opens a modal listing every
                       parameter used to build the active channel: the gene or signature
                       name, leakage-correction α, noise gate, contrast boost, AD-EBIDW
                       σ and iteration count, dropout sensitivity τ (for signatures), and
                       any leakage-correction status. This is the channel's audit trail —
                       paste it into a notebook for reproducibility.")
                    ),
                    tags$p(style="font-size:11px;color:#3a4060;margin:8px 0 0;line-height:1.55;",
                           "Image Controls operate on the rendered image, not on the raw expression
               data — Insights and Find-Similar-Genes always re-read the underlying
               channel. Adjustments here are for visualisation; the analytical readouts
               are unaffected. Each channel keeps its own slider settings, so switching
               channels restores their independent values.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "4"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Run Spatial Insights on ERBB2"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Open ", tags$b("🔭 Insights"), " → click ",
                           tags$b("▶ Run Spatial Analysis"), ". Observed results on this dataset
               (ERBB2, auto α = 0.15, auto mode):"),
                    tags$ul(style="font-size:11px;color:#3a4060;margin:4px 0 0;
                           padding-left:16px;line-height:1.9;",
                            tags$li(tags$b("Compartment Enrichment"), " — near 1× across all three
                       compartments. An epithelial marker like ERBB2 occupies cellular
                       regions generally; the H&E-derived priors (nuclei density, fibrous
                       ECM, soft tissue) are architectural signals that do not separate
                       epithelial identity specifically, so 1× is the honest null.
                       Strong enrichment requires a signature specific to a compartment-
                       bound phenotype."),
                            tags$li(tags$b("Expression Boundaries"), " — gradient map shows sharp
                       ductal perimeters. Boundary alignment with tissue edges is moderate
                       (~0.2–0.3 Pearson r for all three priors), indicating ductal
                       boundaries align with general tissue-architecture transitions
                       rather than with any single compartment edge specifically."),
                            tags$li(tags$b("Spatial Cross-correlation"), " — peaks around r ≈ 0.5 at
                       lag 0 µm for Nuclei (ERBB2+ cells contain nuclei), decays
                       monotonically, crosses zero around 400 µm, and goes negative
                       beyond that. The 400 µm zero-crossing reflects this section's
                       typical duct-to-duct spacing; at longer distances ducts are
                       anti-correlated with each other because they are separated by
                       stroma."),
                            tags$li(tags$b("Expression Islands"), " — 102 islands above the 75th-percentile
                       threshold, median area 16 µm², mean area 5 783 µm², largest
                       island 586 000 µm², median circularity 0.79. The huge mean vs
                       median gap signals a heavily skewed size distribution: many
                       small intra-ductal hotspots plus a few very large connected
                       components that span entire ducts or groups of ducts merged
                       across the threshold."),
                            tags$li(tags$b("Histology Decoupling"), " — 63% (MODERATE). ERBB2's spatial
                       pattern is partly explained by the H&E priors but also carries
                       genuine transcriptional information beyond histology — expected
                       for a gene whose level varies between histologically similar
                       ducts.")
                    )
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "5"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Find spatially similar genes"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "With ERBB2 still selected, click ",
                           tags$b("🔎 Find Similar Genes"), ".
               Correlation runs over ", tags$em("all"), " spots (including zeros), so a
               gene expressed everywhere scores near zero — only genes with matching
               spatial ON/OFF patterns rank high."),
                    tags$p(style="font-size:11px;color:#3a4060;margin:6px 0 4px;line-height:1.55;",
                           "Observed top 10 for ERBB2 (Pearson r, spot-level):"),
                    tags$table(style="font-size:10px;color:#3a4060;border-collapse:collapse;
                              width:100%;margin:0 0 6px;",
                               tags$tr(style="background:#eef3ff;",
                                       tags$th(style="padding:3px 8px 3px 6px;text-align:left;width:28px;", "#"),
                                       tags$th(style="padding:3px 8px;text-align:left;",                    "Gene"),
                                       tags$th(style="padding:3px 8px;text-align:right;width:44px;",         "r"),
                                       tags$th(style="padding:3px 8px;text-align:left;",                    "Role")),
                               tags$tr(tags$td(style="padding:2px 8px 2px 6px;", "1"),
                                       tags$td(style="padding:2px 8px;", tags$code("AZGP1")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.813"),
                                       tags$td(style="padding:2px 8px;", "apocrine/luminal")),
                               tags$tr(style="background:#f8f9fc;",
                                       tags$td(style="padding:2px 8px 2px 6px;", "2"),
                                       tags$td(style="padding:2px 8px;", tags$code("PPDPF")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.774"),
                                       tags$td(style="padding:2px 8px;", "epithelial-enriched")),
                               tags$tr(tags$td(style="padding:2px 8px 2px 6px;", "3"),
                                       tags$td(style="padding:2px 8px;", tags$code("CD24")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.773"),
                                       tags$td(style="padding:2px 8px;", "luminal marker")),
                               tags$tr(style="background:#f8f9fc;",
                                       tags$td(style="padding:2px 8px 2px 6px;", "4"),
                                       tags$td(style="padding:2px 8px;", tags$code("MAL2")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.746"),
                                       tags$td(style="padding:2px 8px;", "apical / trafficking")),
                               tags$tr(tags$td(style="padding:2px 8px 2px 6px;", "5"),
                                       tags$td(style="padding:2px 8px;", tags$code("MUCL1")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.742"),
                                       tags$td(style="padding:2px 8px;", "mammary luminal")),
                               tags$tr(style="background:#f8f9fc;",
                                       tags$td(style="padding:2px 8px 2px 6px;", "6"),
                                       tags$td(style="padding:2px 8px;", tags$code("KRT7")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.739"),
                                       tags$td(style="padding:2px 8px;", "luminal keratin")),
                               tags$tr(tags$td(style="padding:2px 8px 2px 6px;", "7"),
                                       tags$td(style="padding:2px 8px;", tags$code("SPINT2")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.731"),
                                       tags$td(style="padding:2px 8px;", "epithelial protease inh.")),
                               tags$tr(style="background:#f8f9fc;",
                                       tags$td(style="padding:2px 8px 2px 6px;", "8"),
                                       tags$td(style="padding:2px 8px;", tags$code("ATG5")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.722"),
                                       tags$td(style="padding:2px 8px;", "autophagy (co-localised)")),
                               tags$tr(tags$td(style="padding:2px 8px 2px 6px;", "9"),
                                       tags$td(style="padding:2px 8px;", tags$code("KRT8")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.706"),
                                       tags$td(style="padding:2px 8px;", "luminal keratin")),
                               tags$tr(style="background:#f8f9fc;",
                                       tags$td(style="padding:2px 8px 2px 6px;", "10"),
                                       tags$td(style="padding:2px 8px;", tags$code("FOXA1")),
                                       tags$td(style="padding:2px 8px;text-align:right;", "0.694"),
                                       tags$td(style="padding:2px 8px;", "luminal TF"))
                    ),
                    tags$p(style="font-size:11px;color:#3a4060;margin:4px 0 0;line-height:1.55;",
                           "Nine of ten hits are canonical luminal/ductal epithelial markers
               (", tags$code("AZGP1"), ", ", tags$code("CD24"), ", ", tags$code("MUCL1"),
                           ", ", tags$code("KRT7"), ", ", tags$code("KRT8"), ", ", tags$code("FOXA1"),
                           ", etc.) — this is strong biological validation: the spatial correlation
               recovered the HER2+ ductal-epithelial compartment from scratch without
               any prior knowledge. Click any gene name in the results to load it
               into ", tags$b("⬡ Acquire"), " and visualise its spatial pattern.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#2c5be0;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "6"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Acquire a second gene and merge"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Acquire ", tags$code("ACTA2"),
                           " in green — a microenvironmental marker that highlights
               α-smooth-muscle-actin-expressing cells surrounding the ductal
               compartment (cancer-associated fibroblasts, myofibroblasts and
               myoepithelial cells all contribute). The resulting pattern is
               typically ductal-adjacent rather than a sharp perimeter ring —
               use it to visualise the tumour–microenvironment interface.
               Open ", tags$b("⊕ Merge & Interface"), ", set Channel A to ERBB2 and
               Channel B to ACTA2, optionally adjust the reach (μm) slider, then click
               ", tags$b("Compute"), ". Three new channels appear in ", tags$b("📋 Channels"),
                           ": the ", tags$b("composite"), " (blended A+B colors), the ",
                           tags$b("agreement map"), " (bright where A and B overlap), and the new ",
                           tags$b("interface map"), " (bright at the boundary between A's and B's
               territories, in each channel's color). The ", tags$b("Last computation"),
                           " section of the Merge panel reports interface area, A-vs-B asymmetry,
               and mean interface strength. Open ", tags$b("🎚 Image Controls"),
                           " to adjust per-layer noise gate and contrast independently.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#1e7a45;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "7"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Live H&E overlay"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Open ", tags$b("🔆 H&E Overlay"), ". Select ERBB2 from the dropdown.
               Move the opacity and brightness sliders — the preview updates live
               in the panel. When satisfied, click ",
                           tags$b("🔆 Acquire this Overlay"), " to commit it as a permanent channel
               that can be used in further merges and analysis.")
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:12px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#1e7a45;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "8"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Analyze the tissue: three click modalities"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "Click ", tags$b("🔬 Analyze"), " in the toolbar to enter Analyze mode.
               The Analyze panel hosts three distinct interactions, each tied to a
               different mouse action. Press ", tags$b("Escape"), " to exit at any time."),
                    
                    tags$div(style="margin-top:8px; padding:8px 10px; background:#f0f4ff;
                            border-left:3px solid #2c5be0; border-radius:4px;",
                             tags$b(style="font-size:11px;color:#2c5be0;",
                                    "Left-click — inspect a single point"),
                             tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                                    "Click any tissue location to open the spot drawer. It shows: the H&E crop
                 centred on that point, an RGB histogram of the local histology, a tissue
                 composition breakdown (Nuclei / Fibrillar ECM / Soft tissue summing to 100%),
                 the active channel's interpolated intensity at the clicked pixel, and Pearson
                 correlations of the channel image with each tissue prior.")
                    ),
                    
                    tags$div(style="margin-top:8px; padding:8px 10px; background:#fff5ec;
                            border-left:3px solid #ff8844; border-radius:4px;",
                             tags$b(style="font-size:11px;color:#ff8844;",
                                    "Right-click twice — compare two regions"),
                             tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                                    "Right-click anywhere to store region A (a blue numbered marker appears).
                 Right-click a second location to store region B (orange marker) — the
                 drawer flips to a side-by-side comparison: paired H&E crops, paired tissue
                 composition with log₂(B/A) fold-changes, and channel-intensity fold-changes.
                 Positive log₂ means the readout is higher in B than in A; negative means
                 the reverse. Changes within ±0.1 are shown as a gray dash (no meaningful
                 difference). Useful for quick contrasts: duct vs stroma, tumor vs adjacent
                 normal, lesion vs background. Switching channels auto-clears the comparison.")
                    ),
                    
                    tags$div(style="margin-top:8px; padding:8px 10px; background:#eef9f1;
                            border-left:3px solid #1e7a45; border-radius:4px;",
                             tags$b(style="font-size:11px;color:#1e7a45;",
                                    "Shift+drag — draw a region of interest (ROI)"),
                             tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                                    "Hold Shift and drag a rectangle on the tissue. The drawer opens with: an
                 H&E crop of the boxed region, tissue composition averaged inside the
                 ROI, and the top 20 genes by enrichment ratio (mean inside / mean outside).
                 Each gene row has a → button that pre-populates ", tags$b("⬡ Acquire"),
                                    " with that gene — click Acquire to load it as a channel.
                 Try drawing an ROI on a clear DCIS duct: top genes typically include
                 luminal markers (", tags$code("KRT8"), ", ", tags$code("CD24"),
                                    ", ", tags$code("EPCAM"), "). The ROI persists across channel switches —
                 use it to browse multiple genes through the same regional lens.
                 Click ", tags$b("✕ Clear ROI"), " to dismiss it.")
                    )
                  )
              ),
              
              div(style="display:flex;gap:10px;margin-bottom:16px;",
                  div(style="width:24px;height:24px;border-radius:50%;background:#7a3a00;color:#fff;
                     font-size:11px;font-weight:700;display:flex;align-items:center;
                     justify-content:center;flex-shrink:0;margin-top:1px;", "9"),
                  div(
                    tags$b(style="font-size:12px;color:#2a3050;", "Save results"),
                    tags$p(style="font-size:11px;color:#3a4060;margin:3px 0 0;line-height:1.55;",
                           "After running Insights, click ",
                           tags$b("💾 Save all as ZIP"), " to download all five analysis tables
               (compartment enrichment CSV, boundary alignment CSV, cross-correlation CSV,
               island morphology CSV, summary TXT with decoupling score) plus the boundary
               gradient PNG. Individual sections also have per-section 💾 download links.
               Zoom any panel with 🔍 for a larger, print-quality view of each result.")
                  )
              ),
              
              div(class="fp-section", "Navigation"),
              HTML("
          <table style='width:100%;font-size:11px;color:#3a4060;border-collapse:collapse;
                        margin-bottom:10px;'>
            <tr><td style='padding:3px 8px 3px 0;font-weight:600;width:36%;'>Mouse wheel</td>
                <td>Zoom in / out on the stage</td></tr>
            <tr style='background:#f4f6fb;'>
                <td style='padding:3px 8px 3px 0;font-weight:600;'>Click + drag</td>
                <td>Pan the tissue canvas</td></tr>
            <tr><td style='padding:3px 8px 3px 0;font-weight:600;'>⊙ button</td>
                <td>Reset zoom to fit</td></tr>
            <tr style='background:#f4f6fb;'>
                <td style='padding:3px 8px 3px 0;font-weight:600;'>Panel header</td>
                <td>Drag any floating panel to reposition it</td></tr>
            <tr><td style='padding:3px 8px 3px 0;font-weight:600;'>✕</td>
                <td>Close a floating panel</td></tr>
            <tr style='background:#f4f6fb;'>
                <td style='padding:3px 8px 3px 0;font-weight:600;'>Escape</td>
                <td>Exit Analyze / crosshair mode</td></tr>
          </table>
        "),
              
              div(style="padding:10px;background:#f0f8ff;border-radius:6px;
                   border-left:3px solid #2c5be0;font-size:11px;",
                  tags$a(href="https://github.com/Izzilab/NaVis", target="_blank",
                         style="color:#2c5be0;font-weight:600;text-decoration:none;",
                         "📖  Documentation & source → github.com/Izzilab/NaVis"),
                  tags$br(),
                  tags$a(href="https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0",
                         target="_blank",
                         style="color:#2c5be0;text-decoration:none;margin-top:4px;display:inline-block;",
                         "🔗  Demo dataset on 10x Genomics")
              )
              
          )
      ),
      
      div(id="fp-privacy", class="fp",
          style="top:80px; left:50%; transform:translateX(-50%); width:380px; display:none;",
          div(class="fp-header",
              span("🔒", style="font-size:13px;"),
              span("Data & Privacy", class="fp-title"),
              tags$button(class="fp-min", onclick="togglePanelMinimize('fp-privacy')", "−"),
              tags$button(class="fp-close", onclick="togglePanel('fp-privacy')", "✕")
          ),
          div(class="fp-body",
              HTML("
          <p><b>NaVis does not store, log, or transmit your data.</b>
          All processing is session-only and discarded on close.</p>
          <p>The app is hosted on ShinyApps.io. For unpublished or
          patient-derived datasets, use <b>local deployment</b>.</p>
          <p>To deploy locally: clone from
          <a href='https://github.com/Izzilab' target='_blank'>github.com/Izzilab</a>
          and run <code>shiny::runApp()</code>.</p>
        ")
          )
      ),
      
      uiOutput("autocorr_panel_ui"),
      uiOutput("legend_overlay")
      
  ), 
  
  tags$script(HTML("

function clampPanelToViewport(el) {
  if (!el || el.style.display === 'none') return;
  var vw   = window.innerWidth;
  var vh   = window.innerHeight;
  var minX = 10, minY = 80;          // leave room for the top toolbar

  var rect = el.getBoundingClientRect();
  if (rect.width > vw - 20) {
    el.style.width = (vw - 20) + 'px';
  }

  rect = el.getBoundingClientRect();

  var curLeft = rect.left;
  var curTop  = rect.top;

  var maxLeft = Math.max(minX, vw - rect.width  - minX);
  var maxTop  = Math.max(minY, vh - rect.height - 10);

  var newLeft = Math.max(minX, Math.min(curLeft, maxLeft));
  var newTop  = Math.max(minY, Math.min(curTop,  maxTop));

  if (newLeft !== curLeft || newTop !== curTop ||
      el.style.right === '' || el.style.bottom === '') {
    el.style.left      = newLeft + 'px';
    el.style.top       = newTop  + 'px';
    el.style.right     = 'auto';
    el.style.bottom    = 'auto';
    el.style.transform = 'none';
  }
}

function togglePanelMinimize(id) {
  var el = document.getElementById(id);
  if (!el) return;
  el.classList.toggle('minimized');
  var btn = el.querySelector('.fp-min');
  if (btn) btn.textContent = el.classList.contains('minimized') ? '▢' : '−';
  setTimeout(function() { clampPanelToViewport(el); }, 0);
}

function togglePanel(id) {
  var el = document.getElementById(id);
  if (!el) return;
  var vis = el.style.display === 'none' || el.style.display === '';
  el.style.display = vis ? 'flex' : 'none';
  if (vis) {

    setTimeout(function() { clampPanelToViewport(el); }, 0);
  }
}

window.addEventListener('resize', function() {
  document.querySelectorAll('.fp').forEach(function(p) {
    if (p.style.display !== 'none') clampPanelToViewport(p);
  });
});

function makeDraggable(el) {
  var header = el.querySelector('.fp-header');
  if (!header) return;
  var ox=0,oy=0,sx=0,sy=0;
  header.addEventListener('mousedown', function(e) {
    if (e.target.classList.contains('fp-close')) return;
    if (e.target.classList.contains('fp-min')) return;
    e.preventDefault();
    ox = el.offsetLeft; oy = el.offsetTop;
    sx = e.clientX;     sy = e.clientY;
    function onMove(ev) {
      var newLeft = ox + ev.clientX - sx;
      var newTop  = oy + ev.clientY - sy;

      var vw = window.innerWidth, vh = window.innerHeight;
      var w  = el.offsetWidth,    h  = el.offsetHeight;
      newLeft = Math.max(10,  Math.min(newLeft, vw - 80));          // keep 80px handle on-screen
      newTop  = Math.max(80,  Math.min(newTop,  vh - 40));          // below toolbar, above bottom
      el.style.left      = newLeft + 'px';
      el.style.top       = newTop  + 'px';
      el.style.right     = 'auto';
      el.style.bottom    = 'auto';
      el.style.transform = 'none';
    }
    function onUp() {
      document.removeEventListener('mousemove', onMove);
      document.removeEventListener('mouseup',   onUp);
    }
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup',   onUp);
  });
}

document.addEventListener('DOMContentLoaded', function() {
  document.querySelectorAll('.fp').forEach(makeDraggable);

  var observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(m) {
      m.addedNodes.forEach(function(node) {
        if (node.nodeType === 1) {
          if (node.classList && node.classList.contains('fp')) makeDraggable(node);
          node.querySelectorAll && node.querySelectorAll('.fp').forEach(makeDraggable);
        }
      });
    });
  });
  observer.observe(document.body, { childList: true, subtree: true });
});

(function() {
  var scale=1, tx=0, ty=0, dragging=false, sx=0, sy=0;
  function applyT() {
    var el = document.getElementById('zoom-inner');
    if (el) el.style.transform =
      'scale('+scale+') translate('+tx+'px,'+ty+'px)';
  }
  document.addEventListener('DOMContentLoaded', function() {
    var wrap = document.getElementById('zoom-wrap');
    if (!wrap) return;
    wrap.addEventListener('wheel', function(e) {
      e.preventDefault();
      scale = Math.min(Math.max(scale*(e.deltaY<0?1.12:0.89),0.4),16);
      applyT();
    }, {passive:false});
    wrap.addEventListener('mousedown', function(e) {
      if (analyseMode) return;
      dragging=true; sx=e.clientX-tx; sy=e.clientY-ty;
      wrap.style.cursor='grabbing';
    });
    document.addEventListener('mousemove', function(e) {
      if (!dragging) return; tx=e.clientX-sx; ty=e.clientY-sy; applyT();
    });
    document.addEventListener('mouseup', function() {
      dragging=false;
      var w=document.getElementById('zoom-wrap');
      if(w) w.style.cursor = analyseMode ? 'crosshair' : 'grab';
    });
    document.getElementById('btn-zoom-in').onclick    = function(){scale=Math.min(scale*1.25,16);applyT();};
    document.getElementById('btn-zoom-out').onclick   = function(){scale=Math.max(scale*0.8,0.4);applyT();};
    document.getElementById('btn-zoom-reset').onclick = function(){scale=1;tx=0;ty=0;applyT();};
  });
})();

Shiny.addCustomMessageHandler('selectChannel', function(msg) {
  Shiny.setInputValue('selected_img', msg.key);
});

Shiny.addCustomMessageHandler('openPanel', function(msg) {
  var el = document.getElementById(msg.id);
  if (!el) return;
  el.style.display = 'flex';
  setTimeout(function() { clampPanelToViewport(el); }, 0);
});

Shiny.addCustomMessageHandler('showEnrichContent', function(msg) {
  var hint = document.getElementById('enrich-hint');
  var ct   = document.getElementById('enrich-content');
  if (hint) hint.style.display = 'none';
  if (ct)   ct.style.display   = 'block';
});
Shiny.addCustomMessageHandler('hideEnrichContent', function(msg) {
  var hint = document.getElementById('enrich-hint');
  var ct   = document.getElementById('enrich-content');
  if (hint) hint.style.display = 'block';
  if (ct)   ct.style.display   = 'none';

  var sg = document.getElementById('simgenes-content');
  if (sg) sg.style.display = 'none';
});
Shiny.addCustomMessageHandler('showSimgenesContent', function(msg) {
  var sg = document.getElementById('simgenes-content');
  if (sg) sg.style.display = 'block';
});
Shiny.addCustomMessageHandler('hideSimgenesContent', function(msg) {
  var sg = document.getElementById('simgenes-content');
  if (sg) sg.style.display = 'none';
});
Shiny.addCustomMessageHandler('showAnalysisResults', function(msg) {
  var sp = document.getElementById('analysis-spinner');
  var ct = document.getElementById('analysis-content');
  if (sp) sp.style.display = 'none';
  if (ct) ct.style.display = 'block';
});

Shiny.addCustomMessageHandler('showCompareMarker', function(msg) {
  var stage = document.getElementById('scope-stage');
  if (!stage) return;
  var img = stage.querySelector('#zoom-inner img');
  if (!img) return;
  var imgRect = img.getBoundingClientRect();
  var stageRect = stage.getBoundingClientRect();
  var x_px = msg.x_norm * imgRect.width  + (imgRect.left - stageRect.left);
  var y_px = msg.y_norm * imgRect.height + (imgRect.top  - stageRect.top );
  var id = 'compare-marker-' + msg.n;
  var m = document.getElementById(id);
  if (!m) {
    m = document.createElement('div');
    m.id = id;
    m.className = 'compare-marker';
    m.textContent = msg.n;
    stage.appendChild(m);
  }
  m.style.left = x_px + 'px';
  m.style.top  = y_px + 'px';
  m.style.display = 'block';
});
Shiny.addCustomMessageHandler('clearCompareMarkers', function(msg) {
  var ids = ['compare-marker-1', 'compare-marker-2'];
  ids.forEach(function(id){
    var el = document.getElementById(id);
    if (el && el.parentNode) el.parentNode.removeChild(el);
  });

  var single = document.getElementById('analysis-content');
  var compare = document.getElementById('compare-content');
  if (compare) compare.style.display = 'none';

});

Shiny.addCustomMessageHandler('toggleCompareDrawer', function(msg) {
  var single  = document.getElementById('analysis-content');
  var compare = document.getElementById('compare-content');
  var hint    = document.getElementById('analysis-hint');
  var spinner = document.getElementById('analysis-spinner');
  if (msg.compare) {
    if (single)  single.style.display  = 'none';
    if (hint)    hint.style.display    = 'none';
    if (spinner) spinner.style.display = 'none';
    if (compare) compare.style.display = 'block';
  } else {
    if (compare) compare.style.display = 'none';

  }
});

Shiny.addCustomMessageHandler('toggleROIDrawer', function(msg) {
  var single  = document.getElementById('analysis-content');
  var compare = document.getElementById('compare-content');
  var roi     = document.getElementById('roi-content');
  var hint    = document.getElementById('analysis-hint');
  var spinner = document.getElementById('analysis-spinner');
  if (msg.roi) {
    if (single)  single.style.display  = 'none';
    if (compare) compare.style.display = 'none';
    if (hint)    hint.style.display    = 'none';
    if (spinner) spinner.style.display = 'none';
    if (roi)     roi.style.display     = 'block';
  } else {
    if (roi) roi.style.display = 'none';
  }
});

$(document).on('click', '.roi-acquire-btn', function(e) {
  e.preventDefault();
  e.stopPropagation();
  var gene = this.getAttribute('data-gene');
  if (!gene) return;

  Shiny.setInputValue('roi_acquire_gene',
    gene + '__' + Date.now(),
    {priority: 'event'});
});

$(document).on('shiny:recalculating', function(e) {
  if (e.target && e.target.id === 'main_image') {
    var b = document.getElementById('stage-busy-badge');
    if (b) b.style.display = 'flex';
  }
});
$(document).on('shiny:recalculated shiny:value shiny:error', function(e) {
  if (e.target && e.target.id === 'main_image') {
    var b = document.getElementById('stage-busy-badge');
    if (b) b.style.display = 'none';
  }
});

var analyseMode = false;

function toggleAnalyse() {
  analyseMode = !analyseMode;
  var btn    = document.getElementById('tb-analyse');
  var panel  = document.getElementById('fp-analysis');
  var stage  = document.getElementById('scope-stage');
  var marker = document.getElementById('click-marker');
  if (btn)    btn.classList.toggle('active', analyseMode);
  if (panel)  panel.style.display = analyseMode ? 'flex' : 'none';
  if (stage)  stage.style.cursor  = analyseMode ? 'crosshair' : 'grab';
  if (marker && !analyseMode) marker.style.display = 'none';
}

document.addEventListener('DOMContentLoaded', function() {
  var stage = document.getElementById('scope-stage');
  if (!stage) return;
  stage.addEventListener('click', function(e) {
    if (!analyseMode) return;
    if (e.shiftKey) return;   // shift+click is reserved for ROI drag start
    var img = stage.querySelector('#zoom-inner img');
    if (!img) return;
    var imgRect = img.getBoundingClientRect();
    var stageRect = stage.getBoundingClientRect();
    var x_px = e.clientX - imgRect.left;
    var y_px = e.clientY - imgRect.top;
    if (x_px<0||y_px<0||x_px>imgRect.width||y_px>imgRect.height) return;
    var x_norm = x_px/imgRect.width;
    var y_norm = y_px/imgRect.height;
    var marker = document.getElementById('click-marker');
    if (marker) {
      marker.style.display='block';
      marker.style.left=(e.clientX-stageRect.left)+'px';
      marker.style.top =(e.clientY-stageRect.top) +'px';
    }
    Shiny.setInputValue('stage_click',{x:x_norm,y:y_norm,
      px:Math.round(x_px),py:Math.round(y_px),
      iw:Math.round(imgRect.width),ih:Math.round(imgRect.height)
    },{priority:'event'});

    document.getElementById('analysis-hint').style.display='none';
    document.getElementById('analysis-content').style.display='none';
    document.getElementById('analysis-spinner').style.display='flex';
  });

  stage.addEventListener('contextmenu', function(e) {
    if (!analyseMode) return;
    var img = stage.querySelector('#zoom-inner img');
    if (!img) return;
    e.preventDefault();   // suppress browser context menu
    var imgRect = img.getBoundingClientRect();
    var x_px = e.clientX - imgRect.left;
    var y_px = e.clientY - imgRect.top;
    if (x_px<0||y_px<0||x_px>imgRect.width||y_px>imgRect.height) return;
    var x_norm = x_px/imgRect.width;
    var y_norm = y_px/imgRect.height;
    Shiny.setInputValue('stage_right_click',
      {x:x_norm, y:y_norm}, {priority:'event'});
  });

  var roiDrawing = false;
  var roiStartX, roiStartY;
  var roiBox = null;     // the in-progress rubber-band div
  var roiPersist = null; // the committed ROI rectangle (kept until cleared)

  function ensureROIBox() {
    if (!roiBox) {
      roiBox = document.createElement('div');
      roiBox.id = 'roi-rubberband';
      roiBox.style.position = 'absolute';
      roiBox.style.border = '2px dashed #2c5be0';
      roiBox.style.background = 'rgba(44,91,224,0.10)';
      roiBox.style.pointerEvents = 'none';
      roiBox.style.zIndex = '600';
      roiBox.style.display = 'none';
      stage.appendChild(roiBox);
    }
  }
  function commitROIBox(x_px, y_px, w_px, h_px) {
    if (roiPersist && roiPersist.parentNode) roiPersist.parentNode.removeChild(roiPersist);
    roiPersist = document.createElement('div');
    roiPersist.id = 'roi-overlay';
    roiPersist.style.position = 'absolute';
    roiPersist.style.border = '2px solid #2c5be0';
    roiPersist.style.background = 'rgba(44,91,224,0.05)';
    roiPersist.style.pointerEvents = 'none';
    roiPersist.style.zIndex = '500';
    roiPersist.style.left = x_px + 'px';
    roiPersist.style.top  = y_px + 'px';
    roiPersist.style.width  = w_px + 'px';
    roiPersist.style.height = h_px + 'px';
    stage.appendChild(roiPersist);
  }
  Shiny.addCustomMessageHandler('clearROIOverlay', function(msg) {
    if (roiPersist && roiPersist.parentNode) {
      roiPersist.parentNode.removeChild(roiPersist);
      roiPersist = null;
    }
  });

  stage.addEventListener('mousedown', function(e) {
    if (!analyseMode || !e.shiftKey || e.button !== 0) return;
    var img = stage.querySelector('#zoom-inner img');
    if (!img) return;
    e.preventDefault();
    e.stopPropagation();   // prevent the click handler from firing
    var imgRect = img.getBoundingClientRect();
    var stageRect = stage.getBoundingClientRect();
    var x_px = e.clientX - imgRect.left;
    var y_px = e.clientY - imgRect.top;
    if (x_px<0||y_px<0||x_px>imgRect.width||y_px>imgRect.height) return;
    roiDrawing = true;
    roiStartX = e.clientX - stageRect.left;
    roiStartY = e.clientY - stageRect.top;
    ensureROIBox();
    roiBox.style.left = roiStartX + 'px';
    roiBox.style.top  = roiStartY + 'px';
    roiBox.style.width  = '0px';
    roiBox.style.height = '0px';
    roiBox.style.display = 'block';
  }, true);  // capture phase to suppress click-to-inspect

  stage.addEventListener('mousemove', function(e) {
    if (!roiDrawing) return;
    var stageRect = stage.getBoundingClientRect();
    var cx = e.clientX - stageRect.left;
    var cy = e.clientY - stageRect.top;
    var x = Math.min(roiStartX, cx), y = Math.min(roiStartY, cy);
    var w = Math.abs(cx - roiStartX), h = Math.abs(cy - roiStartY);
    roiBox.style.left = x + 'px'; roiBox.style.top = y + 'px';
    roiBox.style.width = w + 'px'; roiBox.style.height = h + 'px';
  });

  stage.addEventListener('mouseup', function(e) {
    if (!roiDrawing) return;
    roiDrawing = false;
    var img = stage.querySelector('#zoom-inner img');
    if (!img) { roiBox.style.display = 'none'; return; }
    var imgRect = img.getBoundingClientRect();
    var stageRect = stage.getBoundingClientRect();
    var cx = e.clientX - stageRect.left;
    var cy = e.clientY - stageRect.top;

    var x1_px = roiStartX - (imgRect.left - stageRect.left);
    var y1_px = roiStartY - (imgRect.top  - stageRect.top );
    var x2_px = cx - (imgRect.left - stageRect.left);
    var y2_px = cy - (imgRect.top  - stageRect.top );

    x1_px = Math.max(0, Math.min(imgRect.width,  x1_px));
    y1_px = Math.max(0, Math.min(imgRect.height, y1_px));
    x2_px = Math.max(0, Math.min(imgRect.width,  x2_px));
    y2_px = Math.max(0, Math.min(imgRect.height, y2_px));

    var w_px = Math.abs(x2_px - x1_px), h_px = Math.abs(y2_px - y1_px);
    roiBox.style.display = 'none';
    if (w_px < 20 || h_px < 20) return;

    var sx = Math.min(roiStartX, cx);
    var sy = Math.min(roiStartY, cy);
    commitROIBox(sx, sy, Math.abs(cx - roiStartX), Math.abs(cy - roiStartY));
    Shiny.setInputValue('stage_roi', {
      x1: Math.min(x1_px, x2_px) / imgRect.width,
      y1: Math.min(y1_px, y2_px) / imgRect.height,
      x2: Math.max(x1_px, x2_px) / imgRect.width,
      y2: Math.max(y1_px, y2_px) / imgRect.height
    }, {priority:'event'});
  });

  document.addEventListener('keydown', function(e) {
    if (e.key==='Escape') {
      analyseMode=false;
      var b=document.getElementById('tb-analyse');
      var s=document.getElementById('scope-stage');
      var m=document.getElementById('click-marker');
      if(b) b.classList.remove('active');
      if(s) s.style.cursor='grab';
      if(m) m.style.display='none';
    }
  });
});

  "))
  
)  

server <- function(input, output, session) {
  
  observe_helpers()
  
  .dt <- 400   
  .tt <- 250   
  he_opacity_d       <- debounce(reactive(input$he_opacity),        .dt)
  gene_brightness_d  <- debounce(reactive(input$gene_brightness),   .dt)
  
  intensity_thr_t   <- throttle(reactive(input$intensity_threshold), .tt)
  contrastup_t      <- throttle(reactive(input$contrastup),          .tt)
  
  comp_thr1_t       <- throttle(reactive(input$comp_thr1),           .tt)
  comp_thr2_t       <- throttle(reactive(input$comp_thr2),           .tt)
  comp_contrast1_t  <- throttle(reactive(input$comp_contrast1),      .tt)
  comp_contrast2_t  <- throttle(reactive(input$comp_contrast2),      .tt)
  
  .resolve_one <- function(nm, known_genes) {
    if (is.null(nm) || is.na(nm) || !nchar(nm)) return(character(0))
    if (nm %in% known_genes) return(nm)
    if (grepl("^SIG_", nm)) {
      parts <- strsplit(sub("^SIG_", "", nm), "\\+", fixed=FALSE)[[1]]
      return(intersect(parts, known_genes))
    }
    known_genes[vapply(known_genes,
                       function(g) grepl(g, nm, fixed=TRUE),
                       logical(1))]
  }
  
  .decontaminate_spots <- function(expr, coords, genes, alpha = 0.10, k = 6L) {
    if (alpha <= 0 || length(genes) == 0) return(expr)
    
    expr_barcodes  <- colnames(expr)
    coord_barcodes <- if ("barcode" %in% colnames(coords)) coords$barcode
    else rownames(coords)
    
    shared <- if (!is.null(expr_barcodes) && !is.null(coord_barcodes))
      intersect(expr_barcodes, coord_barcodes)
    else NULL
    
    if (!is.null(shared) && length(shared) >= 5) {
      expr_idx  <- match(shared, expr_barcodes)   
      coord_idx <- match(shared, coord_barcodes)  
      coords_use <- coords[coord_idx, c("x","y")]
    } else {
      
      n_use      <- min(nrow(coords), ncol(expr))
      expr_idx   <- seq_len(n_use)
      coords_use <- coords[seq_len(n_use), c("x","y")]
    }
    
    n_spots <- nrow(coords_use)
    if (n_spots < k + 1L) return(expr)   
    k <- min(k, n_spots - 1L)
    
    nn      <- RANN::nn2(coords_use, coords_use, k = k + 1L)
    nn_idx  <- nn$nn.idx[, -1L, drop = FALSE]   
    row_i   <- rep(seq_len(n_spots), each = k)
    col_j   <- as.integer(t(nn_idx))
    W       <- Matrix::sparseMatrix(i = row_i, j = col_j,
                                    x = 1 / k,
                                    dims = c(n_spots, n_spots))
    
    sub        <- as.matrix(expr[genes, expr_idx, drop = FALSE])
    tmp        <- W %*% t(sub)          
    neigh_mean <- t(as.matrix(tmp))     
    corrected  <- pmax(sub - alpha * neigh_mean, 0)
    expr[genes, expr_idx] <- corrected
    expr
  }
  
  .auto_decontam_alpha <- function(lst, gene,
                                   alpha_grid = NULL) {
    if (!gene %in% rownames(lst$expression)) return(0.10)
    
    is_hd <- .is_hd_dataset(lst)
    if (is.null(alpha_grid)) {
      alpha_grid <- if (is_hd) seq(0, 0.20, length.out = 4)
      else        seq(0, 0.25, by = 0.05)
    }
    max_px <- if (is_hd) 200L else 300L
    
    r2_vals <- vapply(alpha_grid, function(a) {
      tryCatch({
        lst_tmp <- lst
        lst_tmp$expression <- .decontaminate_spots(
          lst$expression, lst$coords, genes = gene, alpha = a, k = 6L)
        img <- .EBIDW(lst_tmp, gene, sigma = 5, iterations = 1,
                      threshold = 0.2, shp = 20,
                      passive = FALSE, normalize = TRUE)
        
        gmat <- as.matrix(EBImage::normalize(img))
        sc <- max_px / max(dim(gmat))
        if (sc < 1) {
          g_img <- EBImage::resize(
            EBImage::Image(gmat, colormode = "Grayscale"),
            w = round(dim(gmat)[1] * sc),
            h = round(dim(gmat)[2] * sc))
          gmat <- as.matrix(imageData(g_img))
        }
        td <- dim(gmat)
        prior_nms <- c("nuclei_mask", "fibrous_tissue_mask", "soft_tissue_mask")
        pm_list <- lapply(prior_nms, function(nm) .get_prior_mat(nm, lst, td))
        valid   <- !vapply(pm_list, is.null, logical(1))
        if (sum(valid) < 1) return(NA_real_)
        y <- as.numeric(gmat)
        X <- do.call(cbind, lapply(pm_list[valid], as.numeric))
        keep <- is.finite(y) & apply(X, 1, function(r) all(is.finite(r)))
        if (sum(keep) < 50L) return(NA_real_)
        y_k <- y[keep]; X_k <- X[keep, , drop = FALSE]
        ss_tot <- sum((y_k - mean(y_k))^2)
        fit <- lm.fit(cbind(1, X_k), y_k)
        1 - sum(fit$residuals^2) / (ss_tot + 1e-8)
      }, error = function(e) NA_real_)
    }, numeric(1))
    
    ok  <- is.finite(r2_vals)
    if (sum(ok) < 3) return(0.10)
    av  <- alpha_grid[ok];  rv <- r2_vals[ok]
    
    a_n <- (av - av[1])  / (av[length(av)]  - av[1]  + 1e-8)
    r_n <- (rv - min(rv)) / (max(rv) - min(rv) + 1e-8)
    
    chord  <- r_n[1] + a_n * (r_n[length(r_n)] - r_n[1])
    dists  <- abs(r_n - chord)
    best   <- av[which.max(dists)]
    round(best, 2)
  }
  
  auto_decontam_result <- reactiveVal(NULL)
  
  output$auto_decontam_status_ui <- renderUI({
    res <- auto_decontam_result()
    if (is.null(res)) return(NULL)
    div(style=sprintf("font-size:10px;font-weight:700;color:%s;white-space:nowrap;",
                      if (res$ok) "#1e7a45" else "#c85050"),
        if (res$ok) sprintf("→ α = %.2f", res$alpha) else "failed")
  })
  outputOptions(output, "auto_decontam_status_ui", suspendWhenHidden = FALSE)
  
  observeEvent(input$auto_decontam, {
    req(data_list(), input$selected_gene)
    genes_sel <- input$selected_gene
    req(length(genes_sel) > 0)
    
    gene_for_cal <- genes_sel[genes_sel %in% data_list()$rn][1]
    req(!is.na(gene_for_cal))
    
    withProgress(message = sprintf("Auto-detecting α for %s…", gene_for_cal),
                 value = 0, {
                   incProgress(0.2)
                   best_alpha <- tryCatch(
                     .auto_decontam_alpha(isolate(data_list()), gene_for_cal),
                     error = function(e) NULL)
                   incProgress(1)
                   if (!is.null(best_alpha) && is.finite(best_alpha)) {
                     updateSliderInput(session, "decontam_alpha", value = best_alpha)
                     auto_decontam_result(list(ok = TRUE, alpha = best_alpha))
                     showNotification(
                       sprintf("Optimal α = %.2f for %s (architecture-coupling elbow)",
                               best_alpha, gene_for_cal),
                       type = "message", duration = 5)
                   } else {
                     auto_decontam_result(list(ok = FALSE, alpha = NA))
                     showNotification("Auto-detection failed — using current α.",
                                      type = "warning", duration = 4)
                   }
                 })
  })
  
  session$onFlushed(function() {
    session$sendCustomMessage("openPanel", list(id = "fp-help"))
  }, once = TRUE)
  
  observeEvent(input$load_demo, {
    showModal(modalDialog(
      title = tags$div(
        style = "display:flex;align-items:center;gap:10px;",
        tags$span("📦", style="font-size:20px;"),
        tags$span("About the demo dataset",
                  style="font-size:17px;font-weight:700;color:#2a3050;")
      ),
      tags$div(
        style = "font-size:13px;line-height:1.75;color:#212529;",
        tags$p(
          tags$strong("Human breast cancer — ductal carcinoma "),
          tags$em("in situ"), tags$strong(" (DCIS) with invasive carcinoma, FFPE."),
          " 10x Genomics Visium Spatial Gene Expression, Space Ranger 1.3.0."
        ),
        tags$p(
          "This section contains intact DCIS ducts, foci of invasive ductal carcinoma,
           stromal compartments, and patches of immune-infiltrated tissue.
           It is ideal for NaVis because the H&E encodes clearly-resolved compartment
           boundaries — the myoepithelial ring, peri-ductal stroma, and invasive front
           are all visible — providing a meaningful ground truth against which interpolated
           expression maps and image-analytic outputs can be visually verified."
        ),
        tags$div(
          style="background:#F0FAF9;border-left:3px solid #2A9D8F;
                 padding:8px 12px;margin-top:10px;border-radius:4px;",
          tags$a(
            href="https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0",
            target="_blank",
            style="color:#2A9D8F;font-weight:600;text-decoration:none;",
            "↗  Original dataset on 10x Genomics"
          )
        ),
        tags$p(
          style="margin-top:12px;color:#6C757D;font-size:11px;font-style:italic;",
          "Loading in the background — close this dialog when you are ready."
        )
      ),
      footer = tagList(modalButton("Got it")),
      size = "m", easyClose = TRUE, fade = TRUE
    ))
    
    demo_dir <- "www/demo"
    demo_rds <- list.files(demo_dir,
                           pattern = "\\.(rds|Rds|RDS)$",
                           full.names = TRUE)
    if (length(demo_rds) == 0) {
      removeModal()
      showModal(modalDialog(
        title = "Demo not available",
        tags$p("No demo object was found in the ", tags$code("www/demo/"),
               " folder. Please add a NaVis/Seurat ", tags$code(".rds"),
               " file there and restart the app."),
        footer = modalButton("OK"), size = "m", easyClose = TRUE
      ))
      return()
    }
    demo_file <- demo_rds[1]
    
    withProgress(message = "Loading demo dataset…", value = 0, {
      incProgress(0.1, detail = basename(demo_file))
      obj <- tryCatch(readRDS(demo_file), error = function(e) NULL)
      if (is.null(obj)) {
        removeModal()
        showModal(modalDialog(
          title = "Failed to load demo",
          tags$p("The demo .rds file could not be read. Try re-downloading it."),
          footer = modalButton("OK"), size = "m", easyClose = TRUE
        ))
        return()
      }
      
      if (inherits(obj, "Seurat")) {
        incProgress(0.4, detail = "Converting Seurat → NaVis")
        obj <- tryCatch(seurat_to_genesis(obj), error = function(e) NULL)
        if (is.null(obj)) {
          removeModal()
          showModal(modalDialog(
            title = "Conversion failed",
            tags$p("The demo Seurat object could not be converted to the NaVis format."),
            footer = modalButton("OK"), size = "m", easyClose = TRUE
          ))
          return()
        }
        odt("seu")
      } else {
        odt("noseu")
      }
      incProgress(0.7, detail = "Publishing reactive state")
      
      orig(obj)
      incProgress(1, detail = "Done")
    })
  })
  
  
  shinyjs::disable("imgcust")
  shinyjs::show("autocorr_panel")
  
  
  cntr <- reactive({
    v <- contrastup_t()
    if (is.null(v) || is.na(v)) return(1)
    exp(-1.9 * v)
  })
  
  
  observeEvent(input$localtex, {
    if (isTRUE(input$localtex)) {
      shinyjs::enable("imgcust")
    } else {
      shinyjs::disable("imgcust")
    }
  })
  
  
  nmm <- reactiveVal()
  
  
  img_dir <- file.path(tempdir(), paste0("served_images_", session$token))
  dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
  shiny::addResourcePath("served_images", img_dir)
  
  
  imgs <- reactiveVal(NULL)
  
  
  odt <- reactiveVal()
  orig <- reactiveVal()
  
  
  vizdir <- reactiveVal()
  vizdir_is_hd <- reactiveVal(FALSE)
  vizdir_hd_bin <- reactiveVal(NA_character_)
  step <- reactiveVal(1)
  
  observeEvent(input$'next', { step(step() + 1) })
  observeEvent(input$prev, { step(step() - 1) })
  
  output$wizard_ui <- renderUI({
    
    if (step() == 1) {
      tagList(
        h3("Step 1 — Choose spatial raw format"),
        radioButtons("raw_format", "Raw spatial format:",
                     choices = c("10x Visium or compatible (zipped)" = "visium_zip")
        ),
        actionButton("next", "Next")
      )
    }
    
    else if (step() == 2) {
      tagList(
        h3("Step 2 — Upload ZIP containing the raw data and validate"),
        fileInput("raw_zip", "Upload ZIP",
                  accept = ".zip", multiple = FALSE),
        p("Your ZIP must contain the 'spatial/' directory and the
           'filtered_feature_bc_matrix/' directory."),
        actionButton("prev", "Back"),
        actionButton("next", "Validate ZIP")
      )
    }
    
    else if (step() == 3) {
      
      req(input$raw_zip)
      
      tmpdir <- tempfile()
      dir.create(tmpdir)
      unzip(input$raw_zip$datapath, exdir = tmpdir)
      nmm(input$raw_zip)
      
      files <- list.files(tmpdir, recursive = TRUE)
      
      
      has_mtx <- function(files) {
        all(c(
          "filtered_feature_bc_matrix/barcodes.tsv.gz",
          "filtered_feature_bc_matrix/features.tsv.gz",
          "filtered_feature_bc_matrix/matrix.mtx.gz"
        ) %in% files)
      }
      
      has_h5 <- function(files) {
        any(c(
          "filtered_feature_bc_matrix.h5",
          "filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
          "raw_feature_bc_matrix.h5"
        ) %in% files)
      }
      
      spatial_variants <- list(
        
        visium_standard = c(
          "spatial/tissue_hires_image.png",
          "spatial/tissue_lowres_image.png",
          "spatial/scalefactors_json.json",
          "spatial/tissue_positions_list.csv"
        ),
        
        visium_csv = c(
          "spatial/tissue_hires_image.png",
          "spatial/tissue_lowres_image.png",
          "spatial/scalefactors_json.json",
          "spatial/tissue_positions.csv"
        ),
        
        visium_parquet = c(
          "spatial/tissue_hires_image.png",
          "spatial/tissue_lowres_image.png",
          "spatial/scalefactors_json.json",
          "spatial/tissue_positions.parquet"
        ),
        
        cytassist = c(
          "spatial/cytassist_image.tiff",
          "spatial/scalefactors_json.json",
          "spatial/tissue_positions.csv"
        )
      )
      
      .detect_hd_bins <- function(file_list) {
        rx <- "^binned_outputs/square_(\\d+)um/spatial/tissue_positions\\.parquet$"
        m  <- regmatches(file_list, regexpr(rx, file_list))
        if (!length(m)) return(character(0))
        bins <- unique(sub(rx, "\\1", m))
        bins[order(as.integer(bins))]   
      }
      hd_bins <- .detect_hd_bins(files)
      is_hd   <- length(hd_bins) > 0
      
      visium_ok <- FALSE
      detected_type <- NULL
      
      if (is_hd) {
        
        usable_bins <- intersect(hd_bins, c("008", "016"))
        if (!length(usable_bins)) {
          return(tagList(
            h3("Error — Visium HD detected but no usable bin level"),
            p("This dataset only contains the 2 µm bin level (",
              tags$code("square_002um"),
              "), which has ~11 million squares per section and is not ",
              "supported for interactive analysis. Re-run Space Ranger ",
              "with default settings to also produce 8 µm and 16 µm bins, ",
              "or aggregate the 2 µm bins offline before uploading."),
            actionButton("prev", "Back")
          ))
        }
        
        bin_has_expr <- function(b) {
          base <- sprintf("binned_outputs/square_%sum/", b)
          any(c(
            paste0(base, "filtered_feature_bc_matrix.h5"),
            paste0(base, "filtered_feature_bc_matrix/matrix.mtx.gz"),
            paste0(base, "raw_feature_bc_matrix.h5")
          ) %in% files)
        }
        usable_bins <- usable_bins[vapply(usable_bins, bin_has_expr, TRUE)]
        if (!length(usable_bins)) {
          return(tagList(
            h3("Error — Visium HD bins detected but no expression matrix found"),
            p("The 8 µm and 16 µm spatial folders were detected, but no ",
              tags$code("filtered_feature_bc_matrix.h5"),
              " or matrix.mtx.gz was present at the same bin level. ",
              "Please include the expression matrices in the ZIP."),
            actionButton("prev", "Back")
          ))
        }
        
        chosen_bin <- if ("016" %in% usable_bins) "016" else usable_bins[1]
        visium_ok    <- TRUE
        detected_type <- sprintf("visium_hd_%sum", chosen_bin)
      } else {
        has_expression <- has_mtx(files) || has_h5(files)
        if (has_expression) {
          for (name in names(spatial_variants)) {
            if (all(spatial_variants[[name]] %in% files)) {
              visium_ok <- TRUE
              detected_type <- name
              break
            }
          }
        }
      }
      
      if (!visium_ok) {
        return(tagList(
          h3("Error — ZIP is not a valid Visium directory"),
          p("Required Visium spatial structure and/or expression matrix not detected."),
          actionButton("prev", "Back")
        ))
      }
      
      if (is_hd) {
        visium_dir <- file.path(tmpdir, "binned_outputs",
                                sprintf("square_%sum", chosen_bin))
      } else {
        visium_dir <- tmpdir
      }
      vizdir(visium_dir)
      vizdir_is_hd(is_hd)
      vizdir_hd_bin(if (is_hd) chosen_bin else NA_character_)
      
      validated_path <- renderText(visium_dir)
      
      tagList(
        h3("Step 3 — ZIP validated"),
        p("All required Visium files were detected. Please proceed to import data!"),
        verbatimTextOutput("validated_path"),
        actionButton("prev", "Back"),
        actionButton("next", "Load Data")
      )
    }
    
    else if (step() == 4) {
      tagList(
        h3("Step 4 — Importing Raw Data"),
        p("Click Import to create a Navis object."),
        actionButton("prev", "Back"),
        actionButton("import_raw", "Import Raw Data")
      )
    }
  })
  
  observeEvent(input$import_raw, {
    
    is_hd  <- isolate(vizdir_is_hd())
    hd_bin <- isolate(vizdir_hd_bin())
    
    if (is_hd) {
      showModal(modalDialog(sprintf(
        "Loading Visium HD data at %s µm bin level — this may take longer than standard Visium...",
        hd_bin)))
    } else {
      showModal(modalDialog("Loading Spatial Data and tranforming into NaVis format…"))
    }
    
    visium_dir <- isolate(vizdir())
    print(visium_dir)
    print(list.files(visium_dir))
    
    seu <- BuildGenesis(visium_dir)
    
    if (is_hd) {
      bin_um <- as.numeric(hd_bin)
      
      if (!is.null(seu$spot_diameter) && seu$spot_diameter > 0) {
        seu$microns_per_pixel <- bin_um / seu$spot_diameter
      }
      seu$visium_hd_bin_um <- bin_um
    }
    
    removeModal()
    step(1)
    
    odt("noseu")
    orig(seu)
  })
  
  
  
  .has_he_image <- function(obj) {
    if (is.null(obj)) return(FALSE)
    
    if (isS4(obj) && "images" %in% slotNames(obj)) {
      imgs_slot <- tryCatch(obj@images, error = function(e) NULL)
      if (is.null(imgs_slot) || length(imgs_slot) == 0) return(FALSE)
      img_obj <- imgs_slot[[1]]
      if (is.null(img_obj)) return(FALSE)
      if (isS4(img_obj) && "image" %in% slotNames(img_obj)) {
        m <- tryCatch(img_obj@image, error = function(e) NULL)
        return(!is.null(m) && length(m) > 100)
      }
      return(TRUE)
    }
    
    if (is.list(obj)) {
      if (!is.null(obj$image) && length(obj$image) > 100) return(TRUE)
    }
    FALSE
  }
  
  .is_hd_dataset <- function(lst = NULL) {
    if (is.null(lst)) return(FALSE)
    if (!is.list(lst)) return(FALSE)
    bin <- lst$visium_hd_bin_um
    !is.null(bin) && is.finite(bin) && bin > 0
  }
  
  observeEvent(input$proc_file, {
    withProgress(message="Reading data in...", value=0,{
      nmm(input$proc_file)
      q <- readRDS(input$proc_file$datapath)
      
      if (!.has_he_image(q)) {
        showModal(modalDialog(
          title = tags$span(style="color:#c85050;", "⚠ No H&E image found"),
          tags$p("This file does not contain an embedded H&E (histology) image."),
          tags$p(style="font-size:12px;color:#5a4010;",
                 "NaVis requires a Visium object built with the H&E image attached. ",
                 "If you have a Seurat object, rebuild it from the Space Ranger ",
                 "output using ", tags$code("Load10X_Spatial()"),
                 " — or upload the Space Ranger output folder directly via the ",
                 tags$b("Raw"), " tab."),
          easyClose = TRUE, footer = modalButton("Close")
        ))
        return(invisible(NULL))
      }
      
      if(class(q)=="Seurat"){
        v <- as.numeric(strsplit(paste(q@version,collapse=""),"\\.")[[1]][1])
        if(v<4){
          stop("Seurat must be version 4 or higher")
        }else{
          incProgress(.5)
          q <- seurat_to_genesis(q)
          odt("seu")
          orig(q)
        }
      }else{
        if(class(q)=="list" || class(q)=="genesis")
          odt("noseu")
        orig(q)
      }
      incProgress(1)
    })
  })
  
  
  ia <- reactiveVal()
  
  load_id <- reactiveVal(0L)
  
  data_list <- reactive({
    req(orig())
    d <- orig()
    
    
    load_id(isolate(load_id()) + 1L)
    
    if(length(d)>=20){
      imgs(NULL)
      ia(0)
      d
    }else{
      withProgress(message="Performing image analysis...", value=0,{
        if(odt()=="noseu"){
          incProgress(0.5,message = "fixing histogram drift")
          om <- d$image
          mom <- suppressMessages(suppressWarnings(hist_match(om, hv)))
          d$image <- mom
          ia(1)
          incProgress(0.75,message = "image deconvolution")
          dat <- GenesisImageAnalysis(d)
        }else{
          m <- image_read(d$image)
          unsharp_kernel <- matrix(
            c(-1,-1,-1,
              -1,  9,-1,
              -1,-1,-1),
            nrow = 3, byrow = TRUE
          )
          incProgress(0.3,message = "sharpening low-res image")
          m <- image_convolve(m, unsharp_kernel)
          img_array <- as.numeric(image_data(m, channels = "rgb"))
          incProgress(0.5,message = "fixing histogram drift")
          om <- img_array
          mom <- suppressMessages(suppressWarnings(hist_match(om, hv)))
          d$image <- mom
          ia(1)
          incProgress(0.75,message = "image deconvolution")
          dat <- GenesisImageAnalysis(d)
        }
        incProgress(1)
        d <- c(d,dat)
        imgs(NULL)
        d
      })
    }
  })
  
  
  observe({
    d <- data_list()
    updateSelectizeInput(session, "selected_gene",
                         choices = d$rn, server = TRUE)
  })
  
  observe({
    d <- data_list()
    req(d$meta)
    meta_cols <- colnames(d$meta)
    meta_cols <- meta_cols[meta_cols != "no.metadata"]
    updateSelectInput(session, "meta_col", choices = meta_cols)
  })
  
  observeEvent(data_list(), {
    lst <- data_list()
    targets <- c("original_image","entropy_mask", "soft_tissue_mask",
                 "fibrous_tissue_mask", "nuclei_mask")
    available <- intersect(targets, names(lst))
    
    files <- withProgress(message="Preparing images...", value=0, {
      out <- lapply(seq_along(available), function(i) {
        nm <- available[i]
        incProgress(1/length(available), detail=paste("Processing", nm))
        save_prior_png(lst[[nm]], nm, img_dir)
      })
      out
    })
    
    names(files) <- c("H&E Image","Entropy", "Non-fibrillar interstitial tissue",
                      "Fibrillar interstitial tissue", "Nuclei and cells")
    
    imgs(files)
    
    session$sendCustomMessage("openPanel", list(id = "fp-channels"))
    session$sendCustomMessage("selectChannel", list(key = "H&E Image"))
  })
  
  
  prms <- reactiveVal()
  primg <- reactiveVal()
  
  observe({
    req(data_list())
    
    if(isTRUE(input$auto_mode)) {
      if(nrow(orig()$coords)<=20000){
        if(odt() == "seu"){
          prms(c(5, 1, 0.2, 10))
        }else{
          prms(c(10, 3, 0.5, 30))
        }
      }else{
        prms(c(10, 3, 0, 30)) 
      }
    } else {
      prms(c(input$sigma, input$iterations, input$threshold, input$sharpness))
      primg(input$imgcust)
    }
  })
  
  
  gnm <- reactiveVal()
  
  
  observe({
    sel <- input$selected_img
    is_meta <- !is.null(sel) && grepl("^META_", sel)
    if (is_meta) shinyjs::disable("gene_color")
    else         shinyjs::enable("gene_color")
  })
  
  observe({
    d <- data_list()
    req(d, input$meta_col)
    if (!(input$meta_col %in% colnames(d$meta))) return()
    lv <- d$meta[[input$meta_col]]
    if (is.null(lv)) lv <- d$meta[, input$meta_col, drop = TRUE]
    if (!is.null(dim(lv))) lv <- as.vector(lv[, 1, drop = TRUE])
    if (is.numeric(lv)) shinyjs::enable("meta_color")
    else                shinyjs::disable("meta_color")
  })
  
  
  img_legends <- reactiveVal(list())
  lgnd <- reactiveVal(NULL)
  
  
  mll <- reactiveVal()
  
  
  observeEvent(input$gen_gene, {
    req(data_list())
    req(input$selected_gene)
    req(length(input$selected_gene) > 0)
    
    genes_sel  <- input$selected_gene   
    is_sig     <- length(genes_sel) > 1
    sig_label  <- if (is_sig) paste0("SIG_", paste(genes_sel, collapse="+")) else genes_sel[1]
    
    gnm(c(sig_label, gnm()))
    lst <- data_list()
    
    gene_colored <- withProgress(
      message = if (is_sig) paste0("Building signature (", length(genes_sel), " genes)…")
      else "Calculating new variable image...", value=0, {
        
        incProgress(0.1)
        pv <- prms()
        
        valid_genes <- genes_sel[genes_sel %in% lst$rn]
        if (length(valid_genes) == 0) stop("None of the selected genes found in dataset.")
        
        alpha_dc <- if (!is.null(input$decontam_alpha)) input$decontam_alpha else 0
        if (alpha_dc > 0 && !is.null(lst$coords) && length(valid_genes) > 0) {
          incProgress(0.15, detail = sprintf("Leakage correction (α=%.2f)…", alpha_dc))
          lst$expression <- .decontaminate_spots(
            lst$expression,
            lst$coords,
            genes  = valid_genes,
            alpha  = alpha_dc,
            k      = 6L
          )
        }
        
        spot_coexpr_mask <- NULL   
        
        dropout_tau <- if (!is.null(input$dropout_tau)) input$dropout_tau else 0.5
        
        if (is_sig && length(valid_genes) >= 2) {
          coords_df <- as.data.frame(lst$coords)
          if (!is.null(coords_df) && nrow(coords_df) >= 7) {
            k_neigh <- 6L
            
            bc_expr  <- colnames(lst$expression)
            bc_coord <- if ("barcode" %in% colnames(coords_df)) coords_df$barcode
            else rownames(coords_df)
            
            coords_aligned <- NULL
            expr_idx_keep  <- seq_len(ncol(lst$expression))  
            
            if (!is.null(bc_expr) && !is.null(bc_coord)) {
              
              row_of <- match(bc_expr, bc_coord)       
              good   <- !is.na(row_of)
              if (sum(good) >= 7) {
                coords_aligned <- coords_df[row_of[good], c("x", "y"), drop = FALSE]
                expr_idx_keep  <- which(good)          
              }
            }
            
            if (is.null(coords_aligned) && nrow(coords_df) == ncol(lst$expression)) {
              coords_aligned <- coords_df[, c("x", "y"), drop = FALSE]
              expr_idx_keep  <- seq_len(ncol(lst$expression))
            }
            
            if (!is.null(coords_aligned) && nrow(coords_aligned) >= 7) {
              nn     <- RANN::nn2(coords_aligned, coords_aligned, k = k_neigh + 1L)
              nn_idx <- nn$nn.idx[, -1L, drop = FALSE]      
              n_aligned <- nrow(coords_aligned)             
              
              truezero_list <- vector("list", length(valid_genes))
              names(truezero_list) <- valid_genes
              
              for (g in valid_genes) {
                v_full <- as.numeric(lst$expression[g, ])        
                
                v <- v_full[expr_idx_keep]
                is_zero <- v == 0
                
                neigh_vals <- matrix(v[nn_idx], nrow = n_aligned, ncol = k_neigh)
                neigh_zero <- rowSums(neigh_vals == 0) / k_neigh
                true_zero  <- is_zero & (neigh_zero >= dropout_tau)
                is_dropout <- is_zero & !true_zero
                
                if (any(is_dropout)) {
                  drop_rows <- which(is_dropout)
                  imp_val <- vapply(drop_rows, function(i) {
                    nvals <- v[nn_idx[i, ]]
                    nvals <- nvals[nvals > 0]
                    if (length(nvals) > 0) mean(nvals) else 0
                  }, numeric(1))
                  v[drop_rows] <- imp_val
                  
                  v_full[expr_idx_keep] <- v
                  lst$expression[g, ] <- v_full
                }
                
                truezero_list[[g]] <- true_zero
              }
              
              tz_mat   <- do.call(cbind, truezero_list)         
              all_true_zero  <- rowSums(tz_mat) == ncol(tz_mat) 
              coexpr_aligned <- !all_true_zero                   
              
              coexpr_spots <- rep(TRUE, ncol(lst$expression))
              coexpr_spots[expr_idx_keep] <- coexpr_aligned
              
              if (any(!coexpr_spots)) {
                for (g in valid_genes) {
                  lst$expression[g, !coexpr_spots] <- 0
                }
              }
            }
          }
        }
        
        make_map_gated <- function(gene) {
          if (pv[1] == 0) {
            a <- .EBIDW(lst, gene, sigma=10, iterations=pv[2], threshold=pv[3],
                        shp=pv[4], passive=TRUE,
                        datatype=if(odt()=="seu") "seu" else "noseu")
          } else if (isTRUE(input$mlbp)) {
            a <- .EBIDW(lst, gene, sigma=7, iterations=pv[2], threshold=0.2,
                        shp=pv[4], passive=FALSE, normalize=FALSE)  
          } else {
            a <- .EBIDW(lst, gene, sigma=pv[1], iterations=pv[2], threshold=pv[3],
                        shp=pv[4], passive=FALSE, normalize=FALSE)
          }
          
          a_mat <- as.matrix(a)
          
          if (!is.null(spot_coexpr_mask)) {
            
            if (all(dim(a_mat) == dim(spot_coexpr_mask))) {
              a_mat[!spot_coexpr_mask] <- 0
            }
          }
          a <- EBImage::Image(a_mat, colormode = "Grayscale")
          EBImage::normalize(a)
        }
        
        if (is_sig && length(valid_genes) >= 2 && is.null(spot_coexpr_mask)) {
          img_h <- dim(lst$image)[1]; img_w <- dim(lst$image)[2]
          coords_df <- as.data.frame(lst$coords)
          xi <- pmin(pmax(round(coords_df[, ncol(coords_df)-1]), 1), img_w)
          yi <- pmin(pmax(round(coords_df[, ncol(coords_df)]),   1), img_h)
          gene_vecs2 <- lapply(valid_genes, function(g)
            as.numeric(lst$expression[g, ]))
          
          zero_mat_l2  <- do.call(cbind, lapply(gene_vecs2, function(v) v == 0))
          all_zero_l2  <- rowSums(zero_mat_l2) == ncol(zero_mat_l2)
          coexpr2      <- !all_zero_l2
          
          gate_mat <- matrix(0L, nrow = img_h, ncol = img_w)
          gate_mat[cbind(yi[!coexpr2], xi[!coexpr2])] <- 1L
          
          if (pv[1] == 0) {
            base_sigma <- 10; base_iter <- 1          
          } else if (isTRUE(input$mlbp)) {
            base_sigma <- 7;  base_iter <- pv[2]
          } else {
            base_sigma <- pv[1]; base_iter <- pv[2]
          }
          combined_sigma <- sqrt(max(base_iter, 1)) * base_sigma
          dil_r <- max(1L, as.integer(round(combined_sigma * 2)))
          kern  <- EBImage::makeBrush(as.integer(2 * dil_r + 1), shape = "disc")
          gate_dilated <- EBImage::dilate(
            EBImage::Image(gate_mat + 0.0, colormode="Grayscale"), kern)
          
          gate_dilated <- EBImage::rotate(EBImage::Image(as.matrix(imageData(gate_dilated)),
                                                         colormode="Grayscale"), -90)
          gate_dilated <- EBImage::flop(gate_dilated)
          spot_coexpr_mask <- as.matrix(imageData(gate_dilated)) < 0.5   
        }
        
        maps <- lapply(valid_genes, make_map_gated)
        incProgress(0.4)
        
        if (length(maps) == 1) {
          a <- maps[[1]]
        } else {
          
          arr_list <- lapply(maps, function(m) as.matrix(EBImage::normalize(m)))
          a_mat    <- Reduce("+", arr_list) / length(arr_list)
          
          eps <- 0.01
          pix_pos_count <- matrix(0L, nrow = nrow(a_mat), ncol = ncol(a_mat))
          for (m in arr_list) pix_pos_count <- pix_pos_count + (m > eps)
          pix_coexpr <- pix_pos_count >= 1L
          a_mat[!pix_coexpr] <- 0
          
          a <- EBImage::Image(a_mat, colormode = "Grayscale")
        }
        
        incProgress(0.5)
        
        if (pv[1] != 0) {
          z  <- if (isTRUE(input$auto_mode)) EBImage::channel(lst$original_image,"gray")
          else EBImage::channel(EBImage::readImage(imgs()[[primg()]]), "gray")
          
          
          .conform_to <- function(arr_or_img, ref_img) {
            ref_dims <- dim(ref_img)
            cur_dims <- dim(arr_or_img)
            if (is.null(cur_dims) || is.null(ref_dims)) return(arr_or_img)
            if (length(cur_dims) >= 2 && length(ref_dims) >= 2 &&
                (cur_dims[1] != ref_dims[1] || cur_dims[2] != ref_dims[2])) {
              if (inherits(arr_or_img, "Image")) {
                return(EBImage::resize(arr_or_img, w = ref_dims[1], h = ref_dims[2]))
              } else {
                tmp <- EBImage::Image(as.matrix(arr_or_img), colormode = "Grayscale")
                tmp <- EBImage::resize(tmp, w = ref_dims[1], h = ref_dims[2])
                return(EBImage::imageData(tmp))
              }
            }
            arr_or_img
          }
          
          if (isTRUE(input$mlbp)) {
            incProgress(0.6, message = "ML boost and prune")
            learn_and_prune_auto <- function(img_rgb, mask_nuc, mask_fib, mask_soft,
                                             prune_strength = 1, verbose = TRUE) {
              img_intensity <- Image(img_rgb,colormode = "grayscale")
              mask_nuc <- imageData(mask_nuc)
              mask_nuc[mask_nuc>0] <- 1
              X <- cbind(as.numeric(mask_nuc),as.numeric(mask_fib),as.numeric(mask_soft))
              y <- as.numeric(img_intensity)
              X[is.na(X)] <- 0; y[is.na(y)] <- 0
              muX <- colMeans(X); sX <- apply(X,2,sd)+1e-8
              Xs  <- sweep(sweep(X,2,muX,"-"),2,sX,"/")
              loss_fn <- function(beta) mean((y-Xs%*%beta)^2)+0.01*sum(beta^2)
              opt  <- optim(rep(0,ncol(Xs)),loss_fn,method="L-BFGS-B",lower=-5,upper=5)
              beta <- opt$par
              mll(paste("Learned coefficients:",paste(round(beta,3),collapse=", ")))
              pred_map  <- matrix(as.numeric(Xs%*%beta),nrow=nrow(img_intensity))
              pred_map  <- (pred_map-min(pred_map,na.rm=TRUE))/(max(pred_map,na.rm=TRUE)-min(pred_map,na.rm=TRUE)+1e-8)
              res_mask  <- pmax(0,pred_map-mean(pred_map))
              res_mask  <- (res_mask/(max(res_mask)+1e-8))^1.2
              res_mask  <- (1-prune_strength)+prune_strength*res_mask
              pmin(pmax(img_rgb*res_mask,0),1)
            }
            k  <- learn_and_prune_auto(a,lst$nuclei_mask,lst$fibrous_tissue_mask,lst$soft_tissue_mask)
            zs <- gblur(z,sigma=4); local_var <- abs(z-zs)
            decay_mask <- exp(-4*local_var)
            decay_mask <- .conform_to(decay_mask, k)
            k  <- k*(1-decay_mask); k <- normalize(k)^0.45
          } else {
            zs <- gblur(z,sigma=4); local_var <- abs(z-zs)
            decay_mask <- exp(-4*local_var)
            decay_mask <- .conform_to(decay_mask, a)
            k  <- a*(1-decay_mask); k <- normalize(k)^0.45
          }
        } else {
          k <- a
        }
        
        incProgress(0.75)
        
        g   <- col_gene(k, input$gene_color)
        msk <- 1 - lst$white_mask
        d_img <- g$image
        
        if (!is.null(dim(d_img)) && length(dim(d_img)) >= 2 &&
            !is.null(dim(msk)) && length(dim(msk)) >= 2 &&
            (dim(msk)[1] != dim(d_img)[1] || dim(msk)[2] != dim(d_img)[2])) {
          tmp <- EBImage::Image(as.matrix(msk), colormode = "Grayscale")
          tmp <- EBImage::resize(tmp, w = dim(d_img)[1], h = dim(d_img)[2])
          msk <- EBImage::imageData(tmp)
        }
        
        d_img[,,1] <- d_img[,,1]*msk
        d_img[,,2] <- d_img[,,2]*msk
        d_img[,,3] <- d_img[,,3]*msk
        g$image <- d_img
        incProgress(1)
        g
      })
    
    v <- gnm()
    if(length(v[v %in% sig_label]) > 1){
      base_name <- paste0(sig_label,"_",length(v[v %in% sig_label]))
    } else {
      base_name <- sig_label
    }
    
    fname <- file.path(img_dir, paste0(base_name, ".png"))
    
    writeImage(gene_colored$image, fname)
    
    lgnd(gene_colored$legend)
    
    current <- imgs()
    if(is.null(current)) current <- list()
    current[[base_name]] <- fname
    imgs(current)
    
    lgs <- img_legends()
    lgs[[base_name]] <- gene_colored$legend
    img_legends(lgs)
    
    session$sendCustomMessage("openPanel",    list(id  = "fp-channels"))
    session$sendCustomMessage("selectChannel", list(key = base_name))
  })
  
  output$channel_rail_ui <- renderUI({
    req(imgs())
    files <- imgs()
    selected <- input$selected_img
    
    channel_colors <- list(
      "H&E Image" = "#aaaaaa",
      "Entropy"   = "#7766ff",
      "Non-fibrillar interstitial tissue" = "#44bb88",
      "Fibrillar interstitial tissue"     = "#ff8844",
      "Nuclei and cells"                  = "#4488ff"
    )
    
    tagList(lapply(names(files), function(nm) {
      
      rel_path <- paste0("served_images/", basename(files[[nm]]),
                         "?v=", load_id())
      is_active <- !is.null(selected) && selected == nm
      dot_color <- if (nm %in% names(channel_colors)) channel_colors[[nm]] else "#ff3366"
      
      tags$div(
        class = paste("channel-card", if (is_active) "active" else ""),
        onclick = sprintf("Shiny.setInputValue('selected_img','%s')", nm),
        tags$img(src = rel_path, class = "channel-thumb"),
        tags$span(nm, class = "channel-label"),
        tags$div(class = "channel-dot",
                 style = paste0("background:", dot_color, ";"))
      )
    }))
  })
  
  
  
  
  trl <- reactiveVal()
  crispr <- reactiveVal(0)
  last_rendered_pixels <- reactiveVal(NULL)
  
  observeEvent(input$crisp,{
    crispr(1)
  })
  
  
  observeEvent(input$freeze_img, {
    snap <- last_rendered_pixels()
    if (is.null(snap) || is.null(snap$img)) {
      showNotification("Nothing to freeze yet \u2014 select and display a channel first.",
                       type = "warning", duration = 4)
      return()
    }
    
    src_name <- if (!is.null(snap$source)) snap$source else "channel"
    
    base <- sub("^SAVED_", "", src_name)
    base <- gsub("[^A-Za-z0-9_.-]", "_", base)
    frozen_name <- paste0("SAVED_", base)
    
    cur <- imgs(); if (is.null(cur)) cur <- list()
    
    if (frozen_name %in% names(cur)) {
      i <- 2
      while (paste0(frozen_name, "_", i) %in% names(cur)) i <- i + 1
      frozen_name <- paste0(frozen_name, "_", i)
    }
    
    img_dir <- file.path(tempdir(), paste0("served_images_", session$token))
    outpath <- file.path(img_dir, paste0(frozen_name, ".png"))
    fimg <- snap$img
    if (length(dim(fimg)) == 2) fimg <- EBImage::rgbImage(fimg, fimg, fimg)
    EBImage::writeImage(fimg, outpath)
    
    cur[[frozen_name]] <- outpath
    imgs(cur)
    
    
    lgs <- img_legends()
    if (src_name %in% names(lgs)) { lgs[[frozen_name]] <- lgs[[src_name]]; img_legends(lgs) }
    
    showNotification(paste0("Frozen appearance saved as \u201C", frozen_name, "\u201D."),
                     type = "message", duration = 4)
    session$sendCustomMessage("openPanel",     list(id = "fp-channels"))
    session$sendCustomMessage("selectChannel", list(key = frozen_name))
  })
  
  
  selected_path <- reactive({
    req(imgs()); req(input$selected_img)
    imgs()[[input$selected_img]]
  })
  
  output$main_image <- renderImage({
    req(selected_path())
    files <- isolate(imgs())         
    fname <- files[[input$selected_img]]
    
    if(is.character(fname)){
      fname <- EBImage::readImage(fname)
    }
    
    
    if(all(dim(fname)[1:2] < 800)){
      d  <- dim(fname)
      ar <- d[1] / d[2]                       
      target_long <- 2000
      if (ar >= 1) {
        new_w <- target_long
        new_h <- round(target_long / ar)
      } else {
        new_h <- target_long
        new_w <- round(target_long * ar)
      }
      fname <- EBImage::resize(fname, w = new_w, h = new_h)
    }
    
    is_composite <- grepl("composite", input$selected_img, ignore.case = TRUE)
    is_overlay   <- grepl("overlay", input$selected_img, ignore.case = TRUE)
    is_agreement <- grepl("agreement", input$selected_img, ignore.case = TRUE)
    is_interface <- grepl("interface", input$selected_img, ignore.case = TRUE)
    is_frozen    <- grepl("^SAVED_", input$selected_img)
    is_prior <- input$selected_img %in% c("H&E Image","Entropy","Non-fibrillar interstitial tissue",
                                          "Fibrillar interstitial tissue","Nuclei and cells")
    
    if (is_composite || is_agreement || is_interface) {
      
      sources <- composite_sources()[[input$selected_img]]
      if (!is.null(sources) && length(sources) == 2 &&
          all(sources %in% names(imgs()))) {
        
        read_layer <- function(nm) {
          p <- imgs()[[nm]]
          img <- if (is.character(p)) EBImage::readImage(p) else p
          if (all(dim(img)[1:2] < 800)) {
            d  <- dim(img); ar <- d[1] / d[2]
            target_long <- 2000
            new_w <- if (ar >= 1) target_long else round(target_long * ar)
            new_h <- if (ar >= 1) round(target_long / ar) else target_long
            img  <- EBImage::resize(img, w = new_w, h = new_h)
          }
          if (length(dim(img)) == 2) img <- EBImage::rgbImage(img, img, img)
          img
        }
        
        img1 <- tryCatch(read_layer(sources[1]), error = function(e) NULL)
        img2 <- tryCatch(read_layer(sources[2]), error = function(e) NULL)
        
        if (!is.null(img1) && !is.null(img2)) {
          
          if (!all(dim(img1)[1:2] == dim(img2)[1:2]))
            img2 <- EBImage::resize(img2, w=dim(img1)[1], h=dim(img1)[2])
          
          thr1 <- if (!is.null(comp_thr1_t()))      comp_thr1_t()      else 0.30
          thr2 <- if (!is.null(comp_thr2_t()))      comp_thr2_t()      else 0.30
          ct1  <- if (!is.null(comp_contrast1_t())) comp_contrast1_t() else 0
          ct2  <- if (!is.null(comp_contrast2_t())) comp_contrast2_t() else 0
          
          img1[img1 < thr1] <- 0
          img2[img2 < thr2] <- 0
          img1 <- img1 ^ exp(-1.9 * ct1)
          img2 <- img2 ^ exp(-1.9 * ct2)
          
          if (is_agreement) {
            
            gray1 <- EBImage::channel(img1, "gray")
            gray2 <- EBImage::channel(img2, "gray")
            ag <- 1 - abs(gray1 - gray2)
            mask1 <- gray1 > 0
            mask2 <- gray2 > 0
            ag <- ag * mask1 * mask2
            fname <- EBImage::normalize(ag)
          } else if (is_interface) {
            
            saved <- imgs()[[input$selected_img]]
            if (!is.null(saved) && is.character(saved) && file.exists(saved)) {
              fname <- EBImage::readImage(saved)
              
              if (all(dim(fname)[1:2] < 800)) {
                d <- dim(fname); ar <- d[1] / d[2]
                target_long <- 2000
                new_w <- if (ar >= 1) target_long else round(target_long * ar)
                new_h <- if (ar >= 1) round(target_long / ar) else target_long
                fname <- EBImage::resize(fname, w = new_w, h = new_h)
              }
            } else {
              
              fname <- (img1 * .6 + img2 * .4)
              fname <- EBImage::normalize(fname) ^ .4
            }
          } else {
            
            if (any(grepl("composite", c(sources[1], sources[2])))) {
              fname <- ((img1 + img2) / 2) ^ .4
            } else {
              fname <- (img1 * .6 + img2 * .4)
              fname <- EBImage::normalize(fname) ^ .4
            }
          }
        }
        
      }
      
    } else if (is_overlay || is_frozen) {
      fname <- fname   
    } else {
      
      thr <- if (!is.null(intensity_thr_t())) intensity_thr_t() else 0.30
      fname[fname < thr] <- 0
      cp <- crispr()
      if(cp==1){
        fname <- fname ^.4
        fname <- (fname - .5) * 2 + .5
        fname <- pmin(pmax(fname,0),1)
        crispr(0)
      }
      l <- thr
      k <- isolate(trl())
      if(input$selected_img %in% names(k)){
        k[[input$selected_img]] <- l; trl(k)
      } else {
        l <- list(l); names(l) <- input$selected_img; trl(c(k, isolate(l)))
      }
      
    }
    
    if (!is_composite && !is_interface && !is_frozen) {
      ct <- if (!is.null(contrastup_t())) contrastup_t() else 0
      fname <- fname ^ exp(-1.9 * ct)
    }
    
    
    last_rendered_pixels(list(img = fname, source = input$selected_img))
    
    tmpfile <- tempfile(fileext = ".png")
    EBImage::writeImage(fname, tmpfile)
    
    
    if(!file.exists(tmpfile)){
      stop("Temporary file was not created!")
    }
    
    list(
      src = tmpfile,
      contentType = "image/png",
      alt = input$selected_img
    )
  }, deleteFile = TRUE)
  
  
  observe({
    req(imgs())
    
    prior_names_ui <- c("H&E Image","Entropy","Non-fibrillar interstitial tissue",
                        "Fibrillar interstitial tissue","Nuclei and cells")
    mergeable <- setdiff(names(imgs()), prior_names_ui)
    updateSelectInput(session, "img1", choices = mergeable)
    updateSelectInput(session, "img2", choices = mergeable)
    updateSelectInput(session, "he_overlay_gene", choices = mergeable)
    updateSelectInput(session, "imgcust", choices = names(imgs()))
  })
  
  
  composite_sources <- reactiveVal(list())
  interface_last_stats <- reactiveVal(NULL)
  
  output$interface_readout_ui <- renderUI({
    s <- interface_last_stats()
    if (is.null(s)) {
      return(tags$div(
        style="font-size:11px; color:#9aa0b8; font-style:italic;",
        "Run a merge to populate the interface statistics."
      ))
    }
    
    
    a_pct <- round(100 * s$A_side_fraction, 1)
    b_pct <- round(100 * s$B_side_fraction, 1)
    asym_text <- if (abs(s$A_side_fraction - 0.5) < 0.05)
      "balanced"
    else if (s$A_side_fraction > 0.5)
      sprintf("%s-leaning", s$channel_A)
    else
      sprintf("%s-leaning", s$channel_B)
    
    row_div <- function(label, value, tone = "default") {
      vcol <- switch(tone,
                     "good"    = "#1f8a3c",
                     "warning" = "#a86a00",
                     "muted"   = "#7880a0",
                     "#3a4060")
      tags$div(
        style="display:flex; justify-content:space-between; font-size:11px;
               padding:2px 0; border-bottom:1px dotted #e0e3ea;",
        tags$span(style="color:#5c6480;", label),
        tags$span(style=paste0("color:", vcol, "; font-family:Consolas,monospace;"),
                  value)
      )
    }
    
    tagList(
      row_div("Channel A",           s$channel_A,                                  "muted"),
      row_div("Channel B",           s$channel_B,                                  "muted"),
      row_div("Reach",               sprintf("%.1f \u00b5m", s$reach_um)),
      row_div("Threshold mode",      s$threshold_mode,                             "muted"),
      row_div("Threshold A",         sprintf("%.3f", s$threshold_A)),
      row_div("Threshold B",         sprintf("%.3f", s$threshold_B)),
      
      tags$div(style="margin-top:6px; font-size:10px; color:#5c6480;
                     font-weight:600; letter-spacing:0.05em;
                     border-bottom:1px solid #c8d5ee;", "INTERFACE (adjacency, A vs B)"),
      row_div("Interface area",      sprintf("%s \u00b5m\u00b2 (%.2f%%)",
                                             format(s$interface_area_um2,
                                                    big.mark=",",
                                                    scientific=FALSE,
                                                    trim=TRUE),
                                             s$interface_area_pct),
              if (s$interface_area_pct < 0.1) "warning" else "good"),
      row_div("A-side share",        sprintf("%.1f%% (%s)", a_pct, asym_text)),
      row_div("B-side share",        sprintf("%.1f%%", b_pct)),
      row_div("Mean strength",       sprintf("%.4f", s$mean_strength)),
      
      tags$div(style="margin-top:6px; font-size:10px; color:#5c6480;
                     font-weight:600; letter-spacing:0.05em;
                     border-bottom:1px solid #c8d5ee;", "MERGE / OVERLAP (A \u2229 B)"),
      row_div("Overlap area",        sprintf("%s \u00b5m\u00b2",
                                             format(s$overlap_area_um2,
                                                    big.mark=",",
                                                    scientific=FALSE,
                                                    trim=TRUE)),
              if (s$overlap_area_um2 < 1) "warning" else "good"),
      row_div("Overlap / Union",     sprintf("%.1f%% (Jaccard %.3f)",
                                             s$overlap_pct_of_union, s$jaccard_AB)),
      row_div("Union area (A \u222a B)",
              sprintf("%s \u00b5m\u00b2",
                      format(s$union_area_um2,
                             big.mark=",",
                             scientific=FALSE,
                             trim=TRUE)), "muted"),
      row_div("Only A",              sprintf("%s \u00b5m\u00b2",
                                             format(s$onlyA_area_um2,
                                                    big.mark=",",
                                                    scientific=FALSE,
                                                    trim=TRUE)), "muted"),
      row_div("Only B",              sprintf("%s \u00b5m\u00b2",
                                             format(s$onlyB_area_um2,
                                                    big.mark=",",
                                                    scientific=FALSE,
                                                    trim=TRUE)), "muted")
    )
  })
  
  output$layer_controls_ui <- renderUI({
    sel     <- input$selected_img
    
    is_comp <- !is.null(sel) &&
      grepl("^(composite|agreement|interface)_", sel, ignore.case=TRUE)
    is_meta <- !is.null(sel) && grepl("^META_", sel)
    sources <- if (is_comp) composite_sources()[[sel]] else NULL
    
    if (is_meta) {
      
      div(
        div(class="fp-section", "Metadata image"),
        tags$p(style="font-size:11px;color:#7880a0;margin:0 0 6px 0;",
               "Spot size controls the dilation kernel radius applied to each
                metadata point before rainbow/colour mapping."),
        sliderInput("meta_spot_size", "Spot size (pixels)",
                    min=3, max=41, value=11, step=2, width="100%"),
        actionButton("apply_meta_spot", "🔁  Apply new spot size",
                     class="btn-fp-secondary",
                     style="width:100%; margin-top:4px;
                            background:linear-gradient(135deg,#3a006b,#7a20c0) !important;
                            border-color:#7a20c0 !important;
                            color:#ffffff !important;")
      )
    } else if (is_comp && !is.null(sources) && length(sources) == 2) {
      
      lbl <- function(s) if (nchar(s) > 22) paste0(substr(s,1,20),"…") else s
      
      tagList(
        
        div(style="background:#eef3ff;border-radius:6px;padding:8px 10px;margin-bottom:6px;
                   border-left:3px solid #2c5be0;",
            div(class="fp-section",
                style="color:#2c5be0; border-color:#ccd6f8;",
                paste0("Layer A — ", lbl(sources[1]))),
            sliderInput("comp_thr1", "Noise gate",
                        min=0, max=1, value=0.30, step=0.01, width="100%"),
            sliderInput("comp_contrast1", "Contrast",
                        min=0, max=1, value=0, step=0.01, width="100%")
        ),
        
        div(style="background:#edfff4;border-radius:6px;padding:8px 10px;
                   border-left:3px solid #1e7a45;",
            div(class="fp-section",
                style="color:#1e7a45; border-color:#b8e8cb;",
                paste0("Layer B — ", lbl(sources[2]))),
            sliderInput("comp_thr2", "Noise gate",
                        min=0, max=1, value=0.30, step=0.01, width="100%"),
            sliderInput("comp_contrast2", "Contrast",
                        min=0, max=1, value=0, step=0.01, width="100%")
        )
      )
      
    } else {
      
      div(
        div(class="fp-section", "Image adjustments"),
        sliderInput("intensity_threshold", "Noise gate",
                    min=0, max=1, value=0.30, step=0.01, width="100%"),
        sliderInput("contrastup", "Contrast boost",
                    min=0, max=1, value=0, step=0.01, width="100%")
      )
    }
  })
  outputOptions(output, "layer_controls_ui", suspendWhenHidden = FALSE)
  
  observeEvent(input$do_blend, {
    req(input$img1, input$img2, imgs())
    
    withProgress(message = "Processing images...", value = 0, {
      
      f1 <- imgs()[[input$img1]]
      f2 <- imgs()[[input$img2]]
      incProgress(0.1, detail = "Reading images")
      
      img1 <- EBImage::readImage(f1)
      img2 <- EBImage::readImage(f2)
      incProgress(0.2, detail = "Resizing images")
      
      if (!all(dim(img1)[1:2] == dim(img2)[1:2])) {
        
        tw <- max(dim(img1)[1], dim(img2)[1])
        th <- max(dim(img1)[2], dim(img2)[2])
        img1 <- EBImage::resize(img1, w = tw, h = th)
        img2 <- EBImage::resize(img2, w = tw, h = th)
      }
      
      incProgress(0.1, detail = "Converting to RGB")
      if (length(dim(img1)) == 2) img1 <- EBImage::rgbImage(img1, img1, img1)
      if (length(dim(img2)) == 2) img2 <- EBImage::rgbImage(img2, img2, img2)
      if (dim(img1)[3] != dim(img2)[3]) stop("Images have incompatible channels.")
      
      if(any(c(input$img1,input$img2)=="H&E Image")){
        if(min(img1)!=0 & input$img1 == "H&E Image"){
          bg_cutoff <- min(img1) + 0.01
          img2[img2 < bg_cutoff] <- 0
        }else{
          bg_cutoff <- min(img2) + 0.01 
          img1[img1 < bg_cutoff] <- 0
        }
      }
      
      thr_merge <- if (!is.null(input$intensity_threshold)) input$intensity_threshold else 0.30
      
      if(any(c(input$img1,input$img2)=="H&E Image")){
        if(input$img1 %in% "H&E Image"){
          img2[img2 < thr_merge] <- 0
        }else{
          img1[img1 < thr_merge] <- 0
        }
      }else{
        img1[img1 < thr_merge] <- 0
        img2[img2 < thr_merge] <- 0
      }
      
      new_imgs <- list()
      
      
      {
        incProgress(0.2, detail = "Blending images (composite)")
        
        
        if (any(grepl("composite", c(input$img1, input$img2)))) {
          if (grepl("composite", input$img1)) {
            n <- 2
            blended <- (img1 * (n-1) + img2) / n
            blended <- blended ^.4
          } else {
            n <- 2
            blended <- (img2 * (n-1) + img1) / n
            blended <- blended ^.4
          }
        }else{
          blended <- img1*.6 + (img2*.4)
          blended <- EBImage::normalize(blended)
          blended <- blended ^.4
        }
        
        
        
        incProgress(0.5, detail = "Computing agreement map")
        if(input$img1 == input$img2){
          agreement <- img1
          if(is.na(dim(img1)[3])){
            agreement@.Data[agreement@.Data>0] <- 1
          } else {
            for(i in 1:dim(img1)[3]){
              k <- img1[,,i]
              k@.Data[k@.Data>0] <- 1
              agreement[,,i] <- k
            }
          }
        } else {
          gray1 <- EBImage::channel(img1, "gray")
          gray2 <- EBImage::channel(img2, "gray")
          agreement <- 1 - abs(gray1 - gray2)
          mask1 <- gray1 > 0
          mask2 <- gray2 > 0
          agreement <- agreement * mask1 * mask2
          agreement <- normalize(agreement)
        }
        
        
        incProgress(0.75, detail = "Saving images")
        
        v <- names(imgs())
        v <- gsub(".png","",v)
        outname <- paste0("composite_", input$img1, "_", input$img2)
        
        if(length(v[grepl(outname,v,ignore.case = T)])){
          l <- length(v[grepl(outname,v,ignore.case = T)])+1
          outname <- paste0(outname,"_",l)
        }else{
          outname <- outname
        }
        
        outpath <- file.path(img_dir, paste0(outname,".png"))
        EBImage::writeImage(blended, outpath)
        new_imgs[[outname]] <- outpath
        
        cs <- composite_sources()
        cs[[outname]] <- c(input$img1, input$img2)
        composite_sources(cs)
        
        outname2 <- paste0("agreement_", input$img1, "_", input$img2)
        
        if(length(v[grepl(outname2,v,ignore.case = T)])>0){
          l <- length(v[grepl(outname2,v,ignore.case = T)])+1
          outname2 <- paste0(outname2,"_",l)
        }else{
          outname2 <- outname2
        }
        
        outpath2 <- file.path(img_dir, paste0(outname2,".png"))
        EBImage::writeImage(agreement, outpath2)
        new_imgs[[outname2]] <- outpath2
        
        cs <- composite_sources()
        cs[[outname2]] <- c(input$img1, input$img2)
        composite_sources(cs)
        
        
        if (input$img1 != input$img2) {
          incProgress(0.85, detail = "Computing interface map")
          
          
          mpp <- tryCatch(orig()$microns_per_pixel, error = function(e) NULL)
          if (is.null(mpp) || !is.finite(mpp) || mpp <= 0) mpp <- 1.0
          
          
          reach_um <- if (!is.null(input$interface_reach_um))
            as.numeric(input$interface_reach_um) else 15
          if (!is.finite(reach_um) || reach_um <= 0) reach_um <- 15
          
          
          sigma_px <- max(0.6, reach_um / (2 * mpp))
          
          
          
          
          gA <- EBImage::channel(img1, "gray")
          gB <- EBImage::channel(img2, "gray")
          gA_arr <- EBImage::imageData(gA)
          gB_arr <- EBImage::imageData(gB)
          
          
          .scale_unit <- function(x) {
            x[!is.finite(x)] <- 0
            m <- max(x, na.rm = TRUE)
            if (m > 0) x / m else x
          }
          gA_arr <- .scale_unit(gA_arr)
          gB_arr <- .scale_unit(gB_arr)
          
          
          thr_mode  <- if (!is.null(input$interface_thr_mode))
            input$interface_thr_mode else "adaptive"
          thr_fixed <- if (!is.null(input$interface_thr_fixed))
            as.numeric(input$interface_thr_fixed) else 0.10
          
          .pick_thr <- function(arr) {
            if (thr_mode == "fixed") return(thr_fixed)
            
            v <- arr[is.finite(arr) & arr > 0.001]
            if (length(v) < 50) return(0.40)
            otsu_th <- tryCatch(
              EBImage::otsu(EBImage::Image(arr, colormode = "Grayscale")),
              error = function(e) NA_real_
            )
            
            q_lo <- as.numeric(stats::quantile(v, 0.50))
            q_hi <- as.numeric(stats::quantile(v, 0.85))
            if (is.na(otsu_th)) return(q_lo)
            pmin(pmax(otsu_th, q_lo), q_hi)
          }
          thrA <- .pick_thr(gA_arr)
          thrB <- .pick_thr(gB_arr)
          
          
          
          .sharpen_pres <- function(x, thr, k = 24) {
            1 / (1 + exp(-k * (x - thr)))
          }
          A_strong <- .sharpen_pres(gA_arr, thrA)
          B_strong <- .sharpen_pres(gB_arr, thrB)
          
          
          
          margin <- 0.20
          contrast <- A_strong - B_strong
          A_terr <- A_strong * pmax(0, pmin(1, (contrast - margin) / margin))
          B_terr <- B_strong * pmax(0, pmin(1, (-contrast - margin) / margin))
          
          
          
          
          .sobel <- function(arr) {
            kx <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), 3, 3)
            ky <- t(kx)
            gx <- EBImage::imageData(
              EBImage::filter2(EBImage::Image(arr, colormode = "Grayscale"), kx))
            gy <- EBImage::imageData(
              EBImage::filter2(EBImage::Image(arr, colormode = "Grayscale"), ky))
            sqrt(gx * gx + gy * gy)
          }
          
          
          A_terr_smooth <- EBImage::imageData(
            EBImage::gblur(EBImage::Image(A_terr, colormode = "Grayscale"),
                           sigma = max(1, sigma_px / 2)))
          B_terr_smooth <- EBImage::imageData(
            EBImage::gblur(EBImage::Image(B_terr, colormode = "Grayscale"),
                           sigma = max(1, sigma_px / 2)))
          
          edge_A <- .sobel(A_terr_smooth)
          edge_B <- .sobel(B_terr_smooth)
          
          
          ge <- max(c(edge_A, edge_B), na.rm = TRUE)
          if (!is.finite(ge) || ge <= 0) ge <- 1
          edge_A <- edge_A / ge
          edge_B <- edge_B / ge
          
          
          
          reach_sigma <- max(sigma_px * 2, 6)
          B_terr_field <- EBImage::imageData(
            EBImage::gblur(EBImage::Image(B_terr, colormode = "Grayscale"),
                           sigma = reach_sigma))
          A_terr_field <- EBImage::imageData(
            EBImage::gblur(EBImage::Image(A_terr, colormode = "Grayscale"),
                           sigma = reach_sigma))
          
          
          gnorm <- max(c(B_terr_field, A_terr_field), na.rm = TRUE)
          if (!is.finite(gnorm) || gnorm <= 0) gnorm <- 1
          B_terr_field <- B_terr_field / gnorm
          A_terr_field <- A_terr_field / gnorm
          
          
          face_A <- edge_A * B_terr_field
          face_B <- edge_B * A_terr_field
          
          
          maxA <- max(face_A, na.rm = TRUE); if (maxA <= 0) maxA <- 1
          maxB <- max(face_B, na.rm = TRUE); if (maxB <= 0) maxB <- 1
          face_A <- face_A / maxA
          face_B <- face_B / maxB
          
          
          face_A <- face_A ^ 0.45
          face_B <- face_B ^ 0.45
          
          
          face_A <- pmin(face_A * 1.6, 1)
          face_B <- pmin(face_B * 1.6, 1)
          
          
          colA <- if (!is.null(img_legends()[[input$img1]]))
            img_legends()[[input$img1]]$col else "#FF00FF"
          colB <- if (!is.null(img_legends()[[input$img2]]))
            img_legends()[[input$img2]]$col else "#00FFFF"
          .hex2rgb <- function(hex) {
            hex <- gsub("^#", "", hex)
            if (nchar(hex) == 3) hex <- paste0(rep(strsplit(hex,"")[[1]], each=2), collapse="")
            c(strtoi(substr(hex,1,2),16L),
              strtoi(substr(hex,3,4),16L),
              strtoi(substr(hex,5,6),16L)) / 255
          }
          cA <- tryCatch(.hex2rgb(colA), error = function(e) c(1, 0, 1))
          cB <- tryCatch(.hex2rgb(colB), error = function(e) c(0, 1, 1))
          
          h <- nrow(gA_arr); w <- ncol(gA_arr)
          interface_img <- array(0, dim = c(h, w, 3))
          for (k in 1:3) {
            interface_img[,,k] <- pmin(face_A * cA[k] + face_B * cB[k], 1)
          }
          interface_rendered <- EBImage::Image(interface_img, colormode = "Color")
          
          
          support_thr <- 0.05
          support_mask <- (face_A + face_B) > support_thr
          n_support_px <- sum(support_mask, na.rm = TRUE)
          total_px <- length(face_A)
          area_um2 <- n_support_px * mpp * mpp
          area_pct <- 100 * n_support_px / max(total_px, 1)
          
          sumA_in  <- sum(face_A[support_mask], na.rm = TRUE)
          sumB_in  <- sum(face_B[support_mask], na.rm = TRUE)
          total_in <- sumA_in + sumB_in
          A_frac   <- if (total_in > 0) sumA_in / total_in else 0.5
          
          mean_strength <- if (n_support_px > 0)
            mean((face_A + face_B)[support_mask], na.rm = TRUE) else 0
          
          
          pres_thr <- 0.5  
          present_A <- A_strong > pres_thr
          present_B <- B_strong > pres_thr
          overlap_mask <- present_A & present_B
          union_mask   <- present_A | present_B
          onlyA_mask   <- present_A & !present_B
          onlyB_mask   <- present_B & !present_A
          
          n_overlap_px <- sum(overlap_mask, na.rm = TRUE)
          n_union_px   <- sum(union_mask,   na.rm = TRUE)
          n_onlyA_px   <- sum(onlyA_mask,   na.rm = TRUE)
          n_onlyB_px   <- sum(onlyB_mask,   na.rm = TRUE)
          
          overlap_area_um2 <- n_overlap_px * mpp * mpp
          union_area_um2   <- n_union_px   * mpp * mpp
          overlap_pct_of_union <- if (n_union_px > 0)
            100 * n_overlap_px / n_union_px else 0
          jaccard <- if (n_union_px > 0) n_overlap_px / n_union_px else 0
          
          
          outname3 <- paste0("interface_", input$img1, "_", input$img2)
          if (length(v[grepl(outname3, v, ignore.case = TRUE)]) > 0) {
            l <- length(v[grepl(outname3, v, ignore.case = TRUE)]) + 1
            outname3 <- paste0(outname3, "_", l)
          }
          outpath3 <- file.path(img_dir, paste0(outname3, ".png"))
          EBImage::writeImage(interface_rendered, outpath3)
          new_imgs[[outname3]] <- outpath3
          
          cs <- composite_sources()
          cs[[outname3]] <- c(input$img1, input$img2)
          composite_sources(cs)
          
          
          interface_stats <- list(
            channel_A           = input$img1,
            channel_B           = input$img2,
            reach_um            = round(reach_um, 1),
            sigma_px            = round(sigma_px, 2),
            mpp_um_per_px       = round(mpp, 3),
            threshold_mode      = thr_mode,
            threshold_A         = round(thrA, 3),
            threshold_B         = round(thrB, 3),
            interface_area_um2  = round(area_um2, 0),
            interface_area_pct  = round(area_pct, 2),
            overlap_area_um2    = round(overlap_area_um2, 0),
            overlap_pct_of_union= round(overlap_pct_of_union, 2),
            union_area_um2      = round(union_area_um2, 0),
            jaccard_AB          = round(jaccard, 3),
            onlyA_area_um2      = round(n_onlyA_px * mpp * mpp, 0),
            onlyB_area_um2      = round(n_onlyB_px * mpp * mpp, 0),
            A_side_fraction     = round(A_frac, 3),
            B_side_fraction     = round(1 - A_frac, 3),
            mean_strength       = round(mean_strength, 4)
          )
          vs <- isolate(vals_store$results)
          vs[[outname3]] <- interface_stats
          vals_store$results <- vs
          interface_last_stats(interface_stats)
        }
        
        
      }
      
      current <- imgs()
      if (is.null(current)) current <- list()
      imgs(c(current, new_imgs))
      
      session$sendCustomMessage("openPanel",     list(id  = "fp-channels"))
      session$sendCustomMessage("selectChannel", list(key = names(new_imgs)[1]))
      
    })
  })
  
  .blend_he_gene <- function(he_arr, gn_arr, alpha_he, brightness) {
    
    .ensure_rgb <- function(a) {
      if (is.null(a)) return(NULL)
      d <- dim(a)
      if (length(d) == 2) {
        return(array(rep(as.numeric(a), 3), dim = c(d, 3)))
      }
      if (length(d) == 3 && d[3] >= 3) {
        return(a[,,1:3, drop = FALSE])
      }
      if (length(d) == 3 && d[3] == 1) {
        slice <- a[,,1, drop = TRUE]
        return(array(rep(as.numeric(slice), 3), dim = c(dim(slice), 3)))
      }
      
      v <- as.numeric(a)
      side <- floor(sqrt(length(v) / 3))
      array(v[1:(side*side*3)], dim = c(side, side, 3))
    }
    
    he_arr <- .ensure_rgb(he_arr)
    gn_arr <- .ensure_rgb(gn_arr)
    
    he_w <- dim(he_arr)[1]; he_h <- dim(he_arr)[2]
    if (dim(gn_arr)[1] != he_w || dim(gn_arr)[2] != he_h) {
      gn_img  <- EBImage::Image(gn_arr, colormode = "Color")
      gn_img  <- EBImage::resize(gn_img, w = he_w, h = he_h)
      gn_arr  <- imageData(gn_img)
      gn_arr  <- .ensure_rgb(gn_arr)
    }
    
    if (!identical(dim(he_arr), dim(gn_arr))) {
      warning(sprintf(
        "blend_he_gene: dimension mismatch after normalization - he=%s, gene=%s; returning H&E unchanged",
        paste(dim(he_arr), collapse="x"),
        paste(dim(gn_arr), collapse="x")
      ), call. = FALSE)
      return(he_arr)
    }
    
    he_dim      <- he_arr * alpha_he
    gn_boost    <- pmin(gn_arr * brightness, 1)
    blended_arr <- pmin(he_dim + gn_boost * (1 - he_dim * 0.5), 1)
    pmax(blended_arr, 0)
  }
  
  .img_to_arr3 <- function(img) {
    a <- imageData(img)
    if (length(dim(a)) == 2) array(rep(a, 3), dim = c(dim(a), 3))
    else a[,,1:3, drop=FALSE]
  }
  
  he_overlay_sources <- reactive({
    req(input$he_overlay_gene, imgs())
    he_path <- imgs()[["H&E Image"]]
    gn_path <- imgs()[[input$he_overlay_gene]]
    req(!is.null(he_path), !is.null(gn_path))
    
    he_img <- EBImage::readImage(he_path)
    gn_img <- EBImage::readImage(gn_path)
    list(he_full = .img_to_arr3(he_img),
         gn_full = .img_to_arr3(gn_img))
  })
  
  output$he_overlay_preview <- renderImage({
    srcs <- he_overlay_sources()
    alpha_he   <- if (!is.null(he_opacity_d()))      he_opacity_d()      else 0.5
    brightness <- if (!is.null(gene_brightness_d())) gene_brightness_d() else 1.5
    
    he_arr <- srcs$he_full; gn_arr <- srcs$gn_full
    max_px <- 400L
    sc <- max_px / max(dim(he_arr)[1], dim(he_arr)[2])
    if (sc < 1) {
      he_img <- EBImage::Image(he_arr, colormode = "Color")
      he_img <- EBImage::resize(he_img,
                                w = round(dim(he_arr)[1]*sc),
                                h = round(dim(he_arr)[2]*sc))
      he_arr <- imageData(he_img)[,,1:3, drop=FALSE]
      
      gn_img <- EBImage::Image(gn_arr, colormode = "Color")
      gn_img <- EBImage::resize(gn_img,
                                w = round(dim(gn_arr)[1]*sc),
                                h = round(dim(gn_arr)[2]*sc))
      gn_arr <- imageData(gn_img)[,,1:3, drop=FALSE]
    }
    
    blended_arr <- .blend_he_gene(he_arr, gn_arr, alpha_he, brightness)
    blended <- EBImage::Image(blended_arr, colormode = "Color")
    
    tmp <- tempfile(fileext = ".png")
    EBImage::writeImage(blended, tmp)
    list(src = tmp, contentType = "image/png",
         style = "max-width:100%; max-height:100%; display:block;
                  margin:0 auto; border-radius:4px; object-fit:contain;")
  }, deleteFile = TRUE)
  outputOptions(output, "he_overlay_preview", suspendWhenHidden = FALSE)
  outputOptions(output, "channel_rail_ui",   suspendWhenHidden = FALSE)
  
  observeEvent(input$do_he_overlay, {
    req(input$he_overlay_gene, imgs())
    gene_nm <- input$he_overlay_gene
    srcs <- he_overlay_sources()
    
    withProgress(message = "Acquiring H&E overlay at full resolution…", value = 0, {
      incProgress(0.3, detail = "Blending")
      alpha_he   <- input$he_opacity
      brightness <- input$gene_brightness
      blended_arr <- .blend_he_gene(srcs$he_full, srcs$gn_full, alpha_he, brightness)
      blended <- EBImage::Image(blended_arr, colormode = "Color")
      
      incProgress(0.8, detail = "Saving")
      v <- names(imgs()); v <- gsub(".png","",v)
      outname <- paste0("HEoverlay_", gene_nm)
      if (any(grepl(outname, v, ignore.case=TRUE))) {
        l <- sum(grepl(outname, v, ignore.case=TRUE)) + 1
        outname <- paste0(outname, "_", l)
      }
      outpath <- file.path(img_dir, paste0(outname, ".png"))
      EBImage::writeImage(blended, outpath)
      
      current <- imgs(); if (is.null(current)) current <- list()
      imgs(c(current, setNames(list(outpath), outname)))
      
      session$sendCustomMessage("openPanel",     list(id  = "fp-channels"))
      session$sendCustomMessage("selectChannel", list(key = outname))
    })
  })
  
  
  
  
  
  calc_image <- reactiveVal(NULL)
  
  
  observeEvent(input$calcs, {
    calc_image(input$selected_img)
  })
  
  
  observeEvent(input$selected_img, {
    if (!is.null(calc_image()) && input$selected_img != calc_image()) {
      calc_image(NULL)  
    }
    if (!is.null(lgnd_image()) && input$selected_img != lgnd_image()) {
      lgnd_image(NULL)
    }
    
    
    lgs <- img_legends()
    if (!is.null(input$selected_img) && input$selected_img %in% names(lgs)) {
      lgnd(lgs[[ input$selected_img ]])   
    } else {
      lgnd(NULL)
    }
  })
  
  
  output$autocorr_panel_ui <- renderUI({
    req(calc_image())
    req(identical(input$selected_img, calc_image()))
    
    tags$div(
      id = "autocorr_fp",
      class = "fp",
      style = "position:absolute; top:80px; left:20px; width:340px; z-index:9999;",
      tags$div(class = "fp-header",
               tags$span("📈", style="font-size:13px;"),
               tags$span("Spatial Autocorrelation", class="fp-title"),
               tags$button(class="fp-min",
                           onclick="togglePanelMinimize('autocorr_fp')",
                           "−"),
               tags$button(class="fp-close",
                           onclick="document.getElementById('autocorr_fp').style.display='none'",
                           "✕")
      ),
      tags$div(class="fp-body",
               plotOutput("autocorr_plot",  height="180px"),
               tableOutput("autocorr_table")
      )
    )
  })
  
  
  autocorr_cache <- reactiveVal(NULL) 
  
  observeEvent(input$calcs, {
    req(imgs())
    req(input$selected_img)
    
    if (input$selected_img != "") {
      
      
      
      
      withProgress(message = "Processing images...", value = 0, {  
        
        
        files <- imgs()
        fname <- files[[input$selected_img]]
        incProgress(.25,"image reformatting...")
        
        
        if(is.character(fname)){
          img <- EBImage::readImage(fname)
        }
        incProgress(.5,"image analysis...")
        
        
        df <- multi_scale_autocorr(img,micron_per_pixel = orig()$microns_per_pixel)
        incProgress(1,"autocorrelation calculated...")
        df
      })
      
      
      autocorr_cache(df)
    }
  })
  
  
  
  
  output$autocorr_table <- renderTable({
    autocorr_cache()
  }, digits = 3)
  
  
  output$autocorr_plot <- renderPlot({
    req(autocorr_cache())
    p <- plot_spatial_autocorrelation_profile(autocorr_cache())
    p + theme(
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  }, bg = "transparent")
  
  
  lgnd_image <- reactiveVal(NULL)
  
  
  observeEvent(input$lgnd, {
    lgnd_image(input$selected_img)
  })
  
  
  observeEvent(input$selected_img, {
    if (!is.null(lgnd_image()) && input$selected_img != lgnd_image()) {
      lgnd_image(NULL)  
    }
  })
  
  
  observeEvent(input$selected_img, {
    is_composite <- grepl("composite", input$selected_img, ignore.case = TRUE)
    is_overlay   <- grepl("overlay", input$selected_img, ignore.case = TRUE)
    is_agreement   <- grepl("agreement", input$selected_img, ignore.case = TRUE)
    is_interface <- grepl("interface", input$selected_img, ignore.case = TRUE)
    
    if (is_composite || is_overlay || is_agreement || is_interface) {
      shinyjs::disable("lgnd")
      
      if (is_interface) {
        shinyjs::enable("allvals")
      } else {
        shinyjs::disable("allvals")
      }
    } else {
      shinyjs::enable("lgnd")
      if(input$selected_img%in%c("H&E Image","Entropy", "Non-fibrillar interstitial tissue",
                                 "Fibrillar interstitial tissue", "Nuclei and cells", "Stiff and soft tissue")){
        shinyjs::disable("allvals")
      }else{
        shinyjs::enable("allvals")
      }
    }
  })
  
  
  output$legend_overlay <- renderUI({
    req(!is.null(lgnd()))
    req(!is.null(lgnd_image()))
    req(identical(input$selected_img, lgnd_image()))
    
    tags$div(
      id = "legend_fp",
      class = "fp",
      style = "position:absolute; top:80px; right:20px; width:280px; z-index:9999;",
      tags$div(class="fp-header",
               tags$span("🎨", style="font-size:13px;"),
               tags$span("Colour Legend", class="fp-title"),
               tags$button(class="fp-min",
                           onclick="togglePanelMinimize('legend_fp')",
                           "−"),
               tags$button(class="fp-close",
                           onclick="document.getElementById('legend_fp').style.display='none'",
                           "✕")
      ),
      tags$div(class="fp-body",
               plotOutput("legend_plot", height="260px", width="100%")
      )
    )
  })
  
  
  output$legend_plot <- renderPlot({
    req(lgnd())
    g <- lgnd()
    g$level <- factor(g$level, levels = as.character(g$level))
    ggplot(g, aes(x = level, y = 1, fill = level)) +
      geom_tile(color = "white") +
      scale_fill_manual(breaks = g$level, values = as.character(g$color)) +
      theme_void() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, color = "#3a4060", size = 9),
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
      ) +
      labs(x = NULL, y = NULL)
  }, bg = "transparent")
  
  vals_store <- reactiveValues(results = list())
  
  observeEvent(input$gen_gene, {
    
    genes_log <- input$selected_gene
    sig_label_log <- if (length(genes_log) > 1)
      paste0("SIG_", paste(genes_log, collapse="+"))
    else if (length(genes_log) == 1) genes_log[1]
    else return()
    
    v <- gnm()
    if(length(v[v %in% sig_label_log]) > 1){
      base_name <- paste0(sig_label_log, "_", length(v[v %in% sig_label_log]))
    } else {
      base_name <- sig_label_log
    }
    
    key <- isolate(base_name)
    
    a1 <- ifelse(odt() == "seu",
                 "Seurat spatial object (internally converted to NaVis)",
                 "NaVis or coercible object")
    a2 <- ifelse(ia() == 0,
                 "Object contained prior images already",
                 "Prior images calculated after object upload")
    a3 <- if(prms()[1]==0){
      paste0("The following parameters were used (sigma, iterations, threshold, sharpness): ",
             paste(c(0,"not used","not used","not used"), collapse = " "))
    }else{
      paste0("The following parameters were used (sigma, iterations, threshold, sharpness): ",
             paste(prms(), collapse = " "))             
    }
    a4 <- if(prms()[1]!=0){paste0("The following adaptive decay image was used: ",
                                  ifelse(is.null(primg()), "H&E Image", primg()))}else{
                                    paste0("The following adaptive decay image was used: ",
                                           "not used", primg())
                                  }
    a5 <- paste0("Original H&E image had scale: ",
                 orig()$microns_per_pixel, " microns/pixel")
    a6 <- paste0("Variable in the image is: ", sig_label_log)
    a7 <- mll()
    
    info <- list(
      source = a1,
      prior_images = a2,
      parameters = a3,
      adaptive_image = a4,
      scale = a5,
      variable = a6,
      misc = a7
    )
    
    vals_store$results[[key]] <- info
  })
  
  output$stored_info <- renderUI({
    req(input$selected_img)
    key <- input$selected_img
    info <- vals_store$results[[key]]
    if (is.null(info)) {
      return(tags$div(
        class = "alert alert-warning",
        "No saved values for this image yet."
      ))
    }else{
      z <- trl()
      info$threshold <- paste0("Image had threshold: : ",z[[key]])
    }
    
    
    tags$div(
      style = "display:flex; flex-wrap:wrap; gap:10px; justify-content:flex-start;",
      lapply(names(info), function(nm) {
        tags$div(
          class = "card shadow-sm border-0",
          style = "min-width:220px; max-width:350px; flex:1; 
                 background-color:#f9fafb; border-radius:12px;",
          tags$div(
            class = "card-body",
            tags$h6(class = "card-title text-primary",
                    style = "font-weight:600; margin-bottom:6px;",
                    nm),
            tags$p(
              class = "card-text text-muted",
              style = "white-space:pre-wrap; margin:0;",
              as.character(info[[nm]])
            )
          )
        )
      })
    )
  })
  
  observeEvent(input$allvals, {
    req(input$selected_img)
    key <- input$selected_img
    
    showModal(modalDialog(
      title = paste("Saved parameters for", key),
      uiOutput("stored_info"),
      easyClose = TRUE,
      size = "l"
    ))
  })
  
  
  output$dl_session_all <- downloadHandler(
    filename = function() {
      paste0("NaVis_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      files <- imgs()
      if (is.null(files) || length(files) == 0) {
        showNotification("Save all: no channels to export.",
                         type = "warning", duration = 4)
        writeLines("No channels were available at export time.",
                   file)
        return()
      }
      showModal(modalDialog(
        title = "Preparing session export",
        "Bundling channels and parameters into a ZIP...",
        footer = NULL, size = "s"
      ))
      on.exit(removeModal(), add = TRUE)
      
      tdir <- file.path(tempdir(), paste0("navis_session_", as.integer(Sys.time())))
      dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
      out_files <- character()
      
      
      .safe_name <- function(s) {
        s <- gsub("[^A-Za-z0-9._-]+", "_", s)
        s <- gsub("_+", "_", s)
        s <- gsub("^_+|_+$", "", s)
        if (!nzchar(s)) s <- "channel"
        s
      }
      
      manifest_rows <- list()
      vs <- isolate(vals_store$results)
      cs <- isolate(composite_sources())
      thr <- isolate(trl())
      
      for (nm in names(files)) {
        srcfile <- files[[nm]]
        if (is.null(srcfile) || !file.exists(srcfile)) next
        safe <- .safe_name(nm)
        dst <- file.path(tdir, paste0(safe, ".png"))
        ok <- tryCatch({ file.copy(srcfile, dst, overwrite = TRUE); TRUE },
                       error = function(e) FALSE)
        if (!ok) next
        out_files <- c(out_files, dst)
        
        
        srctype <- if (grepl("composite", nm, ignore.case = TRUE)) "composite"
        else if (grepl("overlay", nm, ignore.case = TRUE)) "h_and_e_overlay"
        else if (grepl("agreement", nm, ignore.case = TRUE)) "agreement_map"
        else if (grepl("interface", nm, ignore.case = TRUE)) "interface_map"
        else if (grepl("^META_SR_", nm)) "metadata_super_resolved"
        else if (grepl("^META_", nm)) "metadata_map"
        else if (grepl("^SAVED_", nm)) "saved_view"
        else if (nm %in% c("H&E Image","Entropy",
                           "Non-fibrillar interstitial tissue",
                           "Fibrillar interstitial tissue",
                           "Nuclei and cells")) "prior"
        else "channel"
        
        
        info <- vs[[nm]]
        info_str <- if (!is.null(info)) {
          paste(vapply(names(info), function(k)
            paste0(k, "=", as.character(info[[k]])),
            character(1)), collapse = "; ")
        } else ""
        
        
        sources <- if (!is.null(cs) && nm %in% names(cs)) {
          paste(cs[[nm]], collapse = " + ")
        } else ""
        
        
        thr_val <- if (!is.null(thr) && nm %in% names(thr)) thr[[nm]] else ""
        
        manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
          file = paste0(safe, ".png"),
          channel_name = nm,
          type = srctype,
          composite_sources = sources,
          intensity_threshold = thr_val,
          parameters = info_str,
          stringsAsFactors = FALSE
        )
      }
      
      if (length(manifest_rows) > 0) {
        manifest <- do.call(rbind, manifest_rows)
        write.csv(manifest, file.path(tdir, "manifest.csv"), row.names = FALSE)
        out_files <- c(out_files, file.path(tdir, "manifest.csv"))
      }
      
      
      readme_path <- file.path(tdir, "README.txt")
      writeLines(c(
        "NaVis session export",
        paste0("Exported at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("Number of channels: ", length(files)),
        "",
        "Contents:",
        "  *.png         - rendered channel images at display resolution",
        "  manifest.csv  - per-channel metadata, source type, parameters",
        "",
        "Channel types:",
        "  prior                     - H&E-derived tissue compartments (priors)",
        "  channel                   - single-gene or signature channels from Acquire",
        "  composite                 - merged channels from the Merge & Interface panel",
        "  agreement_map             - agreement maps from the Merge & Interface panel",
        "  interface_map             - juxtaposed-but-non-overlapping interface maps from the Merge & Interface panel",
        "  metadata_map              - metadata-column maps from the Metadata panel (raw/dilated)",
        "  metadata_super_resolved   - metadata-column maps from the Metadata panel rendered via AD-EBIDW",
        "  h_and_e_overlay           - committed channel + H&E blends from the Overlay panel",
        "  saved_view                - channels saved with their on-screen appearance via Save current view",
        "",
        "For reproducibility: the manifest 'parameters' column lists every",
        "parameter used at acquisition time (sigma, iterations, threshold,",
        "sharpness, decontamination alpha, dropout tau, etc.).",
        paste0("Source code and method details: https://github.com/Izzilab/NaVis")
      ), readme_path)
      out_files <- c(out_files, readme_path)
      
      
      old_wd <- getwd()
      setwd(tdir)
      on.exit(setwd(old_wd), add = TRUE)
      zip::zip(file, files = basename(out_files))
    },
    contentType = "application/zip"
  )
  
  
  output$dn1 <- downloadHandler(
    filename = function() {
      paste0(input$selected_img, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      
      req(input$selected_img)
      key <- input$selected_img
      
      showModal(modalDialog(
        "Preparing image and info for download and zipping...",
        footer = NULL
      ))
      
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      files <- imgs()
      fname <- files[[input$selected_img]]
      
      if (is.character(fname)) {
        fname <- EBImage::readImage(fname)
      }
      
      
      if (all(dim(fname)[1:2] < 800)) {
        d  <- dim(fname); ar <- d[1] / d[2]
        new_w <- if (ar >= 1) 2000 else round(2000 * ar)
        new_h <- if (ar >= 1) round(2000 / ar) else 2000
        fname <- EBImage::resize(fname, w = new_w, h = new_h)
      }
      
      fname <- apply_filters(fname, if (!is.null(input$intensity_threshold)) input$intensity_threshold else 0.30, crispr())
      
      
      png_file <- "image.png"
      EBImage::writeImage(fname, png_file)
      
      
      info <- vals_store$results[[key]]
      if (is.null(info)) {
        info <- list()
        info$source <- "No information stored for this image"
      }else{
        z <- trl()
        
        info$threshold <- paste0("Image had threshold: : ",z[[key]])
      }
      
      writeLines(
        c(
          paste("Dataset:", nmm()[[1]]),
          paste("Image:", input$selected_img),
          paste("Time:", Sys.time()),
          "Metadata (info list):",
          capture.output(str(info))
        ),
        "info.txt"
      )
      
      
      zip(file, c("image.png", "info.txt"))
      
      removeModal()
    }
  )
  
  output$dn2 <- downloadHandler(
    filename = function() {
      paste0("ALL_IMAGES_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      
      showModal(modalDialog(
        "Preparing all images and metadata for download (this might take some time)...",
        footer = NULL
      ))
      
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      files <- imgs()
      all_keys <- names(files)
      
      out_files <- c()   
      
      for (key in all_keys) {
        
        fname <- files[[key]]
        
        
        if (is.character(fname)) {
          fname <- EBImage::readImage(fname)
        }
        
        
        if (all(dim(fname)[1:2] < 800)) {
          d  <- dim(fname); ar <- d[1] / d[2]
          new_w <- if (ar >= 1) 2000 else round(2000 * ar)
          new_h <- if (ar >= 1) round(2000 / ar) else 2000
          fname <- EBImage::resize(fname, w = new_w, h = new_h)
        }
        
        fname <- apply_filters(fname, if (!is.null(input$intensity_threshold)) input$intensity_threshold else 0.30, crispr())
        
        
        png_file <- paste0(key, ".png")
        EBImage::writeImage(fname, png_file)
        
        
        info <- vals_store$results[[key]]
        
        if (is.null(info)) {
          info <- list(source = "No information stored for this image")
        } else {
          z <- trl()
          info$threshold <- paste0("Image had threshold: ", z[[key]])
        }
        
        info_file <- paste0(key, "_info.txt")
        writeLines(
          c(
            paste("Dataset:", nmm()[[1]]),
            paste("Image:", key),
            paste("Time:", Sys.time()),
            "Metadata (info list):",
            capture.output(str(info))
          ),
          info_file
        )
        
        
        out_files <- c(out_files, png_file, info_file)
      }
      
      
      zip(file, out_files)
      
      removeModal()
    }
  )
  
  output$dn3 <- downloadHandler(
    filename = function() {
      paste0("NaVis object_", nmm(), ".RDS")
    },
    content = function(file) {
      
      showModal(modalDialog(
        "Preparing the data object for download...",
        footer = NULL
      ))
      
      
      obj <- data_list()
      
      
      saveRDS(obj, file)
      
      removeModal()
    }
  )
  
  
  .crop_region <- function(img, cx, cy, r = 64) {
    if (!inherits(img, "Image")) img <- EBImage::Image(img)
    
    if (length(dim(img)) < 3 || dim(img)[3] < 3)
      img <- EBImage::rgbImage(img, img, img)
    
    iw <- dim(img)[1]   
    ih <- dim(img)[2]   
    x1 <- max(1, cx - r); x2 <- min(iw, cx + r)
    y1 <- max(1, cy - r); y2 <- min(ih, cy + r)
    img[x1:x2, y1:y2, ]
  }
  
  .nearest_spot <- function(data, img_x_norm, img_y_norm) {
    coords <- data$coords
    img_dims <- data$img_dims
    px <- img_x_norm * img_dims[1]
    py <- img_y_norm * img_dims[2]
    dx <- coords$x - px
    dy <- coords$y - py
    which.min(dx^2 + dy^2)
  }
  
  spot_data <- reactiveVal(NULL)
  
  .compute_spot_data <- function(x_n, y_n) {
    req(data_list(), imgs())
    lst   <- data_list()
    
    result <- tryCatch({
      
      he_path <- imgs()[["H&E Image"]]
      if (is.null(he_path)) return(NULL)
      he_img <- if (is.character(he_path)) EBImage::readImage(he_path) else he_path
      
      he_img <- EBImage::Image(he_img, colormode = "Color")
      if (dim(he_img)[3] < 3) he_img <- EBImage::rgbImage(he_img, he_img, he_img)
      if (all(dim(he_img)[1:2] < 800))
      {
        d <- dim(he_img); ar <- d[1] / d[2]
        new_w <- if (ar >= 1) 2000 else round(2000 * ar)
        new_h <- if (ar >= 1) round(2000 / ar) else 2000
        he_img <- EBImage::resize(he_img, w = new_w, h = new_h)
      }
      
      iw <- dim(he_img)[1]   
      ih <- dim(he_img)[2]   
      cx <- max(1L, min(iw, round(x_n * iw)))
      cy <- max(1L, min(ih, round(y_n * ih)))
      
      crop_r  <- max(32L, round(min(iw, ih) * 0.06))
      he_crop <- .crop_region(he_img, cx, cy, r = crop_r)
      
      r_vals <- as.numeric(imageData(he_crop)[,,1])
      g_vals <- as.numeric(imageData(he_crop)[,,2])
      b_vals <- as.numeric(imageData(he_crop)[,,3])
      hist_df <- data.frame(
        intensity = c(r_vals, g_vals, b_vals),
        channel   = rep(c("R","G","B"),
                        c(length(r_vals), length(g_vals), length(b_vals)))
      )
      
      get_patch <- function(nm, binarise = FALSE) {
        tryCatch({
          if (!nm %in% names(lst)) return(NULL)
          pm  <- lst[[nm]]
          mat <- if (inherits(pm, "Image")) imageData(pm) else pm
          if (length(dim(mat)) > 2) mat <- mat[,,1]
          mat <- as.matrix(mat)
          pr  <- nrow(mat); pc <- ncol(mat)
          mcx <- max(1L, min(pc, round(x_n * pc)))
          mcy <- max(1L, min(pr, round(y_n * pr)))
          r   <- max(3L, round(min(pr, pc) * 0.04))
          y1  <- max(1L, mcy - r); y2 <- min(pr, mcy + r)
          x1  <- max(1L, mcx - r); x2 <- min(pc, mcx + r)
          patch <- as.numeric(mat[y1:y2, x1:x2])
          
          if (binarise) patch <- as.numeric(patch > 0)
          patch
        }, error = function(e) NULL)
      }
      
      p_nuc  <- get_patch("nuclei_mask",         binarise = TRUE)    
      p_fib  <- get_patch("fibrous_tissue_mask", binarise = FALSE)   
      p_soft <- get_patch("soft_tissue_mask",    binarise = FALSE)   
      p_ent  <- get_patch("entropy_mask",        binarise = FALSE)
      
      m_nuc  <- if (!is.null(p_nuc))  mean(p_nuc,  na.rm = TRUE) else 0
      m_fib  <- if (!is.null(p_fib))  mean(p_fib,  na.rm = TRUE) else 0
      m_soft <- if (!is.null(p_soft)) mean(p_soft, na.rm = TRUE) else 0
      if (!is.finite(m_nuc))  m_nuc  <- 0
      if (!is.finite(m_fib))  m_fib  <- 0
      if (!is.finite(m_soft)) m_soft <- 0
      
      total <- m_nuc + m_fib + m_soft
      if (total > 1e-6) {
        frac_nuc  <- m_nuc  / total
        frac_fib  <- m_fib  / total
        frac_soft <- m_soft / total
        frac_bg   <- 0
      } else {
        
        frac_nuc <- 0; frac_fib <- 0; frac_soft <- 0; frac_bg <- 1
      }
      
      score_ent <- if (!is.null(p_ent)) mean(p_ent, na.rm=TRUE) else NA_real_
      
      prior_df <- data.frame(
        label = c("Nuclei / Cells", "Fibrillar ECM",
                  "Soft Tissue",    "Background"),
        value = c(frac_nuc, frac_fib, frac_soft, frac_bg),
        color = c("#4488ff", "#ff8844", "#44bb88", "#cccccc"),
        stringsAsFactors = FALSE
      )
      
      attr(prior_df, "entropy") <- score_ent
      
      prior_names  <- c("nuclei_mask","fibrous_tissue_mask","soft_tissue_mask","entropy_mask")
      prior_labels <- c("Nuclei / Cells","Fibrillar ECM","Soft Tissue","Entropy")
      prior_colors <- c("#4488ff","#ff8844","#44bb88","#aa44ff")
      
      sel       <- input$selected_img
      gene_val  <- NA_real_
      gene_label <- if (!is.null(sel)) sel else "—"
      if (!is.null(sel) && sel %in% names(imgs())) {
        gf <- tryCatch({
          p <- imgs()[[sel]]
          img <- if (is.character(p)) EBImage::readImage(p) else p
          if (all(dim(img)[1:2] < 800)) {
            d <- dim(img); ar <- d[1] / d[2]
            new_w <- if (ar >= 1) 2000 else round(2000 * ar)
            new_h <- if (ar >= 1) round(2000 / ar) else 2000
            img <- EBImage::resize(img, w = new_w, h = new_h)
          }
          img
        }, error = function(e) NULL)
        if (!is.null(gf)) {
          gh  <- dim(gf)[1]; gw <- dim(gf)[2]
          gcx <- max(1L, min(gw, round(x_n * gw)))
          gcy <- max(1L, min(gh, round(y_n * gh)))
          gene_val <- tryCatch({
            mat <- imageData(gf)
            if (length(dim(mat)) == 2) as.numeric(mat[gcy, gcx])
            else mean(as.numeric(mat[gcy, gcx, ]), na.rm = TRUE)
          }, error = function(e) NA_real_)
        }
      }
      
      spot_idx <- tryCatch(.nearest_spot(lst, x_n, y_n), error = function(e) NA_integer_)
      spot_name <- NA_character_
      spot_genes_df <- NULL   
      
      if (!is.na(spot_idx)) {
        coords <- lst$coords
        if ("barcode" %in% colnames(coords))
          spot_name <- coords$barcode[spot_idx]
        
        known_genes <- rownames(lst$expression)
        is_composite <- grepl("^(composite|agreement|HEoverlay|interface)_", sel, ignore.case=TRUE)
        is_signature <- grepl("^SIG_", sel)
        
        all_img_genes <- if (is_composite) {
          comp_srcs <- composite_sources()[[sel]]
          if (!is.null(comp_srcs) && length(comp_srcs) >= 1) {
            unique(unlist(lapply(comp_srcs, function(s) .resolve_one(s, known_genes))))
          } else {
            .resolve_one(sel, known_genes)
          }
        } else if (is_signature) {
          .resolve_one(sel, known_genes)
        } else {
          if (sel %in% known_genes) {
            sel
          } else if (length(input$selected_gene) > 0 && nchar(input$selected_gene[1]) > 0) {
            input$selected_gene[1]
          } else {
            character(0)
          }
        }
        
        if (length(all_img_genes) > 0) {
          
          chan_path <- imgs()[[sel]]
          chan_val  <- NA_real_
          if (!is.null(chan_path) && file.exists(chan_path)) {
            chan_img <- tryCatch(EBImage::readImage(chan_path),
                                 error = function(e) NULL)
            if (!is.null(chan_img)) {
              d_chan <- dim(chan_img)
              
              px_chan <- max(1L, min(d_chan[1], round(x_n * d_chan[1])))
              py_chan <- max(1L, min(d_chan[2], round(y_n * d_chan[2])))
              arr <- imageData(chan_img)
              if (length(dim(arr)) == 3 && dim(arr)[3] >= 3) {
                chan_val <- max(arr[px_chan, py_chan, 1],
                                arr[px_chan, py_chan, 2],
                                arr[px_chan, py_chan, 3])
              } else if (length(dim(arr)) == 2) {
                chan_val <- arr[px_chan, py_chan]
              }
            }
          }
          
          spot_genes_df <- data.frame(
            gene = all_img_genes,
            expr = rep(chan_val, length(all_img_genes)),
            cpm  = rep(chan_val, length(all_img_genes)),
            stringsAsFactors = FALSE
          )
        }
      }
      
      spot_expr <- if (!is.null(spot_genes_df)) spot_genes_df$expr[1] else NA_real_
      spot_gene <- if (!is.null(spot_genes_df)) spot_genes_df$gene[1] else NA_character_
      
      corr_df <- NULL
      if (!is.null(sel) && sel %in% names(imgs())) {
        gf2 <- tryCatch({
          p <- imgs()[[sel]]
          if (is.character(p)) EBImage::readImage(p) else p
        }, error = function(e) NULL)
        if (!is.null(gf2)) {
          mat2 <- imageData(gf2)
          gvec <- if (length(dim(mat2)) == 2) as.numeric(mat2)
          else as.numeric(mat2[,,1] + mat2[,,2] + mat2[,,3]) / 3
          
          corrs <- sapply(prior_names, function(nm) {
            tryCatch({
              if (!nm %in% names(lst)) return(NA_real_)
              pm  <- lst[[nm]]
              mat <- if (inherits(pm, "Image")) imageData(pm) else pm
              if (length(dim(mat)) > 2) mat <- mat[,,1]
              pvec <- as.numeric(mat)
              
              if (length(pvec) != length(gvec)) {
                pimg <- EBImage::Image(matrix(pvec, nrow = nrow(mat), ncol = ncol(mat)))
                pimg <- EBImage::resize(pimg, w = dim(gf2)[2], h = dim(gf2)[1])
                pvec <- as.numeric(imageData(pimg))
              }
              suppressWarnings(cor(gvec, pvec, use = "complete.obs"))
            }, error = function(e) NA_real_)
          })
          corr_df <- data.frame(
            label = prior_labels, corr = as.numeric(corrs),
            color = prior_colors, stringsAsFactors = FALSE
          )
        }
      }
      
      list(he_crop=he_crop, hist_df=hist_df, prior_df=prior_df,
           gene_val=gene_val, gene_label=gene_label,
           spot_expr=spot_expr, spot_gene=spot_gene,
           spot_genes_df=spot_genes_df, spot_name=spot_name,
           corr_df=corr_df, x_norm=x_n, y_norm=y_n,
           microns_per_pixel=orig()$microns_per_pixel)
      
    }, error = function(e) {
      warning(sprintf("compute_spot_data failed at (%.3f, %.3f): %s",
                      x_n, y_n, conditionMessage(e)),
              call. = FALSE)
      showNotification(
        paste0("Could not analyse this region: ", conditionMessage(e)),
        type = "error", duration = 6)
      NULL
    })
    
    result
  }
  
  compare_state <- reactiveVal(list(stage = "empty",
                                    first  = NULL,   
                                    second = NULL,   
                                    channel_at_store = NULL))
  
  .reset_compare <- function() {
    compare_state(list(stage = "empty", first = NULL, second = NULL,
                       channel_at_store = NULL))
    session$sendCustomMessage("clearCompareMarkers", list())
  }
  
  observeEvent(input$selected_img, {
    cs <- compare_state()
    if (cs$stage != "empty") {
      .reset_compare()
    }
  }, ignoreInit = TRUE)
  
  roi_state <- reactiveVal(list(active = FALSE, data = NULL,
                                corners = NULL))
  
  .reset_roi <- function() {
    roi_state(list(active = FALSE, data = NULL, corners = NULL))
    session$sendCustomMessage("clearROIOverlay", list())
    session$sendCustomMessage("toggleROIDrawer", list(roi = FALSE))
  }
  
  observeEvent(input$roi_clear, { .reset_roi() })
  
  .compute_roi_data <- function(corners) {
    req(data_list(), imgs())
    lst <- data_list()
    
    iw <- lst$img_dims[1]
    ih <- lst$img_dims[2]
    px1 <- min(corners$x1, corners$x2) * iw
    px2 <- max(corners$x1, corners$x2) * iw
    py1 <- min(corners$y1, corners$y2) * ih
    py2 <- max(corners$y1, corners$y2) * ih
    
    inside_mask <- lst$coords$x >= px1 & lst$coords$x <= px2 &
      lst$coords$y >= py1 & lst$coords$y <= py2
    n_in <- sum(inside_mask)
    if (n_in < 5L) {
      showNotification(sprintf("ROI contains only %d spots — needs at least 5.",
                               n_in),
                       type = "warning", duration = 4)
      return(NULL)
    }
    
    spot_idx_in  <- which(inside_mask)
    spot_idx_out <- which(!inside_mask)
    
    he_path <- imgs()[["H&E Image"]]
    he_img  <- if (is.character(he_path)) EBImage::readImage(he_path) else he_path
    he_img  <- EBImage::Image(he_img, colormode = "Color")
    if (length(dim(he_img)) < 3 || dim(he_img)[3] < 3)
      he_img <- EBImage::rgbImage(he_img, he_img, he_img)
    if (all(dim(he_img)[1:2] < 800)) {
      d <- dim(he_img); ar <- d[1] / d[2]
      new_w <- if (ar >= 1) 2000 else round(2000 * ar)
      new_h <- if (ar >= 1) round(2000 / ar) else 2000
      he_img <- EBImage::resize(he_img, w = new_w, h = new_h)
    }
    crop_iw <- dim(he_img)[1]; crop_ih <- dim(he_img)[2]
    cx1 <- max(1L, round(min(corners$x1, corners$x2) * crop_iw))
    cx2 <- min(crop_iw, round(max(corners$x1, corners$x2) * crop_iw))
    cy1 <- max(1L, round(min(corners$y1, corners$y2) * crop_ih))
    cy2 <- min(crop_ih, round(max(corners$y1, corners$y2) * crop_ih))
    he_crop <- he_img[cx1:cx2, cy1:cy2, ]
    
    .roi_density <- function(mat, binarise = FALSE) {
      a <- if (inherits(mat, "Image")) imageData(mat) else mat
      if (length(dim(a)) > 2) a <- a[,,1]
      a <- as.matrix(a)
      mw <- dim(a)[1]; mh <- dim(a)[2]
      bx1 <- max(1L, round(min(corners$x1, corners$x2) * mw))
      bx2 <- min(mw, round(max(corners$x1, corners$x2) * mw))
      by1 <- max(1L, round(min(corners$y1, corners$y2) * mh))
      by2 <- min(mh, round(max(corners$y1, corners$y2) * mh))
      patch <- a[bx1:bx2, by1:by2]
      if (binarise) patch <- as.numeric(patch > 0)
      mean(patch, na.rm = TRUE)
    }
    rho_nuc <- tryCatch(.roi_density(lst$nuclei_mask, binarise = TRUE),  error=function(e) NA_real_)
    rho_fib <- tryCatch(.roi_density(lst$fibrous_tissue_mask, binarise = FALSE), error=function(e) NA_real_)
    rho_sof <- tryCatch(.roi_density(lst$soft_tissue_mask,    binarise = FALSE), error=function(e) NA_real_)
    rhos <- c(rho_nuc, rho_fib, rho_sof)
    rhos[!is.finite(rhos)] <- 0
    if (sum(rhos) > 0) rhos <- rhos / sum(rhos)
    prior_df <- data.frame(
      label = c("Nuclei and cells", "Fibrillar interstitial tissue",
                "Non-fibrillar interstitial tissue"),
      value = rhos,
      color = c("#4488ff", "#ff8844", "#44bb88"),
      stringsAsFactors = FALSE
    )
    
    expr <- lst$expression
    
    n_expr <- ncol(expr)
    if ("barcode" %in% colnames(lst$coords) &&
        !is.null(colnames(expr)) &&
        any(lst$coords$barcode %in% colnames(expr))) {
      bc_in  <- lst$coords$barcode[spot_idx_in]
      bc_out <- lst$coords$barcode[spot_idx_out]
      col_idx_in  <- which(colnames(expr) %in% bc_in)
      col_idx_out <- which(colnames(expr) %in% bc_out)
    } else {
      
      col_idx_in  <- spot_idx_in[spot_idx_in  <= n_expr]
      col_idx_out <- spot_idx_out[spot_idx_out <= n_expr]
      if (length(col_idx_in) != length(spot_idx_in) ||
          length(col_idx_out) != length(spot_idx_out)) {
        
      }
    }
    if (length(col_idx_in) < 5L || length(col_idx_out) < 1L) {
      showNotification(
        sprintf("ROI: only %d in-tissue spots could be matched to expression data.",
                length(col_idx_in)),
        type = "warning", duration = 4)
      gene_df <- data.frame(gene = character(), ratio = numeric(),
                            mean_in = numeric(), mean_out = numeric(),
                            stringsAsFactors = FALSE)
    } else {
      mean_in  <- Matrix::rowMeans(expr[, col_idx_in,  drop = FALSE])
      mean_out <- Matrix::rowMeans(expr[, col_idx_out, drop = FALSE])
      eps <- 1e-3
      ratio <- (mean_in + eps) / (mean_out + eps)
      keep <- mean_in >= 0.5 & is.finite(ratio)
      if (sum(keep) < 1) {
        gene_df <- data.frame(gene = character(), ratio = numeric(),
                              mean_in = numeric(), mean_out = numeric(),
                              stringsAsFactors = FALSE)
      } else {
        ord <- order(ratio[keep], decreasing = TRUE)[1:min(20, sum(keep))]
        gene_names <- rownames(expr)[keep][ord]
        gene_df <- data.frame(
          gene     = gene_names,
          ratio    = ratio[keep][ord],
          mean_in  = mean_in[keep][ord],
          mean_out = mean_out[keep][ord],
          stringsAsFactors = FALSE
        )
      }
    }
    
    list(
      n_spots   = n_in,
      corners   = corners,
      he_crop   = he_crop,
      prior_df  = prior_df,
      gene_df   = gene_df,
      area_um2  = {
        
        mpp <- if (!is.null(orig()$microns_per_pixel)) orig()$microns_per_pixel else NA_real_
        if (is.finite(mpp))
          (px2 - px1) * (py2 - py1) * mpp^2
        else NA_real_
      }
    )
  }
  
  observeEvent(input$stage_roi, {
    corners <- input$stage_roi
    req(corners)
    withProgress(message = "Computing ROI readouts…", value = 0.5, {
      d <- .compute_roi_data(corners)
    })
    if (is.null(d)) return()
    roi_state(list(active = TRUE, data = d, corners = corners))
    session$sendCustomMessage("toggleROIDrawer", list(roi = TRUE))
  })
  
  observeEvent(input$stage_click, {
    click <- input$stage_click
    
    cs <- compare_state()
    if (cs$stage == "compared") .reset_compare()
    
    result <- .compute_spot_data(click$x, click$y)
    if (!is.null(result)) {
      spot_data(result)
      session$sendCustomMessage("showAnalysisResults", list())
    }
  })
  
  observeEvent(input$stage_right_click, {
    click <- input$stage_right_click
    if (is.null(click) || is.null(click$x) || is.null(click$y)) return()
    
    sd <- .compute_spot_data(click$x, click$y)
    if (is.null(sd)) {
      showNotification("Could not analyse this region — try clicking on tissue.",
                       type = "warning", duration = 4)
      return()
    }
    
    cs <- compare_state()
    
    if (cs$stage == "empty") {
      
      compare_state(list(stage = "first_stored",
                         first  = sd,
                         second = NULL,
                         channel_at_store = isolate(input$selected_img)))
      session$sendCustomMessage("showCompareMarker",
                                list(n = 1, x_norm = click$x, y_norm = click$y))
      showNotification("Region A stored — right-click another region to compare.",
                       type = "message", duration = 4)
    } else if (cs$stage == "first_stored") {
      
      compare_state(list(stage = "compared",
                         first  = cs$first,
                         second = sd,
                         channel_at_store = cs$channel_at_store))
      session$sendCustomMessage("showCompareMarker",
                                list(n = 2, x_norm = click$x, y_norm = click$y))
      session$sendCustomMessage("showAnalysisResults", list())
    } else if (cs$stage == "compared") {
      
      .reset_compare()
      compare_state(list(stage = "first_stored",
                         first  = sd,
                         second = NULL,
                         channel_at_store = isolate(input$selected_img)))
      session$sendCustomMessage("showCompareMarker",
                                list(n = 1, x_norm = click$x, y_norm = click$y))
      showNotification("New comparison started — right-click another region.",
                       type = "message", duration = 3)
    }
  })
  
  observeEvent(input$compare_clear, {
    .reset_compare()
  })
  
  output$spot_coords_ui <- renderUI({
    sd <- spot_data(); req(sd)
    mpp <- sd$microns_per_pixel
    if (is.null(mpp) || is.na(mpp)) mpp <- 1
    div(class = "spot-card",
        div(class = "spot-card-title", "Stage position (normalised)"),
        div(class = "spot-card-value",
            style = "font-size:13px;",
            sprintf("x = %.3f,  y = %.3f", sd$x_norm, sd$y_norm)),
        div(class = "spot-card-sub",
            sprintf("Scale: %.3f µm/px", mpp))
    )
  })
  
  output$spot_crop <- renderImage({
    sd <- spot_data(); req(sd)
    tmp <- tempfile(fileext = ".png")
    EBImage::writeImage(sd$he_crop, tmp)
    list(src = tmp, contentType = "image/png",
         width = "100%", style = "border-radius:4px; border:1px solid #c8ccd8;")
  }, deleteFile = TRUE)
  
  .plot_spot_histogram <- function(df, base_size = 9) {
    ch_cols <- c("R" = "#cc3344", "G" = "#33aa55", "B" = "#3366cc")
    ggplot(df, aes(x = intensity, fill = channel, color = channel)) +
      geom_density(alpha = 0.45, linewidth = 0.5) +
      scale_fill_manual(values = ch_cols) +
      scale_color_manual(values = ch_cols) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      labs(x = "Pixel intensity", y = NULL) +
      theme_void(base_size = base_size) +
      theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position  = "top",
        legend.title     = element_blank(),
        legend.text      = element_text(color = "#4a5070", size = base_size),
        axis.text.x      = element_text(color = "#7880a0", size = base_size - 1),
        axis.line.x      = element_line(color = "#d0d4e0")
      )
  }
  
  .plot_spot_correlations <- function(df, base_size = 9) {
    if (is.null(df) || all(is.na(df$corr))) {
      return(ggplot() + theme_void() +
               labs(title = "No data") +
               theme(plot.background = element_rect(fill = "transparent", color = NA),
                     title = element_text(color = "#7880a0", size = base_size)))
    }
    df$corr[is.na(df$corr)] <- 0
    df$label <- factor(df$label, levels = rev(df$label))
    ggplot(df, aes(x = label, y = corr, fill = color)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0, color = "#d0d4e0", linewidth = 0.5) +
      scale_fill_identity() +
      scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
      coord_flip() +
      labs(x = NULL, y = "Pearson r", title = "Gene channel vs tissue priors") +
      theme_minimal(base_size = base_size) +
      theme(
        plot.background    = element_rect(fill = "transparent", color = NA),
        panel.background   = element_rect(fill = "transparent", color = NA),
        panel.grid.major.x = element_line(color = "#e0e4ef", linewidth = 0.4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text          = element_text(color = "#4a5070", size = base_size),
        axis.title         = element_text(color = "#7880a0", size = base_size),
        plot.title         = element_text(color = "#7880a0", size = base_size, face = "plain")
      )
  }
  
  output$spot_histogram <- renderPlot({
    sd <- spot_data(); req(sd)
    .plot_spot_histogram(sd$hist_df, base_size = 9)
  }, bg = "transparent")
  
  output$spot_expression_ui <- renderUI({
    sd <- spot_data(); req(sd)
    cards <- list(
      div(class = "spot-card",
          div(class = "spot-card-title", "Channel pixel value"),
          div(class = "spot-card-value",
              if (is.na(sd$gene_val)) "—" else sprintf("%.3f", sd$gene_val)),
          div(class = "spot-card-sub", sd$gene_label)
      )
    )
    
    df <- sd$spot_genes_df
    if (!is.null(df) && nrow(df) > 0 && !is.na(sd$spot_name)) {
      gene_cards <- lapply(seq_len(nrow(df)), function(i) {
        raw <- df$expr[i]
        cpm <- if ("cpm" %in% names(df)) df$cpm[i] else NA_real_
        div(class = "spot-card",
            div(class = "spot-card-title",
                paste0("Nearest spot — ", df$gene[i])),
            div(class = "spot-card-value",
                if (is.na(cpm)) {
                  if (is.na(raw)) "—" else sprintf("%.1f counts", raw)
                } else {
                  sprintf("%.1f CPM", cpm)
                }),
            div(class = "spot-card-sub",
                if (is.na(raw)) ""
                else sprintf("raw: %.0f  ·  %s",
                             raw, substr(sd$spot_name, 1, 14)))
        )
      })
      cards <- c(cards, gene_cards)
    }
    do.call(div, cards)
  })
  
  output$spot_priors_ui <- renderUI({
    sd <- spot_data(); req(sd)
    df  <- sd$prior_df
    ent <- attr(sd$prior_df, "entropy")
    
    tagList(
      
      tagList(lapply(seq_len(nrow(df)), function(i) {
        val <- if (is.na(df$value[i])) 0 else df$value[i]
        pct <- round(val * 100)
        div(class = "prior-row",
            div(class = "prior-label", df$label[i]),
            div(class = "prior-bar-bg",
                div(class = "prior-bar-fill",
                    style = sprintf("width:%d%%; background:%s;", pct, df$color[i]))
            ),
            div(class = "prior-value", paste0(pct, "%"))
        )
      })),
      
      if (!is.null(ent) && !is.na(ent)) {
        div(style = "margin-top:8px; padding-top:6px; border-top:1px solid #e8eaf0;",
            div(class = "prior-row",
                div(class = "prior-label",
                    style = "color:#aa44ff; font-weight:600;",
                    "Texture (entropy)"),
                div(class = "prior-bar-bg",
                    div(class = "prior-bar-fill",
                        style = sprintf("width:%d%%; background:#aa44ff;",
                                        round(ent * 100)))
                ),
                div(class = "prior-value",
                    style = "color:#aa44ff;",
                    sprintf("%.2f", ent))
            )
        )
      }
    )
  })
  
  output$spot_correlations <- renderPlot({
    sd <- spot_data(); req(sd)
    .plot_spot_correlations(sd$corr_df, base_size = 9)
  }, bg = "transparent")
  
  output$compare_drawer_ui <- renderUI({
    cs <- compare_state()
    req(cs$stage == "compared", cs$first, cs$second)
    sdA <- cs$first; sdB <- cs$second
    
    fnum <- function(x, fmt = "%.3f")
      if (is.null(x) || !is.finite(x)) "—" else sprintf(fmt, x)
    
    fdelta <- function(a, b, suffix = "") {
      if (is.null(a) || is.null(b) || !is.finite(a) || !is.finite(b))
        return(tags$span(style = "color:#7880a0;", "—"))
      eps <- 1e-6
      if (a + eps <= 0 && b + eps <= 0)
        return(tags$span(style = "color:#7880a0;", "—"))
      lfc <- log2((b + eps) / (a + eps))
      if (abs(lfc) <= 0.1)
        return(tags$span(style = "color:#7880a0; font-weight:600;", "—"))
      col <- if (lfc > 0) "#1e7a45" else "#c85050"
      sym <- if (lfc > 0) "▲" else "▼"
      tags$span(style = sprintf("color:%s; font-weight:600;", col),
                sprintf("%s %+.2f", sym, lfc))
    }
    
    crop_section <- function(sd, label, color) {
      tmp <- tempfile(fileext = ".png")
      EBImage::writeImage(sd$he_crop, tmp)
      img_data <- base64enc::dataURI(file = tmp, mime = "image/png")
      div(style = "text-align:center;",
          div(style = sprintf("font-size:11px; font-weight:700; color:%s;
                             margin-bottom:4px;", color),
              label),
          tags$img(src = img_data,
                   style = "max-width:100%; height:90px; border-radius:4px;
                          border:2px solid #c8ccd8; object-fit:cover;")
      )
    }
    
    pA <- if (!is.null(sdA$prior_df)) sdA$prior_df else
      data.frame(label = character(), value = numeric())
    pB <- if (!is.null(sdB$prior_df)) sdB$prior_df else
      data.frame(label = character(), value = numeric())
    common_priors <- union(pA$label, pB$label)
    tissue_rows <- lapply(common_priors, function(lbl) {
      
      vA_frac <- if (lbl %in% pA$label) pA$value[pA$label == lbl][1] else NA_real_
      vB_frac <- if (lbl %in% pB$label) pB$value[pB$label == lbl][1] else NA_real_
      vA <- if (is.finite(vA_frac)) 100 * vA_frac else NA_real_
      vB <- if (is.finite(vB_frac)) 100 * vB_frac else NA_real_
      div(style = "display:grid; grid-template-columns: 1fr auto 1fr auto;
                   gap:8px; padding:3px 6px; font-size:11px;
                   border-bottom:1px solid #f0f2f6; align-items:center;",
          div(style = "color:#5a6080;",
              if (is.finite(vA)) sprintf("%.1f%%", vA) else "—"),
          div(style = "color:#3a4060; font-weight:500;", lbl),
          div(style = "color:#5a6080; text-align:right;",
              if (is.finite(vB)) sprintf("%.1f%%", vB) else "—"),
          div(style = "text-align:right;", fdelta(vA, vB))
      )
    })
    
    gA <- sdA$spot_genes_df; gB <- sdB$spot_genes_df
    gene_rows <- if (!is.null(gA) && nrow(gA) && !is.null(gB) && nrow(gB)) {
      genes <- union(gA$gene, gB$gene)
      lapply(genes, function(g) {
        vA <- if (g %in% gA$gene) gA$cpm[gA$gene == g][1] else NA_real_
        vB <- if (g %in% gB$gene) gB$cpm[gB$gene == g][1] else NA_real_
        div(style = "display:grid; grid-template-columns: 1fr auto 1fr auto;
                     gap:8px; padding:3px 6px; font-size:11px;
                     border-bottom:1px solid #f0f2f6; align-items:center;
                     font-family:'Consolas',monospace;",
            div(style = "color:#5a6080;", fnum(vA, "%.3f")),
            div(style = "color:#3a4060; font-weight:500;
                       font-family:-apple-system,system-ui,sans-serif;", g),
            div(style = "color:#5a6080; text-align:right;", fnum(vB, "%.3f")),
            div(style = "text-align:right;", fdelta(vA, vB))
        )
      })
    } else NULL
    
    tagList(
      
      div(style = "padding:10px 12px; background:#f5f7fb;
                   border-radius:6px; margin-bottom:12px;
                   display:flex; align-items:center; justify-content:space-between;",
          div(
            div(style = "font-size:12px; font-weight:700; color:#2a3050;",
                "🔬 Comparing two regions"),
            div(style = "font-size:11px; color:#7880a0; margin-top:2px;",
                sprintf("Channel: %s", cs$channel_at_store %||% "—"))
          ),
          actionButton("compare_clear", "✕ Clear",
                       class = "btn-fp-icon",
                       style = "padding:4px 10px; font-size:11px;")
      ),
      
      div(style = "display:grid; grid-template-columns:1fr 1fr; gap:10px;
                   margin-bottom:14px;",
          crop_section(sdA, "Region A", "#2c5be0"),
          crop_section(sdB, "Region B", "#ff8844")
      ),
      
      div(style = "display:grid; grid-template-columns: 1fr 1fr; gap:10px;
                   margin-bottom:12px; font-size:10px; color:#7880a0;
                   text-align:center;",
          div(sprintf("(%.3f, %.3f)", sdA$x_norm, sdA$y_norm)),
          div(sprintf("(%.3f, %.3f)", sdB$x_norm, sdB$y_norm))
      ),
      
      div(class = "fp-section", "🧱 Tissue Context"),
      div(style = "font-size:10px; color:#7880a0; padding:4px 6px;
                   display:grid; grid-template-columns:1fr auto 1fr auto; gap:8px;",
          div(style = "color:#2c5be0;", "A"),
          div(""),
          div(style = "color:#ff8844; text-align:right;", "B"),
          div(style = "text-align:right;", "log₂(B/A)")
      ),
      tagList(tissue_rows),
      
      if (!is.null(gene_rows) && length(gene_rows)) {
        tagList(
          div(class = "fp-section", style = "margin-top:12px;",
              "🧬 Expression at Spot"),
          div(style = "font-size:10px; color:#7880a0; padding:4px 6px;
                       display:grid; grid-template-columns:1fr auto 1fr auto; gap:8px;",
              div(style = "color:#2c5be0;", "A"),
              div(""),
              div(style = "color:#ff8844; text-align:right;", "B"),
              div(style = "text-align:right;", "log₂(B/A)")
          ),
          tagList(gene_rows)
        )
      },
      
      div(style = "margin-top:14px; text-align:center;",
          downloadButton("dl_compare_all", "💾 Save comparison as ZIP",
                         class = "btn-fp", style = "width:auto;")
      )
    )
  })
  
  output$dl_compare_all <- downloadHandler(
    filename = function() {
      sprintf("navis_compare_%s.zip",
              .safe_name(isolate(compare_state()$channel_at_store)))
    },
    content = function(zipfile) {
      cs <- compare_state(); req(cs$stage == "compared")
      sdA <- cs$first; sdB <- cs$second
      tdir <- tempfile(); dir.create(tdir)
      on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
      
      EBImage::writeImage(sdA$he_crop, file.path(tdir, "regionA_crop.png"))
      EBImage::writeImage(sdB$he_crop, file.path(tdir, "regionB_crop.png"))
      write.csv(sdA$prior_df %||% data.frame(),
                file.path(tdir, "regionA_tissue_context.csv"), row.names = FALSE)
      write.csv(sdB$prior_df %||% data.frame(),
                file.path(tdir, "regionB_tissue_context.csv"), row.names = FALSE)
      
      gA <- sdA$spot_genes_df; gB <- sdB$spot_genes_df
      if (!is.null(gA) && !is.null(gB) && nrow(gA) && nrow(gB)) {
        genes <- union(gA$gene, gB$gene)
        out <- data.frame(
          gene  = genes,
          regA_cpm  = vapply(genes, function(g) {
            if (g %in% gA$gene) gA$cpm[gA$gene == g][1] else NA_real_
          }, numeric(1)),
          regB_cpm  = vapply(genes, function(g) {
            if (g %in% gB$gene) gB$cpm[gB$gene == g][1] else NA_real_
          }, numeric(1)),
          regA_raw  = vapply(genes, function(g) {
            if (g %in% gA$gene) gA$expr[gA$gene == g][1] else NA_real_
          }, numeric(1)),
          regB_raw  = vapply(genes, function(g) {
            if (g %in% gB$gene) gB$expr[gB$gene == g][1] else NA_real_
          }, numeric(1)),
          stringsAsFactors = FALSE)
        eps <- 1e-6
        out$log2_FC_B_vs_A <- log2((out$regB_cpm + eps) / (out$regA_cpm + eps))
        write.csv(out, file.path(tdir, "expression_comparison.csv"), row.names = FALSE)
      }
      
      writeLines(c(
        "NaVis two-region comparison bundle",
        sprintf("  Channel:    %s", cs$channel_at_store %||% "—"),
        sprintf("  Region A:   (%.3f, %.3f)", sdA$x_norm, sdA$y_norm),
        sprintf("  Region B:   (%.3f, %.3f)", sdB$x_norm, sdB$y_norm),
        "",
        "Files:",
        "  regionA_crop.png            H&E crop centred on region A",
        "  regionB_crop.png            H&E crop centred on region B",
        "  regionA_tissue_context.csv  per-compartment fractions at A",
        "  regionB_tissue_context.csv  per-compartment fractions at B",
        "  expression_comparison.csv   gene values at A, at B, and log₂(B/A) fold-change"
      ), file.path(tdir, "summary.txt"))
      
      zip::zipr(zipfile, files = list.files(tdir, full.names = TRUE))
    }
  )
  
  observeEvent(compare_state(), {
    cs <- compare_state()
    if (cs$stage == "compared") {
      session$sendCustomMessage("toggleCompareDrawer",
                                list(compare = TRUE))
    } else {
      session$sendCustomMessage("toggleCompareDrawer",
                                list(compare = FALSE))
    }
  })
  
  output$roi_drawer_ui <- renderUI({
    rs <- roi_state()
    req(rs$active, !is.null(rs$data))
    d <- rs$data
    
    fnum <- function(x, fmt = "%.2f")
      if (is.null(x) || !is.finite(x)) "—" else sprintf(fmt, x)
    
    tmp <- tempfile(fileext = ".png")
    he_for_save <- d$he_crop
    if (max(dim(he_for_save)[1:2]) > 600) {
      sc <- 600 / max(dim(he_for_save)[1:2])
      he_for_save <- EBImage::resize(he_for_save,
                                     w = round(dim(he_for_save)[1] * sc),
                                     h = round(dim(he_for_save)[2] * sc))
    }
    EBImage::writeImage(he_for_save, tmp)
    img_data <- base64enc::dataURI(file = tmp, mime = "image/png")
    
    pp <- d$prior_df
    tissue_rows <- lapply(seq_len(nrow(pp)), function(i) {
      pct <- 100 * pp$value[i]
      div(style = "display:grid; grid-template-columns:1fr auto; gap:8px;
                   padding:3px 6px; font-size:11px;
                   border-bottom:1px solid #f0f2f6; align-items:center;",
          div(style = sprintf("color:%s; font-weight:500;", pp$color[i]),
              pp$label[i]),
          div(style = "color:#3a4060; font-weight:600;",
              if (is.finite(pct)) sprintf("%.1f%%", pct) else "—")
      )
    })
    
    gd <- d$gene_df
    gene_grid <- if (nrow(gd) > 0) {
      cells <- list()
      
      cells[[length(cells)+1]] <- div(
        style = "font-size:9px; color:#7880a0; padding:6px 8px 4px;
                 text-transform:uppercase; letter-spacing:0.5px;
                 border-bottom:1px solid #d0d4e0;",
        "gene")
      cells[[length(cells)+1]] <- div(
        style = "font-size:9px; color:#7880a0; padding:6px 8px 4px;
                 text-transform:uppercase; letter-spacing:0.5px;
                 border-bottom:1px solid #d0d4e0; text-align:right;",
        "ratio")
      cells[[length(cells)+1]] <- div(
        style = "font-size:9px; color:#7880a0; padding:6px 8px 4px;
                 text-transform:uppercase; letter-spacing:0.5px;
                 border-bottom:1px solid #d0d4e0; text-align:right;",
        "in")
      cells[[length(cells)+1]] <- div(
        style = "font-size:9px; color:#7880a0; padding:6px 8px 4px;
                 text-transform:uppercase; letter-spacing:0.5px;
                 border-bottom:1px solid #d0d4e0; text-align:right;",
        "out")
      cells[[length(cells)+1]] <- div(
        style = "font-size:9px; color:#7880a0; padding:6px 8px 4px;
                 text-transform:uppercase; letter-spacing:0.5px;
                 border-bottom:1px solid #d0d4e0; text-align:center;",
        "→")
      
      for (i in seq_len(nrow(gd))) {
        gene_name <- gd$gene[i]
        cells[[length(cells)+1]] <- div(
          style = "padding:4px 8px; font-size:11px; color:#3a4060;
                   font-weight:600; border-bottom:1px solid #f0f2f6;
                   align-self:center;",
          gene_name)
        cells[[length(cells)+1]] <- div(
          style = "padding:4px 8px; font-size:11px; color:#1e7a45;
                   font-weight:700; text-align:right;
                   font-family:'Consolas',monospace;
                   border-bottom:1px solid #f0f2f6; align-self:center;",
          fnum(gd$ratio[i], "%.2f"))
        cells[[length(cells)+1]] <- div(
          style = "padding:4px 8px; font-size:11px; color:#5a6080;
                   text-align:right; font-family:'Consolas',monospace;
                   border-bottom:1px solid #f0f2f6; align-self:center;",
          fnum(gd$mean_in[i], "%.2f"))
        cells[[length(cells)+1]] <- div(
          style = "padding:4px 8px; font-size:11px; color:#7880a0;
                   text-align:right; font-family:'Consolas',monospace;
                   border-bottom:1px solid #f0f2f6; align-self:center;",
          fnum(gd$mean_out[i], "%.2f"))
        
        cells[[length(cells)+1]] <- div(
          style = "padding:4px 8px; border-bottom:1px solid #f0f2f6;
                   text-align:center; align-self:center;",
          tags$button(
            "→",
            type = "button",
            class = "btn btn-default btn-fp-icon roi-acquire-btn",
            `data-gene` = gene_name,
            style = "padding:2px 10px; font-size:13px; min-width:auto;
                     line-height:1; cursor:pointer;"
          )
        )
      }
      
      do.call(div, c(
        list(style = "display:grid;
                      grid-template-columns: 1.4fr 0.7fr 0.7fr 0.7fr auto;
                      align-items:stretch;"),
        cells))
    } else {
      div(style = "padding:12px; color:#7880a0; font-size:11px;",
          "No genes passed the enrichment filter (mean expression < 0.5 inside).")
    }
    
    tagList(
      
      div(style = "padding:10px 12px; background:#f5f7fb;
                   border-radius:6px; margin-bottom:12px;
                   display:flex; align-items:center; justify-content:space-between;",
          div(
            div(style = "font-size:12px; font-weight:700; color:#2a3050;",
                "🟦 Region of Interest"),
            div(style = "font-size:11px; color:#7880a0; margin-top:2px;",
                sprintf("%d spots / bins · area %s µm²",
                        d$n_spots,
                        if (is.finite(d$area_um2))
                          format(round(d$area_um2), big.mark = ",")
                        else "—"))
          ),
          actionButton("roi_clear", "✕ Clear ROI",
                       class = "btn-fp-icon",
                       style = "padding:4px 10px; font-size:11px;")
      ),
      
      div(style = "text-align:center; margin-bottom:12px;",
          tags$img(src = img_data,
                   style = "max-width:100%; max-height:200px;
                          border-radius:6px; border:2px solid #c8ccd8;
                          object-fit:contain;")
      ),
      
      div(class = "fp-section", "🧱 Tissue Context"),
      tagList(tissue_rows),
      
      div(class = "fp-section", style = "margin-top:14px;",
          "🧬 Top enriched genes"),
      gene_grid,
      
      div(style = "margin-top:14px; text-align:center;",
          downloadButton("dl_roi_all", "💾 Save ROI as ZIP",
                         class = "btn-fp", style = "width:auto;")
      )
    )
  })
  
  observeEvent(input$roi_acquire_gene, {
    raw <- input$roi_acquire_gene
    if (is.null(raw) || !nzchar(raw)) return()
    gene_name <- sub("__\\d+$", "", raw)
    if (!nzchar(gene_name)) return()
    
    session$sendCustomMessage("openPanel", list(id = "fp-acquire"))
    
    d <- data_list()
    updateSelectizeInput(session, "selected_gene",
                         choices = if (!is.null(d) && !is.null(d$rn)) d$rn else gene_name,
                         selected = gene_name,
                         server = TRUE)
  }, ignoreInit = TRUE)
  
  output$dl_roi_all <- downloadHandler(
    filename = function() {
      sprintf("navis_roi_%s.zip", .safe_name(input$selected_img))
    },
    content = function(zipfile) {
      rs <- roi_state(); req(rs$active, !is.null(rs$data))
      d <- rs$data
      tdir <- tempfile(); dir.create(tdir)
      on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
      
      EBImage::writeImage(d$he_crop, file.path(tdir, "roi_crop.png"))
      write.csv(d$prior_df, file.path(tdir, "roi_tissue_context.csv"), row.names = FALSE)
      write.csv(d$gene_df,  file.path(tdir, "roi_top_genes.csv"),      row.names = FALSE)
      
      writeLines(c(
        "NaVis ROI bundle",
        sprintf("  Channel:   %s", input$selected_img %||% "—"),
        sprintf("  Spots/bins inside: %d", d$n_spots),
        sprintf("  Area:      %s µm²",
                if (is.finite(d$area_um2))
                  format(round(d$area_um2), big.mark = ",") else "—"),
        sprintf("  Corners:   (%.3f, %.3f) → (%.3f, %.3f)",
                d$corners$x1, d$corners$y1, d$corners$x2, d$corners$y2),
        "",
        "Files:",
        "  roi_crop.png            H&E crop of the ROI bounding box",
        "  roi_tissue_context.csv  per-compartment fractions inside the ROI",
        "  roi_top_genes.csv       top-20 genes by enrichment ratio (mean_in / mean_out)"
      ), file.path(tdir, "summary.txt"))
      
      zip::zipr(zipfile, files = list.files(tdir, full.names = TRUE))
    }
  )
  
  observeEvent(roi_state(), {
    rs <- roi_state()
    session$sendCustomMessage("toggleROIDrawer",
                              list(roi = isTRUE(rs$active)))
  })
  
  outputOptions(output, "spot_coords_ui",    suspendWhenHidden = FALSE)
  outputOptions(output, "spot_crop",         suspendWhenHidden = FALSE)
  outputOptions(output, "spot_histogram",    suspendWhenHidden = FALSE)
  outputOptions(output, "spot_expression_ui",suspendWhenHidden = FALSE)
  outputOptions(output, "spot_priors_ui",    suspendWhenHidden = FALSE)
  outputOptions(output, "spot_correlations", suspendWhenHidden = FALSE)
  outputOptions(output, "compare_drawer_ui", suspendWhenHidden = FALSE)
  outputOptions(output, "roi_drawer_ui",     suspendWhenHidden = FALSE)
  
  .gen_meta_image <- function(meta_nm, spot_brush = NULL, super_resolve = FALSE) {
    lst <- data_list(); req(lst)
    pv  <- prms()
    label_prefix <- if (super_resolve) "Super-resolving metadata: " else "Generating metadata map: "
    withProgress(message = paste0(label_prefix, meta_nm), value=0, {
      incProgress(0.3)
      a <- .EBIDW(lst, meta_nm,
                  sigma=pv[1], iterations=pv[2], threshold=pv[3], shp=pv[4],
                  passive=!super_resolve,
                  normalize=FALSE,
                  spot_brush=spot_brush,
                  datatype=if(odt()=="seu") "seu" else "noseu")
      
      a_mat <- if (inherits(a, "Image")) imageData(a) else as.matrix(a)
      if (length(dim(a_mat)) > 2) a_mat <- a_mat[,,1]
      a_mat[!is.finite(a_mat)] <- 0
      
      
      
      if (!super_resolve) {
        rng <- range(a_mat, na.rm = TRUE)
        if (diff(rng) > 1e-12) a_mat <- (a_mat - rng[1]) / diff(rng)
        else                   a_mat <- a_mat * 0
      } else {
        
        a_mat <- pmin(pmax(a_mat, 0), 1)
      }
      a <- EBImage::Image(a_mat, colormode = "grayscale")
      
      
      
      if (super_resolve && pv[1] != 0) {
        z <- if (isTRUE(input$auto_mode))
          EBImage::channel(lst$original_image, "gray")
        else if (!is.null(primg()) && !is.null(imgs()[[primg()]]))
          EBImage::channel(EBImage::readImage(imgs()[[primg()]]), "gray")
        else
          EBImage::channel(lst$original_image, "gray")
        
        zs         <- EBImage::gblur(z, sigma = 4)
        local_var  <- abs(z - zs)
        decay_mask <- exp(-4 * local_var)
        
        
        .conform_to <- function(arr_or_img, ref_img) {
          ref_dims <- dim(ref_img); cur_dims <- dim(arr_or_img)
          if (is.null(cur_dims) || is.null(ref_dims)) return(arr_or_img)
          if (length(cur_dims) >= 2 && length(ref_dims) >= 2 &&
              (cur_dims[1] != ref_dims[1] || cur_dims[2] != ref_dims[2])) {
            if (inherits(arr_or_img, "Image"))
              return(EBImage::resize(arr_or_img, w = ref_dims[1], h = ref_dims[2]))
            tmp <- EBImage::Image(as.matrix(arr_or_img), colormode = "Grayscale")
            tmp <- EBImage::resize(tmp, w = ref_dims[1], h = ref_dims[2])
            return(EBImage::imageData(tmp))
          }
          arr_or_img
        }
        decay_mask <- .conform_to(decay_mask, a)
        a <- a * (1 - decay_mask)
        a <- EBImage::normalize(a) ^ 0.45
      }
      
      incProgress(0.7)
      
      lv <- lst$meta[[meta_nm]]
      if (is.null(lv)) lv <- lst$meta[, meta_nm, drop = TRUE]
      if (!is.null(dim(lv))) lv <- as.vector(lv[, 1, drop = TRUE])
      is_cat <- !is.numeric(lv)
      g   <- col_gene(a, color = if (is_cat) NULL else input$meta_color)
      
      di <- g$image
      di_dim <- dim(di)
      if (length(di_dim) == 3 && di_dim[3] >= 3 &&
          !is.null(lst$white_mask)) {
        wm  <- lst$white_mask
        wm_arr <- if (inherits(wm, "Image")) imageData(wm) else as.matrix(wm)
        if (length(dim(wm_arr)) > 2) wm_arr <- wm_arr[,,1]
        if (all(dim(wm_arr) == di_dim[1:2])) {
          msk <- 1 - wm_arr
          di[,,1] <- di[,,1] * msk
          di[,,2] <- di[,,2] * msk
          di[,,3] <- di[,,3] * msk
          g$image <- di
        }
      }
      
      base_nm  <- if (super_resolve) paste0("META_SR_", meta_nm) else paste0("META_", meta_nm)
      outpath  <- file.path(img_dir, paste0(base_nm, ".png"))
      EBImage::writeImage(di, outpath)
      cur <- imgs(); if (is.null(cur)) cur <- list()
      cur[[base_nm]] <- outpath
      imgs(cur)
      lgs <- img_legends(); lgs[[base_nm]] <- g$legend; img_legends(lgs)
      incProgress(1)
      session$sendCustomMessage("openPanel",     list(id="fp-channels"))
      session$sendCustomMessage("selectChannel", list(key=base_nm))
    })
  }
  
  observeEvent(input$gen_meta, {
    req(data_list(), input$meta_col)
    .gen_meta_image(input$meta_col, spot_brush = NULL)
  })
  
  
  meta_col_is_numeric <- reactive({
    d <- data_list()
    if (is.null(d) || is.null(d$meta) || is.null(input$meta_col)) return(FALSE)
    if (!(input$meta_col %in% colnames(d$meta))) return(FALSE)
    lv <- d$meta[[input$meta_col]]
    if (is.null(lv)) lv <- d$meta[, input$meta_col, drop = TRUE]
    if (!is.null(dim(lv))) lv <- as.vector(lv[, 1, drop = TRUE])
    is.numeric(lv)
  })
  
  
  output$meta_col_typehint <- renderUI({
    d <- data_list()
    if (is.null(d) || is.null(d$meta) || is.null(input$meta_col) ||
        !(input$meta_col %in% colnames(d$meta))) {
      return(tags$div(style = "font-size:10px; color:#9aa0b8;
                              font-style:italic; margin-top:-6px; margin-bottom:6px;",
                      "Pick a column to see its type."))
    }
    lv <- d$meta[[input$meta_col]]
    if (is.null(lv)) lv <- d$meta[, input$meta_col, drop = TRUE]
    if (!is.null(dim(lv))) lv <- as.vector(lv[, 1, drop = TRUE])
    
    is_num <- is.numeric(lv)
    n_total <- length(lv)
    n_na    <- sum(is.na(lv))
    
    if (is_num) {
      vfin <- lv[is.finite(lv)]
      if (length(vfin) > 0) {
        rng <- range(vfin)
        med <- stats::median(vfin)
        n_uniq <- length(unique(vfin))
        body <- sprintf("numeric \u2022 n=%d \u2022 NA=%d \u2022 unique=%d \u2022 range [%.3g, %.3g] \u2022 median %.3g",
                        n_total, n_na, n_uniq, rng[1], rng[2], med)
      } else {
        body <- sprintf("numeric \u2022 all values NA or non-finite (n=%d)", n_total)
      }
      tags$div(style = "font-size:10px; color:#1f8a3c; margin-top:-6px;
                       margin-bottom:6px; font-family:Consolas,monospace;",
               body)
    } else {
      lv_chr <- as.character(lv)
      n_levels <- length(unique(lv_chr[!is.na(lv_chr)]))
      body <- sprintf("categorical \u2022 n=%d \u2022 NA=%d \u2022 levels=%d",
                      n_total, n_na, n_levels)
      tags$div(style = "font-size:10px; color:#a86a00; margin-top:-6px;
                       margin-bottom:6px; font-family:Consolas,monospace;",
               body)
    }
  })
  
  
  output$meta_param_preview <- renderUI({
    if (!isTRUE(meta_col_is_numeric())) {
      return(tags$div(
        style = "font-size:11px; color:#a86a00; padding:6px 8px;
                 background:#fff5e6; border-radius:4px; border:1px solid #f0c98a;",
        tags$b("Disabled."), " The selected column is not numeric.
        Super-resolution applies only to continuous values."
      ))
    }
    
    
    auto_on <- isTRUE(input$auto_mode)
    if (auto_on) {
      sigma   <- if (isTRUE(odt() == "seu")) 5  else 10
      iters   <- if (isTRUE(odt() == "seu")) 1  else 3
      thr     <- if (isTRUE(odt() == "seu")) 0.2 else 0.5
      shp     <- if (isTRUE(odt() == "seu")) 10 else 30
      src_lbl <- "from Acquire auto-mode"
    } else {
      sigma <- if (!is.null(input$sigma))      input$sigma      else 10
      iters <- if (!is.null(input$iterations)) input$iterations else 3
      thr   <- if (!is.null(input$threshold))  input$threshold  else 0.5
      shp   <- if (!is.null(input$sharpness))  input$sharpness  else 30
      src_lbl <- "from Acquire manual settings"
    }
    
    row_fn <- function(label, value) {
      tags$div(
        style = "display:flex; justify-content:space-between; font-size:11px;
                 padding:1px 0;",
        tags$span(style="color:#5c6480;", label),
        tags$span(style="color:#3a4060; font-family:Consolas,monospace;", value)
      )
    }
    
    tags$div(
      style = "font-size:11px; padding:6px 8px;
               background:#eef4ff; border-radius:4px; border:1px solid #c8d5ee;",
      tags$div(style="font-size:10px; color:#5c7;
                     margin-bottom:3px; font-style:italic;",
               sprintf("Parameters %s:", src_lbl)),
      row_fn("Sigma",      sprintf("%g", sigma)),
      row_fn("Iterations", sprintf("%d", as.integer(iters))),
      row_fn("Threshold",  sprintf("%.2f", thr)),
      row_fn("Sharpness",  sprintf("%g", shp))
    )
  })
  
  
  observe({
    if (isTRUE(meta_col_is_numeric())) {
      shinyjs::enable("super_meta")
    } else {
      shinyjs::disable("super_meta")
    }
  })
  
  
  observeEvent(input$super_meta, {
    req(data_list(), input$meta_col)
    if (!isTRUE(meta_col_is_numeric())) {
      showNotification(
        "Super-resolve: the selected column is not numeric. Operation skipped.",
        type = "warning", duration = 4)
      return()
    }
    .gen_meta_image(input$meta_col, spot_brush = NULL, super_resolve = TRUE)
  })
  
  
  outputOptions(output, "meta_col_typehint",  suspendWhenHidden = FALSE)
  outputOptions(output, "meta_param_preview", suspendWhenHidden = FALSE)
  
  observeEvent(input$apply_meta_spot, {
    sel <- input$selected_img
    req(sel, grepl("^META_", sel))
    
    is_sr   <- grepl("^META_SR_", sel)
    meta_nm <- sub("^META_(SR_)?", "", sel)
    .gen_meta_image(meta_nm,
                    spot_brush = input$meta_spot_size,
                    super_resolve = is_sr)
  })
  
  insights_data <- reactiveVal(NULL)
  
  .get_gene_mat <- function(sel, imgs_list) {
    tryCatch({
      p <- imgs_list[[sel]]
      img <- if (is.character(p)) EBImage::readImage(p) else p
      if (all(dim(img)[1:2] < 800)) {
        d <- dim(img); ar <- d[1] / d[2]
        new_w <- if (ar >= 1) 2000 else round(2000 * ar)
        new_h <- if (ar >= 1) round(2000 / ar) else 2000
        img <- EBImage::resize(img, w = new_w, h = new_h)
      }
      mat <- imageData(img)
      if (length(dim(mat)) == 2) mat
      else (mat[,,1] + mat[,,2] + mat[,,3]) / 3
    }, error=function(e) NULL)
  }
  
  .get_prior_mat <- function(nm, lst, target_dim) {
    tryCatch({
      pm  <- lst[[nm]]
      mat <- if (inherits(pm,"Image")) imageData(pm) else pm
      if (length(dim(mat)) > 2) mat <- mat[,,1]
      mat <- as.matrix(mat)
      if (!all(dim(mat) == target_dim)) {
        pimg <- EBImage::Image(mat, colormode="Grayscale")
        
        pimg <- EBImage::resize(pimg, w=target_dim[2], h=target_dim[1])
        mat  <- as.matrix(imageData(pimg))
        
        if (!all(dim(mat) == target_dim) && all(rev(dim(mat)) == target_dim)) {
          mat <- t(mat)
        }
        if (!all(dim(mat) == target_dim)) {
          
          out <- matrix(0, nrow=target_dim[1], ncol=target_dim[2])
          r_copy <- min(nrow(mat), target_dim[1])
          c_copy <- min(ncol(mat), target_dim[2])
          out[1:r_copy, 1:c_copy] <- mat[1:r_copy, 1:c_copy]
          mat <- out
        }
      }
      mat / (max(mat, na.rm=TRUE) + 1e-8)
    }, error=function(e) NULL)
  }
  
  observeEvent(input$run_insights, {
    sel <- input$selected_img
    lst <- isolate(data_list())
    if (is.null(sel) || is.null(lst)) { insights_data(NULL); return() }
    
    is_prior <- sel %in% c("H&E Image","Entropy","Non-fibrillar interstitial tissue",
                           "Fibrillar interstitial tissue","Nuclei and cells")
    is_meta  <- grepl("^META_", sel)
    if (is_prior || is_meta) { insights_data(NULL); return() }
    
    gmat_full <- .get_gene_mat(sel, isolate(imgs()))
    if (is.null(gmat_full)) { insights_data(NULL); return() }
    
    withProgress(message = "Computing spatial insights…", value = 0, {
      incProgress(0.1, detail = "Downsampling")
      
      max_analysis_px <- 400L
      scale_f <- max_analysis_px / max(dim(gmat_full))
      if (scale_f < 1) {
        g_img <- EBImage::Image(gmat_full, colormode="Grayscale")
        g_img <- EBImage::resize(g_img, w=round(dim(gmat_full)[1]*scale_f),
                                 h=round(dim(gmat_full)[2]*scale_f))
        gmat <- as.matrix(imageData(g_img))
      } else {
        gmat <- gmat_full
      }
      
      td   <- dim(gmat)
      comp_nms  <- c("nuclei_mask","fibrous_tissue_mask","soft_tissue_mask")
      comp_lbls <- c("Nuclei","Fibrous ECM","Soft Tissue")
      comp_cols <- c("#4488ff","#ff8844","#44bb88")
      
      gvec  <- as.numeric(gmat)
      gmean <- mean(gvec[gvec > 0], na.rm=TRUE)   
      enrich <- lapply(seq_along(comp_nms), function(i) {
        pm <- .get_prior_mat(comp_nms[i], lst, td)
        if (is.null(pm)) return(list(ratio=NA, p=NA))
        pm_vec <- as.numeric(pm)
        
        nz <- pm_vec[pm_vec > 0]
        if (length(nz) < 10) return(list(ratio=NA, p=NA))
        thresh <- median(nz)
        mask   <- pm_vec >= thresh
        if (sum(mask) < 10) return(list(ratio=NA, p=NA))
        
        g_in   <- gvec[mask]
        g_in   <- g_in[g_in > 0]
        if (length(g_in) < 5) return(list(ratio=1, p=1))
        obs    <- mean(g_in, na.rm=TRUE)
        ratio  <- obs / (gmean + 1e-8)
        
        n_perm <- if (.is_hd_dataset(lst)) 100L else 200L
        n_mask <- sum(mask & gvec > 0)
        g_expr <- gvec[gvec > 0]
        null_r <- replicate(n_perm, {
          idx <- sample.int(length(g_expr), min(n_mask, length(g_expr)))
          mean(g_expr[idx], na.rm=TRUE) / (gmean + 1e-8)
        })
        p <- mean(abs(null_r - 1) >= abs(ratio - 1))
        list(ratio=ratio, p=p)
      })
      enrich_df <- data.frame(
        label = comp_lbls,
        ratio = sapply(enrich, `[[`, "ratio"),
        p     = sapply(enrich, `[[`, "p"),
        color = comp_cols,
        stringsAsFactors = FALSE
      )
      
      kx <- matrix(c(-1,0,1,-2,0,2,-1,0,1), nrow=3)
      ky <- t(kx)
      gx_map <- filter2(gmat, kx)
      gy_map <- filter2(gmat, ky)
      grad   <- sqrt(gx_map^2 + gy_map^2)
      grad   <- EBImage::normalize(grad)
      
      prior_edge_corrs <- sapply(comp_nms, function(nm) {
        pm <- .get_prior_mat(nm, lst, td)
        if (is.null(pm)) return(NA_real_)
        pe <- filter2(pm, kx)^2 + filter2(pm, ky)^2
        pe <- sqrt(pe)
        suppressWarnings(cor(as.numeric(grad), as.numeric(pe), use="complete.obs"))
      })
      boundary_align <- data.frame(label=comp_lbls, corr=as.numeric(prior_edge_corrs),
                                   color=comp_cols, stringsAsFactors=FALSE)
      
      shifts_px <- round(c(0, 25, 50, 100, 200, 400) /
                           max(lst$microns_per_pixel, 0.1))
      shifts_px <- pmax(shifts_px, 0L)
      xcorr_list <- lapply(comp_nms, function(nm) {
        pm <- .get_prior_mat(nm, lst, td)
        if (is.null(pm) || !all(dim(pm) == td))
          return(rep(NA_real_, length(shifts_px)))
        pv2 <- as.numeric(pm)
        sapply(shifts_px, function(s) {
          tryCatch({
            if (s == 0) return(suppressWarnings(cor(gvec, pv2, use="complete.obs")))
            nr <- td[1]; nc <- td[2]
            if (s >= nc) return(NA_real_)
            pshift <- matrix(NA_real_, nr, nc)
            pshift[, (s+1):nc] <- pm[, 1:(nc-s)]
            suppressWarnings(cor(gvec, as.numeric(pshift), use="complete.obs"))
          }, error = function(e) NA_real_)
        })
      })
      xcorr_df <- do.call(rbind, lapply(seq_along(comp_nms), function(i) {
        data.frame(
          label    = comp_lbls[i],
          lag_um   = shifts_px * max(lst$microns_per_pixel, 0.1),
          corr     = xcorr_list[[i]],
          color    = comp_cols[i],
          stringsAsFactors = FALSE
        )
      }))
      
      thr_val  <- quantile(gvec, 0.75, na.rm=TRUE)
      bin_mask <- gmat > thr_val
      labeled  <- EBImage::bwlabel(EBImage::Image(bin_mask+0, colormode="Grayscale"))
      lab_mat  <- imageData(labeled)
      n_labs   <- max(lab_mat)
      island_df <- if (n_labs > 0 && n_labs <= 5000) {
        incProgress(0.75, detail = "Island morphology")
        mpp    <- max(lst$microns_per_pixel, 0.1)
        
        nz_idx <- which(lab_mat > 0, arr.ind=TRUE)
        nz_lab <- lab_mat[lab_mat > 0]
        if (length(nz_lab) > 0) {
          areas <- tabulate(nz_lab, nbins=n_labs)
          
          row_vec <- nz_idx[,1]; col_vec <- nz_idx[,2]
          rmin <- tapply(row_vec, nz_lab, min)
          rmax <- tapply(row_vec, nz_lab, max)
          cmin <- tapply(col_vec, nz_lab, min)
          cmax <- tapply(col_vec, nz_lab, max)
          
          lb_keys <- as.integer(names(areas_keys <- tabulate(nz_lab, nbins=n_labs)))
          hr <- rmax - rmin + 1
          wr <- cmax - cmin + 1
          circs <- 4*pi*as.numeric(areas[as.integer(names(rmin))]) /
            ((2*(hr+wr))^2 + 1e-6)
          
          full_circs <- rep(NA_real_, n_labs)
          full_circs[as.integer(names(rmin))] <- circs
          areas_um2 <- areas * mpp^2
          data.frame(area_um2   = areas_um2,
                     circularity = full_circs,
                     stringsAsFactors = FALSE)
        } else data.frame(area_um2=numeric(0), circularity=numeric(0))
      } else data.frame(area_um2=numeric(0), circularity=numeric(0))
      
      pm_list <- lapply(comp_nms, function(nm) .get_prior_mat(nm, lst, td))
      valid   <- !sapply(pm_list, is.null)
      if (sum(valid) >= 1) {
        X_cols  <- do.call(cbind, lapply(pm_list[valid], as.numeric))
        y_dec   <- gvec
        keep    <- is.finite(y_dec) & apply(X_cols, 1, function(r) all(is.finite(r)))
        X_k     <- X_cols[keep,, drop=FALSE]; y_k <- y_dec[keep]
        ss_tot  <- sum((y_k - mean(y_k))^2)
        fit_dec <- tryCatch(lm.fit(cbind(1, X_k), y_k), error=function(e) NULL)
        r2 <- if (!is.null(fit_dec)) {
          ss_res <- sum(fit_dec$residuals^2)
          1 - ss_res / (ss_tot + 1e-8)
        } else NA_real_
        decouple <- 1 - max(0, min(1, r2))
      } else {
        decouple <- NA_real_
      }
      
      incProgress(1, detail = "Done")
      insights_data(list(
        enrich_df     = enrich_df,
        grad_img      = grad,
        boundary_align= boundary_align,
        xcorr_df      = xcorr_df,
        island_df     = island_df,
        decouple      = decouple,
        n_islands     = nrow(island_df)
      ))
    })   
  })
  
  simgenes_data <- reactiveVal(NULL)
  
  observeEvent(input$run_simgenes, {
    sel <- input$selected_img
    lst <- isolate(data_list())
    if (is.null(sel) || is.null(lst)) return()
    is_prior <- sel %in% c("H&E Image","Entropy","Non-fibrillar interstitial tissue",
                           "Fibrillar interstitial tissue","Nuclei and cells")
    is_meta  <- grepl("^META_", sel)
    is_interface <- grepl("^interface_", sel)
    if (is_prior || is_meta) {
      showNotification("Similar genes are not available for priors or metadata images.",
                       type="warning", duration=4)
      return()
    }
    withProgress(message="Finding similar spatial profiles\u2026", value=0, {
      
      if (is_interface) {
        incProgress(0.1, detail="Sampling interface image at spot locations")
        img_path <- isolate(imgs())[[sel]]
        if (is.null(img_path) || !is.character(img_path) || !file.exists(img_path)) {
          showNotification("Interface image file is missing; recompute the merge.",
                           type="error", duration=5)
          return()
        }
        
        
        img <- tryCatch(EBImage::readImage(img_path),
                        error = function(e) NULL)
        if (is.null(img)) {
          warning("run_simgenes: failed to read interface image at ", img_path)
          showNotification("Could not read interface image.",
                           type="error", duration=5)
          return()
        }
        
        
        if (length(dim(img)) == 3) {
          dat <- EBImage::imageData(img)
          intf <- (dat[,,1] + dat[,,2] + dat[,,3]) / 3
        } else {
          intf <- EBImage::imageData(img)
        }
        ih <- nrow(intf); iw <- ncol(intf)
        
        
        coords <- lst$coords
        if (is.null(coords) || nrow(coords) == 0) {
          showNotification("No spot coordinates available.",
                           type="error", duration=5)
          return()
        }
        cx <- pmax(1L, pmin(iw, as.integer(round(coords$x))))
        cy <- pmax(1L, pmin(ih, as.integer(round(coords$y))))
        
        
        rad_px <- max(2L, as.integer(round(
          (if (!is.null(lst$microns_per_pixel))
            27.5 / lst$microns_per_pixel else 6) / 2
        )))
        q_vec <- vapply(seq_along(cx), function(i) {
          x1 <- max(1L, cx[i] - rad_px); x2 <- min(iw, cx[i] + rad_px)
          y1 <- max(1L, cy[i] - rad_px); y2 <- min(ih, cy[i] + rad_px)
          
          mean(intf[y1:y2, x1:x2], na.rm = TRUE)
        }, numeric(1))
        
        
        bcs_expr <- colnames(lst$expression)
        if (!is.null(coords$barcode) && !is.null(bcs_expr)) {
          ord <- match(bcs_expr, coords$barcode)
          q_vec <- q_vec[ord]
          
          if (any(is.na(q_vec))) {
            q_vec[is.na(q_vec)] <- 0
          }
        }
        
        if (sum(q_vec > 0.001) < 5) {
          showNotification(
            paste0("Interface footprint at spots is too sparse (", sum(q_vec > 0.001),
                   " non-zero spots). Try adjusting reach or threshold."),
            type="warning", duration=6)
          return()
        }
        
        
        qg <- character(0)
        kg <- rownames(lst$expression)
        
      } else {
        incProgress(0.1, detail="Resolving query gene")
        kg <- rownames(lst$expression)
        query_genes <- .resolve_one(sel, kg)
        if (is.null(query_genes) || length(query_genes) == 0) query_genes <- sel
        qg <- query_genes[query_genes %in% kg]
        if (length(qg) == 0) {
          showNotification("Query gene not found in expression matrix.", type="error", duration=5)
          return()
        }
        q_mat <- lst$expression[qg, , drop=FALSE]
        q_vec <- if (nrow(q_mat) == 1) as.numeric(q_mat) else Matrix::colMeans(q_mat)
        
        if (sum(q_vec > 0) < 5) {
          showNotification("Too few expressed spots to compute correlations.", type="warning", duration=5)
          return()
        }
      }
      
      test_genes <- setdiff(kg, qg)
      incProgress(0.3, detail = sprintf("Correlating %d genes over all spots\u2026", length(test_genes)))
      
      X <- as.matrix(lst$expression[test_genes, , drop = FALSE])
      
      nz_per_gene <- rowSums(X > 0)
      X <- X[nz_per_gene >= 10L, , drop = FALSE]
      if (nrow(X) == 0) stop("no genes with sufficient expression")
      
      q_c  <- q_vec - mean(q_vec)
      X_c  <- X - rowMeans(X)
      
      dots    <- as.numeric(X_c %*% q_c)
      X_norms <- sqrt(rowSums(X_c^2))
      q_norm  <- sqrt(sum(q_c^2))
      cors    <- dots / (X_norms * q_norm + 1e-8)
      names(cors) <- rownames(X)
      cors <- cors[is.finite(cors)]
      
      cors <- cors[!names(cors) %in% qg]
      top_idx <- order(cors, decreasing = TRUE)[seq_len(min(10, length(cors)))]
      df <- data.frame(gene = names(cors)[top_idx],
                       r    = round(cors[top_idx], 3),
                       stringsAsFactors = FALSE)
      incProgress(1, detail="Done")
      simgenes_data(df)
    })
    session$sendCustomMessage("showSimgenesContent", list())
  })
  
  observeEvent(input$selected_img, {
    simgenes_data(NULL)
    session$sendCustomMessage("hideSimgenesContent", list())
  }, ignoreNULL=TRUE, priority=-1)   
  
  observeEvent(input$selected_img, {
    insights_data(NULL)
    session$sendCustomMessage("hideEnrichContent", list())
  }, ignoreNULL = TRUE)
  
  output$enrich_target_ui <- renderUI({
    sel <- input$selected_img
    if (is.null(sel) || nchar(sel) == 0) {
      return(div(style="text-align:center;font-size:11px;color:#c85050;
                        background:#fff4f4;border-radius:5px;padding:6px;margin-bottom:8px;",
                 "No channel selected."))
    }
    is_prior <- sel %in% c("H&E Image","Entropy","Non-fibrillar interstitial tissue",
                           "Fibrillar interstitial tissue","Nuclei and cells")
    is_meta  <- grepl("^META_", sel)
    if (is_prior) {
      return(div(style="text-align:center;font-size:11px;color:#c07010;
                        background:#fff8f0;border-radius:5px;padding:6px;margin-bottom:8px;",
                 "Insights only apply to gene / signature / composite channels, not tissue priors."))
    }
    if (is_meta) {
      return(div(style="text-align:center;font-size:11px;color:#7030a0;
                        background:#f5f0fa;border-radius:5px;padding:6px;margin-bottom:8px;",
                 "Insights are not available for metadata images — they apply to
                  gene, signature, or composite expression channels."))
    }
    div(style="text-align:center;font-size:11px;color:#2a5030;
               background:#f0f8f2;border-radius:5px;padding:6px;margin-bottom:8px;",
        "Active channel: ",
        tags$b(style="color:#1e5030;", sel))
  })
  
  observe({
    id <- insights_data()
    if (!is.null(id)) {
      session$sendCustomMessage("showEnrichContent", list())
    }
  })
  
  output$enrich_plot <- renderPlot({
    id <- insights_data(); req(id)
    df <- id$enrich_df; req(!all(is.na(df$ratio)))
    .plot_enrich(df, base_size = 9) +
      theme(plot.background=element_rect(fill="transparent",color=NA),
            panel.background=element_rect(fill="transparent",color=NA),
            axis.text=element_text(color="#3a4060",size=9),
            axis.title=element_text(color="#7880a0",size=8))
  }, bg="transparent")
  
  output$boundary_img <- renderImage({
    id <- insights_data(); req(id)
    tmp <- tempfile(fileext=".png")
    EBImage::writeImage(EBImage::Image(id$grad_img, colormode="Grayscale"), tmp)
    list(src=tmp, contentType="image/png",
         style="max-width:100%; max-height:100%; display:block;
                margin:0 auto; border-radius:4px; border:1px solid #c8ccd8;
                object-fit:contain;")
  }, deleteFile=TRUE)
  
  output$boundary_align_ui <- renderUI({
    id <- insights_data(); req(id)
    df <- id$boundary_align
    tagList(
      tags$p(style="font-size:10px;color:#7880a0;margin:4px 0 2px;",
             "Boundary alignment with tissue edges:"),
      tagList(lapply(seq_len(nrow(df)), function(i) {
        v <- if (is.na(df$corr[i])) 0 else df$corr[i]
        pct <- round(abs(v)*100)
        div(class="prior-row",
            div(class="prior-label", df$label[i]),
            div(class="prior-bar-bg",
                div(class="prior-bar-fill",
                    style=sprintf("width:%d%%;background:%s;", pct, df$color[i]))),
            div(class="prior-value", sprintf("r=%.2f", v))
        )
      }))
    )
  })
  
  output$xcorr_plot <- renderPlot({
    id <- insights_data(); req(id)
    df <- id$xcorr_df
    df <- df[is.finite(df$corr),]
    req(nrow(df) > 0)
    .plot_xcorr(df, base_size = 9) +
      theme(plot.background=element_rect(fill="transparent",color=NA),
            panel.background=element_rect(fill="transparent",color=NA),
            panel.grid=element_line(color="#e0e4ef",linewidth=0.3),
            axis.text=element_text(color="#3a4060",size=8),
            legend.text=element_text(size=8),
            legend.key.height=unit(0.4,"cm"))
  }, bg="transparent")
  
  output$island_plot <- renderPlot({
    id <- insights_data(); req(id)
    df <- id$island_df
    req(nrow(df) > 0)
    df2 <- df[is.finite(df$circularity) & is.finite(df$area_um2),]
    req(nrow(df2) > 0)
    .plot_island(df, base_size = 9) +
      theme(plot.background=element_rect(fill="transparent",color=NA),
            panel.background=element_rect(fill="transparent",color=NA),
            panel.grid=element_line(color="#e0e4ef",linewidth=0.3),
            axis.text=element_text(color="#3a4060",size=8),
            plot.title=element_text(color="#7880a0",size=9))
  }, bg="transparent")
  
  output$island_summary_ui <- renderUI({
    id <- insights_data(); req(id)
    df <- id$island_df
    if (nrow(df) == 0) return(tags$p(style="font-size:11px;color:#7880a0;","No islands above threshold."))
    med_area <- median(df$area_um2, na.rm=TRUE)
    med_circ <- median(df$circularity, na.rm=TRUE)
    interp   <- .island_interpreter(df, isolate(data_list()))
    div(class="spot-card",
        div(class="spot-card-title","Island summary"),
        div(style="font-size:12px;color:#2a3050;line-height:1.8;",
            sprintf("Count: %d", nrow(df)), tags$br(),
            sprintf("Median area: %.0f µm²", med_area), tags$br(),
            sprintf("Median circularity: %.2f", med_circ)),
        
        div(style=sprintf("margin-top:8px; padding-top:6px;
                         border-top:1px solid #eceff4; font-size:11px;
                         line-height:1.5; color:#3a4060;"),
            tags$span(style=sprintf("color:%s; font-weight:600;", interp$color),
                      interp$label),
            tags$br(),
            tags$span(style="color:#5a6080;", interp$desc),
            if (nzchar(interp$hint)) tagList(
              tags$br(),
              tags$span(style="color:#7880a0; font-style:italic;", interp$hint)
            )
        )
    )
  })
  
  .island_interpreter <- function(df, lst = NULL) {
    if (is.null(df) || !nrow(df))
      return(list(label = "No islands",
                  color = "#7880a0",
                  desc  = "Channel has no high-expression regions above threshold.",
                  hint  = ""))
    
    bin_um <- if (!is.null(lst) && isTRUE(.is_hd_dataset(lst)))
      lst$visium_hd_bin_um else 100
    pitch_ratio   <- 100 / max(as.numeric(bin_um), 1)
    k_n           <- max(1, pitch_ratio)
    k_sparse_area <- max(1, pitch_ratio)
    
    a  <- df$area_um2
    cc <- df$circularity[is.finite(df$circularity)]
    n  <- nrow(df)
    A_med <- median(a, na.rm = TRUE)
    A_max <- max(a, na.rm = TRUE)
    A_tot <- sum(a, na.rm = TRUE)
    c_med <- if (length(cc)) median(cc, na.rm = TRUE) else NA_real_
    f_lrg <- if (A_tot > 0) A_max / A_tot else 1
    
    out <- list(label = "Mixed / ambiguous",
                color = "#7880a0",
                desc  = sprintf("%d islands, median area %s µm², median circularity %s.",
                                n,
                                format(round(A_med), big.mark = ","),
                                if (is.finite(c_med)) sprintf("%.2f", c_med) else "—"),
                hint  = "Pattern does not match a clean morphological category.")
    
    if (n <= 3 * k_n && f_lrg >= 0.85) {
      out <- list(
        label = "Diffuse / pan-tissue",
        color = "#9aabcf",
        desc  = "One or a few large regions cover most of the expressing area.",
        hint  = "e.g. broadly active gene, fibrotic bridging, or large necrotic region.")
    } else if (n <= 5 * k_n && f_lrg >= 0.6 && is.finite(c_med) && c_med < 0.5) {
      out <- list(
        label = "Confluent / fused",
        color = "#c87844",
        desc  = "A few large irregular regions dominate, with no fine compartmentalisation.",
        hint  = "e.g. confluent disease, large tumour mass, or coalesced inflammation.")
    } else if (n >= 8 * k_n && is.finite(c_med) && c_med >= 0.6 && f_lrg < 0.4) {
      out <- list(
        label = "Compartmental, regular",
        color = "#1e7a45",
        desc  = sprintf("%d discrete round regions of consistent shape.", n),
        hint  = "e.g. glandular units (ducts, acini, crypts), follicles, tertiary lymphoid structures, well-circumscribed nodules.")
    } else if (n >= 8 * k_n && is.finite(c_med) && c_med < 0.5 && f_lrg < 0.4) {
      out <- list(
        label = "Compartmental, irregular",
        color = "#7a5c1e",
        desc  = sprintf("%d discrete regions with irregular boundaries.", n),
        hint  = "e.g. inflammatory infiltrates, infiltrative tumour growth, ill-defined immune aggregates.")
    } else if (n >= 5 * k_n && is.finite(c_med) && c_med < 0.35 && f_lrg < 0.5) {
      out <- list(
        label = "Linear / vascular",
        color = "#2c5be0",
        desc  = sprintf("%d elongated regions of low circularity.", n),
        hint  = "e.g. vasculature, tubular structures in cross-section, fibre tracts.")
    } else if (n >= 15 * k_n && A_med <= 200 * k_sparse_area && f_lrg < 0.2) {
      out <- list(
        label = "Sparse / scattered",
        color = "#7a4ab0",
        desc  = sprintf("%d small regions distributed across the tissue.", n),
        hint  = "e.g. rare positive cells, salt-and-pepper marker, isolated lineage-specific cells.")
    }
    out
  }
  
  .decouple_band <- function(v) {
    if (is.null(v) || !is.finite(v)) {
      return(list(band="—", color="#7880a0",
                  msg="Decoupling could not be computed (insufficient finite values)."))
    }
    if (v < 0.40)
      list(band="LOW", color="#1e7a45",
           msg="Gene spatial pattern is largely explained by tissue architecture.")
    else if (v < 0.80)
      list(band="MODERATE", color="#d07010",
           msg="Gene spatial pattern is partly independent of tissue architecture.")
    else
      list(band="HIGH", color="#e03060",
           msg="Gene spatial pattern is largely independent of tissue architecture —
                a candidate for novel biology not reducible to histology.")
  }
  
  output$decouple_ui <- renderUI({
    id <- insights_data(); req(id)
    v  <- id$decouple
    if (is.na(v)) return(tags$p(style="font-size:11px;color:#7880a0;","Insufficient tissue priors."))
    b  <- .decouple_band(v)
    pct <- round(v * 100)
    div(
      div(class="decouple-score", style=paste0("color:",b$color,";"), paste0(pct,"%")),
      div(class="decouple-label", "Histology Decoupling Score (1−R²)"),
      div(style=paste0("font-size:11px;color:",b$color,";text-align:center;padding:4px;"), b$msg)
    )
  })
  
  .plot_enrich <- function(df, base_size = 9) {
    df$ratio[is.na(df$ratio)] <- 1
    df$sig <- ifelse(df$p < 0.01, "★★", ifelse(df$p < 0.05, "★", ""))
    df$label2 <- paste0(df$label, " ", df$sig)
    df$label2 <- factor(df$label2, levels=rev(df$label2))
    df$ratio_c <- pmin(pmax(df$ratio, 0), 3)
    ggplot(df, aes(x=label2, y=ratio_c, fill=color)) +
      geom_col(width=0.55) +
      geom_hline(yintercept=1, linetype="dashed", color="#888888", linewidth=0.5) +
      scale_fill_identity() +
      scale_y_continuous(limits=c(0,3), breaks=c(0,1,2,3),
                         labels=c("0","1×\n(expected)","2×","3×")) +
      coord_flip() +
      labs(x=NULL, y="Enrichment ratio", title=NULL) +
      theme_minimal(base_size=base_size) +
      theme(panel.grid.major.x=element_line(color="#e0e4ef",linewidth=0.3),
            panel.grid.major.y=element_blank())
  }
  .plot_xcorr <- function(df, base_size = 9) {
    df <- df[is.finite(df$corr),]
    df$label <- factor(df$label, levels=unique(df$label))
    ggplot(df, aes(x=lag_um, y=corr, color=color, group=label)) +
      geom_line(linewidth=1) + geom_point(size=1.5) +
      geom_hline(yintercept=0, linetype="dashed", color="#cccccc", linewidth=0.4) +
      scale_color_identity(guide="legend", labels=unique(df$label)) +
      scale_x_continuous(breaks=unique(round(df$lag_um))) +
      scale_y_continuous(limits=c(-1,1)) +
      labs(x="Spatial lag (µm)", y="Cross-correlation", color=NULL) +
      theme_minimal(base_size=base_size) +
      theme(legend.position="bottom")
  }
  .plot_island <- function(df, base_size = 9) {
    df2 <- df[is.finite(df$circularity) & is.finite(df$area_um2),]
    ggplot(df2, aes(x=area_um2, y=circularity)) +
      geom_point(alpha=0.55, size=1.5, color="#2c5be0") +
      scale_x_log10() +
      labs(x="Island area (µm²)", y="Circularity",
           title=paste0(nrow(df), " expression islands")) +
      theme_minimal(base_size=base_size)
  }
  
  zoom_target <- reactiveVal(NULL)
  showZoomModal <- function(target, title) {
    zoom_target(target)
    body <- switch(target,
                   enrich   = plotOutput("zoom_plot", height = "550px"),
                   xcorr    = plotOutput("zoom_plot", height = "550px"),
                   boundary = tagList(
                     div(style="background:#000;border-radius:6px;padding:8px;margin-bottom:16px;
                   display:flex;align-items:center;justify-content:center;min-height:380px;",
                         imageOutput("zoom_image", height="400px", inline=TRUE, width="auto")
                     ),
                     div(style="padding:12px 16px;background:#fafbfc;border-radius:6px;
                   border:1px solid #e0e4ef;",
                         tags$h4("Boundary alignment with tissue compartments",
                                 style="margin:0 0 4px 0;color:#3a4060;font-size:14px;"),
                         tags$p(style="font-size:11px;color:#7880a0;margin:0 0 10px;",
                                "Pearson correlation between gene gradient and each compartment's
                  gradient. Higher absolute values = gene boundaries align with
                  that tissue architecture."),
                         uiOutput("zoom_boundary_align")
                     )
                   ),
                   island = tagList(
                     plotOutput("zoom_plot", height = "420px"),
                     div(style="margin-top:16px;padding:12px 16px;background:#fafbfc;
                   border-radius:6px;border:1px solid #e0e4ef;",
                         tags$h4("Island summary", style="margin:0 0 4px 0;color:#3a4060;font-size:14px;"),
                         uiOutput("zoom_island_summary")
                     )
                   ),
                   decouple = uiOutput("zoom_decouple_display"),
                   NULL
    )
    showModal(modalDialog(
      title = title, size = "l", easyClose = TRUE,
      footer = modalButton("Close"), body
    ))
  }
  observeEvent(input$zoom_enrich,   { showZoomModal("enrich",   "Compartment Enrichment") })
  observeEvent(input$zoom_boundary, { showZoomModal("boundary", "Expression Boundaries — Gradient Map & Alignment") })
  observeEvent(input$zoom_xcorr,    { showZoomModal("xcorr",    "Spatial Cross-correlation") })
  observeEvent(input$zoom_island,   { showZoomModal("island",   "Expression Islands — Distribution & Summary") })
  observeEvent(input$zoom_decouple, { showZoomModal("decouple", "Histology Decoupling Score") })
  
  output$zoom_plot <- renderPlot({
    tgt <- zoom_target(); req(tgt)
    id  <- insights_data(); req(id)
    gene_lbl <- input$selected_img
    base_plot <- switch(tgt,
                        enrich = .plot_enrich(id$enrich_df, base_size = 14),
                        xcorr  = .plot_xcorr (id$xcorr_df,  base_size = 14),
                        island = .plot_island(id$island_df, base_size = 14),
                        NULL
    )
    req(!is.null(base_plot))
    base_plot + labs(subtitle = paste0("Channel: ", gene_lbl))
  }, bg = "white", res = 100)
  
  output$zoom_image <- renderImage({
    req(zoom_target() == "boundary")
    id <- insights_data(); req(id)
    tmp <- tempfile(fileext = ".png")
    EBImage::writeImage(EBImage::Image(id$grad_img, colormode = "Grayscale"), tmp)
    list(src = tmp, contentType = "image/png",
         style = "max-width:100%; max-height:400px; display:block; margin:0 auto;
                  border-radius:6px; border:1px solid #444;
                  object-fit:contain;")
  }, deleteFile = TRUE)
  
  output$zoom_boundary_align <- renderUI({
    id <- insights_data(); req(id)
    df <- id$boundary_align
    tagList(lapply(seq_len(nrow(df)), function(i) {
      v   <- if (is.na(df$corr[i])) 0 else df$corr[i]
      pct <- round(abs(v)*100)
      div(style="display:flex; align-items:center; gap:12px; padding:8px 4px;
                 border-bottom:1px solid #eceff4;",
          div(style=sprintf("width:140px; font-size:13px; font-weight:600; color:%s;",
                            df$color[i]), df$label[i]),
          div(style="flex:1; background:#e8ecf5; border-radius:6px; height:22px;
                   overflow:hidden; position:relative;",
              div(style=sprintf("width:%d%%; height:100%%; background:%s;
                             border-radius:6px;", pct, df$color[i]))
          ),
          div(style="width:90px; text-align:right; font-size:13px; font-variant-numeric:tabular-nums;
                   color:#3a4060;",
              sprintf("r = %.3f", v))
      )
    }))
  })
  
  output$zoom_island_summary <- renderUI({
    id <- insights_data(); req(id)
    df <- id$island_df
    if (is.null(df) || nrow(df) == 0) {
      return(tags$p(style="font-size:13px;color:#7880a0;",
                    "No expression islands detected above the 75th percentile threshold."))
    }
    a  <- df$area_um2
    c  <- df$circularity[is.finite(df$circularity)]
    stat_row <- function(lbl, val)
      div(style="display:flex; justify-content:space-between; padding:4px 0;
                 border-bottom:1px solid #eceff4;",
          div(style="font-size:12px;color:#7880a0;", lbl),
          div(style="font-size:13px;font-weight:600;color:#3a4060;
                     font-variant-numeric:tabular-nums;", val))
    interp <- .island_interpreter(df, isolate(data_list()))
    tagList(
      stat_row("Number of islands",       format(nrow(df), big.mark=",")),
      stat_row("Median area (µm²)",        sprintf("%s", format(round(median(a)), big.mark=","))),
      stat_row("Mean area (µm²)",          sprintf("%s", format(round(mean(a)),   big.mark=","))),
      stat_row("Largest island (µm²)",     sprintf("%s", format(round(max(a)),    big.mark=","))),
      stat_row("Smallest island (µm²)",    sprintf("%s", format(round(min(a)),    big.mark=","))),
      stat_row("Median circularity",       if (length(c)) sprintf("%.3f", median(c)) else "—"),
      stat_row("Mean circularity",         if (length(c)) sprintf("%.3f", mean(c))   else "—"),
      stat_row("Total expressing area (µm²)", format(round(sum(a)), big.mark=",")),
      
      div(style=sprintf("margin-top:14px; padding:12px;
                         border-radius:6px; background:%s22;
                         border-left:3px solid %s;",
                        interp$color, interp$color),
          div(style=sprintf("font-size:13px; font-weight:600; color:%s; margin-bottom:4px;",
                            interp$color),
              interp$label),
          div(style="font-size:12px; color:#3a4060; line-height:1.55;",
              interp$desc),
          if (nzchar(interp$hint))
            div(style="font-size:11px; color:#7880a0; font-style:italic;
                       margin-top:6px; line-height:1.5;",
                interp$hint)
      )
    )
  })
  
  output$zoom_decouple_display <- renderUI({
    id <- insights_data(); req(id)
    v <- id$decouple
    gene_lbl <- input$selected_img
    if (is.null(v) || !is.finite(v)) {
      return(div(style="padding:40px;text-align:center;color:#7880a0;",
                 "Decoupling could not be computed (insufficient finite values)."))
    }
    b <- .decouple_band(v)
    div(style="padding:30px 20px; text-align:center;",
        tags$h3(style="margin:0 0 6px 0;color:#7880a0;font-weight:500;font-size:15px;",
                paste0("Channel: ", gene_lbl)),
        tags$h4(style="margin:0 0 30px 0;color:#7880a0;font-weight:400;font-size:13px;",
                "Histology Decoupling Score (1 − R²)"),
        
        div(style=sprintf("font-size:96px; font-weight:700; font-variant-numeric:tabular-nums;
                         color:%s; line-height:1;", b$color),
            sprintf("%.2f", v)),
        div(style=sprintf("display:inline-block; margin:16px 0 4px; padding:6px 16px;
                         background:%s; color:#fff; border-radius:20px; font-weight:600;
                         font-size:13px; letter-spacing:1px;", b$color),
            paste0(b$band, " DECOUPLING")),
        div(style="max-width:520px; margin:20px auto 12px; font-size:14px;
                 color:#3a4060; line-height:1.5;",
            b$msg),
        
        div(style="max-width:400px; margin:24px auto 12px;",
            div(style="height:10px; background:linear-gradient(to right,
                   #1e7a45 0%, #1e7a45 40%,
                   #d07010 40%, #d07010 80%,
                   #e03060 80%, #e03060 100%);
                   border-radius:5px; position:relative;",
                
                div(style=sprintf("position:absolute; left:%.1f%%; top:-6px; transform:translateX(-50%%);
                             width:0; height:0; border-left:7px solid transparent;
                             border-right:7px solid transparent; border-bottom:10px solid #222;",
                                  min(max(v,0),1)*100))
            ),
            div(style="display:flex; justify-content:space-between; margin-top:4px;
                   font-size:11px; color:#7880a0;",
                span("0.0 — Architecture-bound"),
                span("1.0 — Architecture-free"))
        )
    )
  })
  
  .safe_name <- function(x) gsub("[^A-Za-z0-9_+.-]", "_", if (is.null(x)) "channel" else x)
  
  output$dl_enrich <- downloadHandler(
    filename = function() sprintf("navis_%s_enrichment.csv", .safe_name(input$selected_img)),
    content  = function(file) {
      id <- insights_data(); req(id)
      out <- id$enrich_df[, c("label","ratio","p","color")]
      write.csv(out, file, row.names = FALSE)
    }
  )
  output$dl_boundary <- downloadHandler(
    filename = function() sprintf("navis_%s_boundary_gradient.png", .safe_name(input$selected_img)),
    content  = function(file) {
      id <- insights_data(); req(id)
      EBImage::writeImage(EBImage::Image(id$grad_img, colormode = "Grayscale"), file)
    }
  )
  output$dl_xcorr <- downloadHandler(
    filename = function() sprintf("navis_%s_cross_correlation.csv", .safe_name(input$selected_img)),
    content  = function(file) {
      id <- insights_data(); req(id)
      write.csv(id$xcorr_df, file, row.names = FALSE)
    }
  )
  output$dl_island <- downloadHandler(
    filename = function() sprintf("navis_%s_islands.csv", .safe_name(input$selected_img)),
    content  = function(file) {
      id <- insights_data(); req(id)
      write.csv(id$island_df, file, row.names = FALSE)
    }
  )
  output$dl_all_insights <- downloadHandler(
    filename = function() sprintf("navis_%s_insights.zip", .safe_name(input$selected_img)),
    content  = function(zipfile) {
      id <- insights_data(); req(id)
      tmpdir <- tempfile("insights_"); dir.create(tmpdir)
      on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
      write.csv(id$enrich_df,     file.path(tmpdir, "compartment_enrichment.csv"),   row.names=FALSE)
      write.csv(id$boundary_align,file.path(tmpdir, "boundary_tissue_alignment.csv"),row.names=FALSE)
      write.csv(id$xcorr_df,      file.path(tmpdir, "spatial_cross_correlation.csv"),row.names=FALSE)
      write.csv(id$island_df,     file.path(tmpdir, "expression_islands.csv"),       row.names=FALSE)
      writeLines({
        interp <- .island_interpreter(id$island_df, isolate(data_list()))
        c(sprintf("Channel: %s", input$selected_img),
          sprintf("Decoupling score (1 - R^2): %.4f", id$decouple),
          sprintf("Number of islands: %d", id$n_islands),
          sprintf("Island pattern: %s", interp$label),
          sprintf("  %s", interp$desc),
          if (nzchar(interp$hint)) sprintf("  %s", interp$hint) else "")
      }, file.path(tmpdir, "summary.txt"))
      EBImage::writeImage(EBImage::Image(id$grad_img, colormode="Grayscale"),
                          file.path(tmpdir, "boundary_gradient.png"))
      files_to_zip <- list.files(tmpdir, full.names = TRUE)
      
      zip::zipr(zipfile, files = files_to_zip, recurse = FALSE)
    }
  )
  
  outputOptions(output, "enrich_plot",      suspendWhenHidden=FALSE)
  outputOptions(output, "boundary_img",     suspendWhenHidden=FALSE)
  outputOptions(output, "boundary_align_ui",suspendWhenHidden=FALSE)
  outputOptions(output, "xcorr_plot",       suspendWhenHidden=FALSE)
  outputOptions(output, "island_plot",      suspendWhenHidden=FALSE)
  outputOptions(output, "island_summary_ui",suspendWhenHidden=FALSE)
  outputOptions(output, "decouple_ui",      suspendWhenHidden=FALSE)
  outputOptions(output, "enrich_target_ui", suspendWhenHidden=FALSE)
  
  output$simgenes_ui <- renderUI({
    df <- simgenes_data(); req(df)
    if (nrow(df) == 0)
      return(tags$p(style="font-size:11px;color:#7880a0;",
                    "Could not compute similar genes — check that the query gene is in the expression matrix."))
    tagList(lapply(seq_len(nrow(df)), function(i) {
      pct <- round(max(df$r[i], 0) * 100)
      div(style="display:flex;align-items:center;gap:8px;padding:4px 0;
                 border-bottom:1px solid #f0f2f7;",
          
          div(style=sprintf("width:18px;height:18px;border-radius:50%%;
                           background:%s;color:#fff;font-size:9px;font-weight:700;
                           display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                            if(i<=3) "#2c5be0" else "#9aabcf"),
              i),
          
          tags$span(df$gene[i],
                    style="font-size:11px;font-weight:600;color:#2c3050;
                         flex:1;cursor:pointer;text-decoration:underline dotted #9aabcf;",
                    onclick=sprintf("Shiny.setInputValue('sim_gene_click','%s',{priority:'event'})",
                                    df$gene[i]),
                    title="Click to load in Acquire panel"),
          
          div(style="width:80px;background:#e8ecf5;border-radius:3px;height:8px;overflow:hidden;",
              div(style=sprintf("width:%d%%;height:100%%;background:%s;border-radius:3px;",
                                pct,
                                if(i<=3) "#2c5be0" else "#9aabcf"))),
          div(style="width:36px;text-align:right;font-size:10px;color:#7880a0;
                   font-variant-numeric:tabular-nums;flex-shrink:0;",
              sprintf("%.3f", df$r[i]))
      )
    }))
  })
  outputOptions(output, "simgenes_ui", suspendWhenHidden=FALSE)
  
  observeEvent(input$sim_gene_click, {
    gene <- input$sim_gene_click
    req(nchar(gene) > 0)
    updateSelectizeInput(session, "selected_gene",
                         choices  = isolate(data_list())$rn,
                         selected = gene,
                         server   = TRUE)
    session$sendCustomMessage("openPanel", list(id="fp-acquire"))
  })
  
  observeEvent(input$zoom_simgenes, {
    zoom_target("simgenes")
    showModal(modalDialog(
      title = "Similar Spatial Profiles — Top 10 genes",
      size = "m", easyClose = TRUE, footer = modalButton("Close"),
      uiOutput("zoom_simgenes_display")
    ))
  })
  output$zoom_simgenes_display <- renderUI({
    df <- simgenes_data(); req(df)
    gene_lbl <- input$selected_img
    if (nrow(df) == 0)
      return(div(style="padding:30px;text-align:center;color:#7880a0;",
                 "No similar genes computed."))
    tagList(
      tags$h4(style="margin:0 0 4px 0;font-size:14px;color:#3a4060;",
              paste0("Query: ", gene_lbl)),
      tags$p(style="font-size:12px;color:#7880a0;margin:0 0 16px;",
             "Pearson r on spot-level expression (query-expressed spots only).
              Click a gene name to load it in the Acquire panel."),
      lapply(seq_len(nrow(df)), function(i) {
        pct <- round(max(df$r[i], 0) * 100)
        div(style="display:flex;align-items:center;gap:12px;padding:10px 4px;
                   border-bottom:1px solid #eceff4;",
            div(style=sprintf("width:24px;height:24px;border-radius:50%%;
                             background:%s;color:#fff;font-size:11px;font-weight:700;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                              if(i<=3) "#2c5be0" else "#9aabcf"), i),
            tags$span(df$gene[i],
                      style="width:90px;font-size:13px;font-weight:600;color:#2c3050;
                           cursor:pointer;text-decoration:underline dotted #9aabcf;",
                      onclick=sprintf(
                        "Shiny.setInputValue('sim_gene_click','%s',{priority:'event'});
                       $('.modal').modal('hide');", df$gene[i]),
                      title="Click to load in Acquire panel"),
            div(style="flex:1;background:#e8ecf5;border-radius:4px;height:12px;overflow:hidden;",
                div(style=sprintf("width:%d%%;height:100%%;background:%s;border-radius:4px;",
                                  pct, if(i<=3) "#2c5be0" else "#9aabcf"))),
            div(style="width:60px;text-align:right;font-size:13px;color:#3a4060;
                     font-variant-numeric:tabular-nums;font-weight:600;",
                sprintf("r = %.3f", df$r[i]))
        )
      })
    )
  })
  
  output$dl_simgenes <- downloadHandler(
    filename = function() {
      nm <- if (!is.null(input$selected_img) && nchar(input$selected_img) > 0)
        input$selected_img else "channel"
      sprintf("navis_%s_similar_genes.csv", gsub("[^A-Za-z0-9_+.-]","_", nm))
    },
    content = function(file) {
      df <- simgenes_data(); req(df)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  spot_zoom_target <- reactiveVal(NULL)
  
  showSpotZoomModal <- function(target, title) {
    spot_zoom_target(target)
    body <- switch(target,
                   crop = div(style="background:#fafbfc;border-radius:6px;padding:12px;
                        display:flex;align-items:center;justify-content:center;
                        min-height:500px;",
                              imageOutput("zoom_spot_image", height="500px",
                                          inline=TRUE, width="auto")),
                   hist = plotOutput("zoom_spot_plot", height="500px"),
                   corr = plotOutput("zoom_spot_plot", height="500px"),
                   NULL
    )
    showModal(modalDialog(
      title = title, size = "l", easyClose = TRUE,
      footer = modalButton("Close"), body
    ))
  }
  observeEvent(input$zoom_spot_crop, { showSpotZoomModal("crop", "H&E Crop at Clicked Spot") })
  observeEvent(input$zoom_spot_hist, { showSpotZoomModal("hist", "H&E Colour Histogram") })
  observeEvent(input$zoom_spot_corr, { showSpotZoomModal("corr", "Channel vs Tissue Correlations") })
  
  output$zoom_spot_image <- renderImage({
    req(spot_zoom_target() == "crop")
    sd <- spot_data(); req(sd)
    tmp <- tempfile(fileext = ".png")
    EBImage::writeImage(sd$he_crop, tmp)
    list(src = tmp, contentType = "image/png",
         style = "max-width:100%; max-height:500px; display:block; margin:0 auto;
                  border-radius:6px; border:1px solid #c8ccd8; object-fit:contain;")
  }, deleteFile = TRUE)
  
  output$zoom_spot_plot <- renderPlot({
    tgt <- spot_zoom_target(); req(tgt)
    sd  <- spot_data(); req(sd)
    if (tgt == "hist") {
      .plot_spot_histogram(sd$hist_df, base_size = 14)
    } else if (tgt == "corr") {
      .plot_spot_correlations(sd$corr_df, base_size = 14)
    }
  }, bg = "white", res = 100)
  
  .spot_suffix <- function() {
    sd <- spot_data()
    if (is.null(sd)) return("spot")
    sprintf("spot_%s_%s",
            .safe_name(input$selected_img),
            if (!is.null(sd$spot_name)) .safe_name(sd$spot_name) else "location")
  }
  
  output$dl_spot_crop <- downloadHandler(
    filename = function() sprintf("navis_%s_crop.png", .spot_suffix()),
    content  = function(file) {
      sd <- spot_data(); req(sd)
      EBImage::writeImage(sd$he_crop, file)
    }
  )
  
  output$dl_spot_hist <- downloadHandler(
    filename = function() sprintf("navis_%s_histogram.csv", .spot_suffix()),
    content  = function(file) {
      sd <- spot_data(); req(sd)
      write.csv(sd$hist_df, file, row.names = FALSE)
    }
  )
  
  output$dl_spot_context <- downloadHandler(
    filename = function() sprintf("navis_%s_tissue_context.csv", .spot_suffix()),
    content  = function(file) {
      sd <- spot_data(); req(sd)
      
      out <- if (!is.null(sd$prior_df)) sd$prior_df else
        data.frame(label = character(), percent = numeric())
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  output$dl_spot_corr <- downloadHandler(
    filename = function() sprintf("navis_%s_correlations.csv", .spot_suffix()),
    content  = function(file) {
      sd <- spot_data(); req(sd)
      write.csv(sd$corr_df, file, row.names = FALSE)
    }
  )
  
  output$dl_spot_all <- downloadHandler(
    filename = function() sprintf("navis_%s_bundle.zip", .spot_suffix()),
    content  = function(zipfile) {
      sd <- spot_data(); req(sd)
      tdir <- tempfile(); dir.create(tdir)
      on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
      suffix <- .spot_suffix()
      f_crop <- file.path(tdir, sprintf("%s_crop.png",         suffix))
      f_hist <- file.path(tdir, sprintf("%s_histogram.csv",    suffix))
      f_ctx  <- file.path(tdir, sprintf("%s_tissue_context.csv", suffix))
      f_corr <- file.path(tdir, sprintf("%s_correlations.csv", suffix))
      f_sum  <- file.path(tdir, sprintf("%s_summary.txt",      suffix))
      
      EBImage::writeImage(sd$he_crop, f_crop)
      write.csv(sd$hist_df, f_hist, row.names = FALSE)
      write.csv(
        if (!is.null(sd$prior_df)) sd$prior_df
        else data.frame(label = character(), percent = numeric()),
        f_ctx, row.names = FALSE)
      write.csv(sd$corr_df, f_corr, row.names = FALSE)
      
      writeLines(c(
        sprintf("NaVis spot analysis bundle"),
        sprintf("  Channel:   %s", input$selected_img),
        sprintf("  Spot:      %s", if (!is.null(sd$spot_name)) sd$spot_name else "(custom point)"),
        sprintf("  Gene value @ spot: %s",
                if (is.null(sd$gene_val) || is.na(sd$gene_val)) "—"
                else sprintf("%.4f", sd$gene_val)),
        "",
        "Files in this bundle:",
        "  *_crop.png             H&E image crop centred on the clicked spot",
        "  *_histogram.csv        RGB intensity densities of the crop",
        "  *_tissue_context.csv   % composition across the three tissue priors",
        "  *_correlations.csv     Pearson r of channel image vs each tissue prior"
      ), f_sum)
      
      zip::zipr(zipfile, files = c(f_crop, f_hist, f_ctx, f_corr, f_sum))
    }
  )
  
  session$onSessionEnded(function() {
    unlink(img_dir, recursive = TRUE, force = TRUE)
  })
}

shinyApp(ui, server)