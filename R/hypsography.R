#' generate_hypsography
#' @importFrom stats coefficients lm loess predict
#' @export
#' @param nlayers integer number of layers
#' @param max_depth numeric max depth
#' @param surface_area numeric surface area
#' @param layer_thickness numeric optional overrides nlayers#'
#' @examples
#' nlayers <- 14
#' max_depth <- 18.288000
#' surface_area <- 637641.6
#' res <- generate_hypsography(nlayers, max_depth, surface_area)
#' plot(res[,c("area", "depth_linear")],
#'     ylim = c(max(res$depth_linear), min(res$depth_linear)))
#' points(res[,c("area", "depth_ellipsoid")], col = "red", pch = 19)
#'
generate_hypsography <- function(nlayers, max_depth, surface_area, layer_thickness = NA){

  estimate_linear <- function(nlayers, surface_area){
    layer_thickness <- max_depth / nlayers
    depth_vec <- seq(0, max_depth, by = layer_thickness)

    df_pred <- data.frame(
      areas = c(surface_area, 0),
      layer_ind = c(1, nlayers + 1))

    fit <- lm(areas ~ layer_ind, data = df_pred)
    area_vec <- predict(fit,
                  newdata = data.frame(
                  layer_ind = c(seq_len(nlayers), nlayers + 1)))

    data.frame(cbind(area_vec, depth_vec))
  }

  hypso_linear <- estimate_linear(nlayers = nlayers,
                                  surface_area = surface_area)

  estimate_ellipsoid <- function(x, y, nlayers,
                                 knot_points, max_deform_per, deform_level){

    generate_perpendicular <- function(x, y, knot_points, max_deform_per){
      lapply(knot_points, function(i) {
        slope <- coefficients(lm(y ~ x))[2]
        intercept <- coefficients(lm(y ~ x))[1]

        center_point_y <- slope * -1 * x
        shifted_x <- x + (i * max(x)) - mean(x)
        cross_point_x <- (i * max(x)) + ((mean(x) - (i * max(x))) / 2)
        cross_point_y <- cross_point_x * slope + intercept

        res <- data.frame(cbind(shifted_x, center_point_y))
        names(res) <- c("x", "y")

        res <- res[res$x > cross_point_x & res$y > cross_point_y,]
        res <- res[
          (nrow(res) - floor(max_deform_per * length(x)) + 1):nrow(res),]

        res
      })
    }

    perp_lines <- generate_perpendicular(x = hypso_linear$area_vec,
                                         y = hypso_linear$depth_vec,
                                         knot_points = knot_points,
                                         max_deform_per = max_deform_per)

    n_shortest <- nrow(perp_lines[[
        which.min(lapply(perp_lines, function(x) min(length(x))))]])

    fit_ellipsoid <- function(x, perp_lines, pred_length, deform_level, nlayers, n_shortest){
      perp_lines <- lapply(seq_len(n_shortest), function(i)
        rbind(cbind(x, y)[c(length(y), 1),],
              perp_lines[[1]][i,],
              perp_lines[[2]][i,])[
                c(1, 3, 4, 2),])

      new_dt <- data.frame(seq(min(x), max(x), length.out = length(x)))
      names(new_dt) <- c("x")
      fit <- loess(y ~ x, data = data.frame(perp_lines[[deform_level]]), span = 0.8)
      cbind(new_dt, predict(fit, newdata = new_dt))
    }

    res <- fit_ellipsoid(x = x, perp_lines = perp_lines,
                         nlayers = nlayers,
                         deform_level = deform_level, n_shortest = n_shortest)
    res
  }

 res <- estimate_ellipsoid(hypso_linear$area_vec, hypso_linear$depth_vec,
                     knot_points = c(0.1, 0.7), max_deform_per = 0.25,
                     deform_level = 2, nlayers = nlayers)
 res <- cbind(hypso_linear, rev(res[,2]))
 names(res) <- c("area", "depth_linear", "depth_ellipsoid")
 res
}
