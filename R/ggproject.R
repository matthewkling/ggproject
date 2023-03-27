
add_limits <- function(prj, dlim){

        if(class(prj) != "proj3d") stop("`prj` must be a `proj3d` object but is not.")

        if(is.null(prj$xlim)) prj$xlim <- range(dlim[, 1])
        if(is.null(prj$ylim)) prj$ylim <- range(dlim[, 2])
        if(is.null(prj$zlim)) prj$zlim <- range(dlim[, 3])

        breaks <- prj[c("xlim", "ylim", "zlim")] %>%
                lapply(function(x) labeling::extended(x[1], x[2], 5, only.loose = T) )
        prj$xbreaks <- breaks$xlim
        prj$ybreaks <- breaks$ylim
        prj$zbreaks <- breaks$zlim

        prj
}

#' Define a projection
#'
#' This function
#' @export
#'
#' @param shear
#' @param yaw degrees C
#' @param pitch degrees C
#' @param roll degrees C
#' @param persp degrees C
#' @param dist degrees C
#' @param hjust degrees C
#' @param vjust degrees C
#' @param xlim degrees C
#' @param ylim degrees C
#' @param zlim degrees C
#' @param xbreaks degrees C
#' @param ybreaks degrees C
#' @param zbreaks degrees C
#' @param dlim A three-column dataset from which to calculate axis limits.
#' @return A `proj3d` projection object
projection <- function(shear = 0, yaw = 0, pitch = 0, roll = 0,
                       persp = F, dist = 1, hjust = .5, vjust = .5,
                       xlim = NULL, ylim = NULL, zlim = NULL,
                       xbreaks = NULL, ybreaks = NULL, zbreaks = NULL,
                       dlim = NULL){

        # placeholder for parameter checks and messages

        prj <- list(shear = shear, yaw = yaw, pitch = pitch, roll = roll,
                    persp = persp, dist = dist, hjust = hjust, vjust = vjust,
                    xlim = xlim, ylim = ylim, zlim = zlim)
        class(prj) <- "proj3d"
        if(!is.null(dlim)) prj <- add_limits(prj, dlim)
        prj
}

project <- function(data, scales, prj = projection(), expand = 0) {

        # rescale data to unit cube with origin at midpoint
        if(is.null(prj$xlim)) prj$xlim <- range(data[, "x"])
        if(is.null(prj$ylim)) prj$ylim <- range(data[, "y"])
        if(is.null(prj$zlim)) prj$zlim <- range(data[, "z"])
        rscl <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)) - .5
        data[, "x"] <- scales::rescale(data[, "x"], c(-.5, .5), prj$xlim)
        data[, "y"] <- scales::rescale(data[, "y"], c(-.5, .5), prj$ylim)
        data[, "z"] <- scales::rescale(data[, "z"], c(-.5, .5), prj$zlim)
        # data[,c("x", "y", "z")] <- apply(data[,c("x", "y", "z")], 2, rscl)
        if(expand != 0) data[, c("x", "y", "z")] <- data[, c("x", "y", "z")] * (1 + expand)

        # shear
        shr <- diag(3)
        shr[2, 1] <- prj$shear

        # rotate
        yaw <- prj$yaw / 180 * pi
        pitch <- prj$roll / 180 * pi
        roll <- prj$pitch / 180 * pi
        rot <- matrix(c(cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll),
                        sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll),
                        -sin(pitch),         cos(pitch)*sin(roll),                             cos(pitch)*cos(roll)),
                      nrow = 3, byrow = T)

        y <- as.matrix(data[,c("x", "y", "z")]) %*% shr %*% rot

        # perspective
        if(prj$persp){
                p1 <- matrix(0, 4, 4)
                n <- prj$dist
                f <- prj$dist + 1
                lr <- .5
                tb <- .5
                p1[1, 1] <- n / lr
                p1[2, 2] <- n / tb
                p1[3, 3] <- (f + n) / (n - f)
                p1[3, 4] <- 2 * f * n / (n - f)
                p1[4, 3] <- -1
                y[,3] <- y[,3] - .5 - prj$dist
                y[,2] <- y[,2] + .5 - prj$vjust
                y[,1] <- y[,1] + .5 - prj$hjust
                y <- cbind(y, 1) %*% p1
                y <- apply(y, 1, function(x) x / x[4]) %>% t()
        }

        data[,c("x", "y", "z")] <- y[,1:3]
        data
}



StatProject <- ggproto("StatProject", Stat,
                       compute_panel = function(data, scales, prj = projection()) {
                               project(data, scales, prj)
                       },
                       required_aes = c("x", "y", "z")
)


proj_data <- function(mapping = NULL, data = NULL,
                      geom = "polygon", position = "identity",
                      na.rm = FALSE, show.legend = NA,
                      inherit.aes = TRUE,
                      ...) {
        layer(
                stat = StatProject,
                data = data,
                mapping = mapping,
                geom = geom,
                position = position,
                show.legend = show.legend,
                inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}



StatProjectMargin <- ggproto("StatProjectMargin", Stat,
                             compute_panel = function(data, scales, prj = projection(),
                                                      facets = c("xmin", "ymin", "zmin")) {

                                     d <- map_df(facets, function(f){
                                             if(f == "xmin") return(mutate(data, x = min(prj$xbreaks)))
                                             if(f == "xmax") return(mutate(data, x = max(prj$xbreaks)))
                                             if(f == "ymin") return(mutate(data, y = min(prj$ybreaks)))
                                             if(f == "ymax") return(mutate(data, y = max(prj$ybreaks)))
                                             if(f == "zmin") return(mutate(data, z = min(prj$zbreaks)))
                                             if(f == "zmax") return(mutate(data, z = max(prj$zbreaks)))
                                     })
                                     project(d, scales, prj)

                             },
                             required_aes = c("x", "y", "z")
)

proj_margin <- function(mapping = NULL, data = NULL,
                        geom = "point", position = "identity",
                        na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE,
                        ...) {
        layer(
                stat = StatProjectMargin,
                data = data,
                mapping = mapping,
                geom = geom,
                position = position,
                show.legend = show.legend,
                inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}




gridlines <- function(prj, facets = "auto"){

        breaks <- expand_grid(x = prj$xbreaks,
                              y = prj$ybreaks,
                              z = prj$zbreaks)

        x <- bind_rows(filter(breaks, x == min(x), y %in% range(y) | z %in% range(z)),
                       filter(breaks, x == max(x), y %in% range(y) | z %in% range(z))) %>%
                mutate(piece = paste("x", y, z, sep = "_"))
        y <- bind_rows(filter(breaks, y == min(y), x %in% range(x) | z %in% range(z)),
                       filter(breaks, y == max(y), x %in% range(x) | z %in% range(z))) %>%
                mutate(piece = paste("y", x, z, sep = "_"))
        z <- bind_rows(filter(breaks, z == min(z), y %in% range(y) | x %in% range(x)),
                       filter(breaks, z == max(z), y %in% range(y) | x %in% range(x))) %>%
                mutate(piece = paste("z", x, y, sep = "_"))
        g <- as.data.frame(bind_rows(x, y, z)) %>%
                mutate(axis = str_sub(piece, 1, 1))

        # convert to facets
        g <- bind_rows(g %>% filter(x == min(x)) %>% mutate(facet = "xmin"),
                       g %>% filter(x == max(x)) %>% mutate(facet = "xmax"),
                       g %>% filter(y == min(y)) %>% mutate(facet = "ymin"),
                       g %>% filter(y == max(y)) %>% mutate(facet = "ymax"),
                       g %>% filter(z == min(z)) %>% mutate(facet = "zmin"),
                       g %>% filter(z == max(z)) %>% mutate(facet = "zmax"))

        if(facets[1] == "auto"){
                # placeholder for code to identify background facets
                facets <- c("xmin", "ymin", "zmin")
        }
        if(facets[1] != "all") g <- g %>% filter(facet %in% facets)

        g %>%
                group_by(piece, facet) %>%
                mutate(n = n()) %>%
                filter(n > 1) %>%
                dplyr::select(-n) %>%
                ungroup() %>%
                as.data.frame()

        # determine whether outside of facet faces viewer
        # currently does not work for perspective projections
        # normals <- data.frame(x = c(prj$xlim, rep(mean(prj$xlim), 2), rep(mean(prj$xlim), 2)),
        #                       y = c(rep(mean(prj$ylim), 2), prj$ylim, rep(mean(prj$ylim), 2)),
        #                       z = c(rep(mean(prj$zlim), 2), rep(mean(prj$zlim), 2), prj$zlim),
        #                       facet = c("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"))
        # np <- project(normals, NA, prj) %>%
        #         rename(nx = x, ny = y, nz = z)
        # view <- c(0, 0, 100)
        # g <- left_join(g, np) %>%
        #         mutate(visible = x*nx + y*ny + z*nz <= view[1]*nx + view[2]*ny + view[3]*nz)

}

StatProjectGridlines <- ggproto("StatProjectGridlines", Stat,
                                compute_panel = function(data, scales, prj, facets = "auto", expand = 0) {
                                        if(is.null(prj$xlim) | is.null(prj$ylim) | is.null(prj$zlim)){
                                                prj <- add_limits(prj, data[,c("x", "y", "z")])
                                        }
                                        g <- gridlines(prj, facets)
                                        g$group <- paste(g$facet, g$piece)
                                        project(g, scales, prj, expand)
                                },
                                required_aes = c("x", "y", "z")
)

proj_gridlines <- function(mapping = NULL, data = NULL,
                           position = "identity",
                           na.rm = FALSE, show.legend = NA,
                           inherit.aes = TRUE,
                           ...) {
        layer(
                stat = StatProjectGridlines,
                data = data,
                mapping = mapping,
                geom = "path",
                position = position,
                show.legend = show.legend,
                inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}



StatProjectLabels <-
        ggproto("StatProjectLabels", Stat,
                compute_panel = function(data, scales, prj, expand = 0,
                                         facets = "all",
                                         edges = list(c("ymin", "zmax"), c("xmin", "zmax"), c("ymin", "xmax")),
                                         xtitle = "x", ytitle = "y", ztitle = "z"
                ) {
                        if(is.null(prj$xlim) | is.null(prj$ylim) | is.null(prj$zlim)){
                                prj <- add_limits(prj, data[,c("x", "y", "z")])
                        }

                        g <- gridlines(prj, facets)
                        gp <- project(g, scales, prj, expand)

                        # axis text
                        make_labs <- function(edge = c("ymin", "zmax")){
                                ax <- case_when(!any(str_detect(edge, "xm")) ~ "x",
                                                !any(str_detect(edge, "ym")) ~ "y",
                                                !any(str_detect(edge, "zm")) ~ "z")
                                xx <- g %>% filter(axis != ax, facet == edge[2])
                                xx <- g %>% filter(axis != ax,
                                                   facet == edge[1],
                                                   paste(x, y, z) %in% paste(xx$x, xx$y, xx$z))
                                xx <- gp %>%
                                        bind_cols(g %>% dplyr::select(xr = x, yr = y, zr = z)) %>%
                                        filter(piece %in% xx$piece,
                                               facet %in% edge) %>%
                                        mutate(lab = xr %in% unique(xx$x) &
                                                       yr %in% unique(xx$y) &
                                                       zr %in% unique(xx$z)) %>%
                                        group_by(piece) %>%
                                        arrange(lab) %>%
                                        mutate(angle = atan2(diff(y), diff(x)) / pi * 180,
                                               angle = ifelse(between(angle, -90, 90), angle, angle + 180),
                                               label = case_when(ax == "x" ~ xr,
                                                                 ax == "y" ~ yr,
                                                                 ax == "z" ~ zr),
                                               label = paste0(" ", label, " ")) %>%
                                        slice(2) %>%
                                        ungroup() %>%
                                        mutate(hjust = ifelse(mean(x) > 0, 0, 1),
                                               lab_axis = ax)
                                xx
                        }
                        text <- map_dfr(edges, make_labs)

                        # axis titles
                        # browser()
                        text <-
                                text %>%
                                group_by(lab_axis) %>%
                                summarize(angle = atan2(diff(y[1:2]), diff(x[1:2])) / pi * 180,
                                          length = sqrt(diff(range(x))^2 + diff(range(y))^2),
                                          x = median(x),
                                          y = median(y),
                                          # y = y + sin((90-angle) / 180 * pi) * .15 * length, # shift outward
                                          # x = x - cos((90-angle) / 180 * pi) * .15 * length,
                                          angle = ifelse(between(angle, -90, 90), angle, angle + 180)
                                ) %>%
                                mutate(label = c(xtitle, ytitle, ztitle)) %>%
                                bind_rows(text)
                        text
                },
                required_aes = c("x", "y", "z")
        )

proj_labels <- function(mapping = NULL, data = NULL,
                        position = "identity",
                        na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE,
                        ...) {
        layer(
                stat = StatProjectLabels,
                data = data,
                mapping = mapping,
                geom = "text",
                position = position,
                show.legend = show.legend,
                inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}

