palettes <- list(
BrBG =     rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",
                 "#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e",
                 "#003c30")),
PiYG =     rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef",
                 "#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221",
                "#276419")),
PRGn =     rev(c("#40004b","#762a83","#9970ab","#c2a5cf","#e7d4e8",
                 "#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837",
                "#00441b")),
PuOr =     rev(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6",
                 "#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788",
                 "#2d004b")),
RdBu =     rev(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
                 "#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac",
                 "#053061")),
RdYlBu =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695")),
BuOr =     c("#002bff","#1a66ff","#3399ff","#66CCff","#99eeff",
            "#ccffff","#ffffcc","#ffee99","#ffee66","#ff9933",
            "#ff661a","#ff2b00"),
BuOrRd =   c("#085aff","#3377ff","#5991ff","#8cb2ff","#bfd4FF",
             "#e6eeff","#f7faff","#ffffcc","#ffff99","#ffff00",
             "#ffcc00","#ff9900","#ff6600","#ff0000"),
DkRdBu =   c("#2a0bd9","#264eff","#40a1ff","#73daff","#abf8ff",
             "#e0ffff","#ffffbf","#ffe099","#ffad73","#f76e5e",
             "#d92632","#a60021"),
BuDkOr =   c("#1f8f99","#52c4cc","#99faff","#b2fcff","#ccfeff",
             "#e6ffff","#ffe6cc","#ffca99","#ffad66","#ff8f33",
             "#cc5800","#994000"),
GnPu =     c("#005100","#008600","#00bc00","#00f100","#51ff51",
             "#86ff86","#bcffbc","#ffffff","#fff1ff","#ffbcff",
             "#ff86ff","#ff51ff","#f100f1","#bc00bc","#860086",
             "#510051"),
RdYlGn =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee08b",
                 "#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850",
                 "#006837")),
Spectral = rev(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b",
             "#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd",
             "#5e4fa2")),
grayscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",
        "#252525",
        "#000000"),
revgray = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",
        "#252525",
        "#000000")),
greyscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",
        "#252525",
        "#000000"),
revgrey = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",
        "#252525",
        "#000000"))
);
.colorEnv <- new.env();
assign("palettes", palettes, env = .colorEnv);


# for backwards compatibility

##' @title Rich color palette
##'
##' @description Deprecated function. Please use
##'     \code{\link[gplots]{rich.colors}} instead.
##'
##' @param n The number of desired colors.
##'
##' @seealso \code{\link[gplots]{rich.colors}}
##' @export
richColors <- function (n) {
    cat("NOTE: function `richColors` is deprecated. Please use `rich.colors` instead.\n")
    return(gplots::rich.colors(n))
}

