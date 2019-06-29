library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Heebo")
showtext_auto()## Automatically use showtext to render text for future devices

## use the ggplot2 example
library(hexSticker)
img_url <- normalizePath("~/mapper/man/figures/mapper_logo.png")
sticker(img_url, package="Mapper", p_size=9, 
        s_x = 1, s_y = 0.70, s_width = 0.6, s_height = 0.4,
        p_color = rgb(240/256, 240/256, 240/256), 
        h_size = 2, h_color = rgb(240/256, 240/256, 240/256), h_fill = "red",
        spotlight = FALSE, 
        dpi = 380, 
        filename="~/mapper/hexbin/mapper_sticker.png")