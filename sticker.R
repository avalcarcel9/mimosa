library(hexSticker)

setwd("~/Box/mimosa")

imgurl = "icon2.png"
sticker(imgurl, package="mimosa",p_size=8,
        p_color="slategray",
        s_x=1,s_y=.9,
        h_fill = "orange",
        h_color = "slategray",
        s_width=.33,s_height=.1,
        filename="sticker.png")