library(magick)
# import figures (tiffs) and convert them to have cmyk colorspace, and then write to file

# import
fig1 <- image_read("figs/fig_schisto.tiff")
fig1_cmyk <- image_convert(fig1, colorspace = "cmyk")
image_write(fig1_cmyk, path = "figs/fig_schisto_cmyk.tiff")

fig2 <- image_read("figs/fig_cam.tiff")
fig2_cmyk <- image_convert(fig2, colorspace = "cmyk")
image_write(fig2_cmyk, path = "figs/fig_cam_cmyk.tiff")
