.libPaths("D:/R/library")
suppressMessages({library(pdftools); library(magick)})
im <- image_read(pdf_render_page("C:/Users/lenovo/Desktop/Science2023_人卵黄囊细胞图谱_Science.pdf",
                                 page = 37, dpi = 300))
geo <- sprintf("%dx%d+%d+%d", round(2550*0.345), round(3300*0.30),
               round(2550*0.205), round(3300*0.425))
d <- image_crop(im, geo)
image_write(d, "fig2d/results/fig2d_panelD.png", format = "png")
cat("cropped:", image_info(d)$width, "x", image_info(d)$height, "\n")
