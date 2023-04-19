library(devtools)
library(usethis)


usethis::use_tidy_description()
devtools::document()
devtools::install()
devtools::check(cran = TRUE)
devtools::check_win_devel()
devtools::check_rhub()

devtools::release()

citation(package = 'PSAgraphics')

# Convert to roxygen2 and Github for maintaining the package
# Rd2roxygen::Rd2roxygen(getwd())


##### Hex Logo #################################################################
library(hexSticker)
library(PSAgraphics)

data("lindner")
lindner.ps <- glm(abcix ~ stent + height + female +	diabetic + acutemi + ejecfrac + ves1proc,
				  data = lindner,
				  family = binomial(link = 'logit'))

hexSticker::sticker(~loess.psa(log(lindner$cardbill),
							   lindner$abcix,
							   lindner.ps$fitted,
							   int = c(.37, .56, .87, 1),
							   cex = 0.1,
							   pch = c(1,1),
							   lines = TRUE,
							   rg = FALSE,
							   xlab = '',
							   ylab = ''),
					filename = 'man/figures/PSAgraphics.png',
					p_size = 16,
					package = 'PSAgraphics',
					url = "jbryer.github.io/PSAgraphics/",
					u_size = 5.0,
					s_width = 2.75, s_height = 2.75,
					s_x = 1, s_y = 1,
					p_x = 1, p_y = 1.5,
					p_color = "#00008B",
					h_fill = '#FFFFFF',
					h_color = '#00008B',
					u_color = 'dark green',
					white_around_sticker = TRUE)

