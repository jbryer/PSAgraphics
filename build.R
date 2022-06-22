library(devtools)
library(usethis)


devtools::document()
usethis::use_tidy_description()
devtools::install()
devtools::check()


# Convert to roxygen2 and Github for maintaining the package
# Rd2roxygen::Rd2roxygen(getwd())


##### Hex Logo #################################################################
library(hexSticker)
library(PSAgraphics)

data("lindner")

# p <- "man/figures/speed_icon.png"
?loess.psa
lindner.ps <- glm(abcix ~ stent + height + female +	diabetic + acutemi + ejecfrac + ves1proc,
				  data = lindner,
				  family = binomial(link = 'logit'))

loess.psa(log(lindner$cardbill),
		  lindner$abcix,
		  lindner.ps$fitted,
		  int = c(.37, .56, .87, 1),
		  lines = TRUE,
		  xlab = '',
		  ylab =)

hexSticker::sticker(~loess.psa(log(lindner$cardbill),
							   lindner$abcix,
							   lindner.ps$fitted,
							   int = c(.37, .56, .87, 1),
							   lines = TRUE),
					filename = 'man/figures/PSAgraphics.png',
					p_size = 8,
					package = 'PSAgraphics',
					url = "github.com/jbryer/PSAgraphics",
					u_size = 2.5,
					s_width = .75, s_height = .75,
					s_x = 1, s_y = 1,
					p_x = 1, p_y = .60,
					p_color = "#00008B",
					h_fill = '#FFFFFF',
					h_color = '#00008B',
					u_color = '#fa9fb5',
					white_around_sticker = FALSE)

