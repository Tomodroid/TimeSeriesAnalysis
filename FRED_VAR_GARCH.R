# FREDの時系列データ（GDPやPPIなど）を用いたVAR・GARCH分析

library(tidyverse)
library(vars)
library(fredr)

fredr_set_key("724bc62bb11871a434fc640df13d7a7c")

GDPC1 <- fredr_series_observations(
 series_id = "GDPC1",
 observation_start = as.Date("1959-04-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "pch" 
)

GDPDEF <- fredr_series_observations(
 series_id = "GDPDEF",
 observation_start = as.Date("1959-04-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "pch"
)

PPIACO <- fredr_series_observations(
 series_id = "PPIACO",
 observation_start = as.Date("1959-04-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "pch" 
)

FEDFUNDS <- fredr_series_observations(
 series_id = "FEDFUNDS",
 observation_start = as.Date("1959-04-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "lin"
)

M2SL <- fredr_series_observations(
 series_id = "M2SL",
 observation_start = as.Date("1959-04-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "pch"
)

df <- data.frame(
 Time = as.Date(GDPC1$date, "%Y/%m/%d"),
 GDPC1 = GDPC1$value,
 GDPDEF = GDPDEF$value,
 PPIACO = PPIACO$value,
 FEDFUNDS = FEDFUNDS$value,
 M2SL = M2SL$value
)

g1 <- ggplot(df, aes(x = Time, y = GDPC1)) +
 geom_line(col="navy") +
 xlab("Time") + ylab("GDPC1") +
 theme_bw()
g1
# ggsave("GDPC1.png", device = "png", width = 15, height = 5, units = "cm")

g2 <- ggplot(df, aes(x = Time, y = GDPDEF)) +
 geom_line(col="red") +
 xlab("Time") + ylab("GDPDEF") +
 theme_bw()
g2
# ggsave("GDPDEF.png", device = "png", width = 15, height = 5, units = "cm")

g3 <- ggplot(df, aes(x = Time, y = PPIACO)) +
 geom_line(col="darkgreen") +
 xlab("Time") + ylab("PPIACO") +
 theme_bw()
g3
# ggsave("PPIACO.png", device = "png", width = 15, height = 5, units = "cm")

g4 <- ggplot(df, aes(x = Time, y = FEDFUNDS)) +
 geom_line(col="orange") +
 xlab("Time") + ylab("FEDFUNDS") +
 theme_bw()
g4
# ggsave("FEDFUNDS.png", device = "png", width = 15, height = 5, units = "cm")

g5 <- ggplot(df, aes(x = Time, y = M2SL)) +
 geom_line(col="purple") +
 xlab("Time") + ylab("M2SL") +
 theme_bw()
g5
# ggsave("M2SL.png", device = "png", width = 15, height = 5, units = "cm")

VARselect(df[, 2:6], lag.max = 20, type = "const")

VARModel <- VAR(df[, 2:6], p = 3)
summary(VARModel)

bmat <- diag(5)
bmat[lower.tri(bmat)] <- NA

SVARmodel <- SVAR(
 VARModel,
 Amat = NULL,
 Bmat = bmat,
 estmethod = c("scoring", "direct"),
 hessian = TRUE,
 max.iter = 100  # 1000
)

irfs <- irf(
 SVARmodel,
 impulse = "FEDFUNDS",
 response = c("GDPC1", "GDPDEF", "PPIACO", "FEDFUNDS", "M2SL"),
 ortho = TRUE,
 cumulative = FALSE,
 n.ahead = 12,
 runs = 100  # 1000
)

# png("irfs.png", width = 2000, height = 3000, res = 200, pointsize = 18)
par(mfrow = c(3, 2))
plot(irfs, plot.type = "single", main = "")
par(mfrow = c(1, 1))
# dev.off()


GDPC1 <- fredr_series_observations(
 series_id = "GDPC1",
 observation_start = as.Date("1948-01-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "pch" 
)

UNRATE <- fredr_series_observations(
 series_id = "UNRATE",
 observation_start = as.Date("1948-01-01"),
 observation_end = as.Date("2019-12-31"),
 frequency = "q",
 units = "lin" 
)

df2 <- data.frame(
 GDPC1 = GDPC1$value,
 UNRATE = UNRATE$value
)

VARselect(df2, lag.max = 20, type = "const")
VARModel2 <- VAR(df2, p = 3)
summary(VARModel2)

BQmodel <- BQ(VARModel2)
summary(BQmodel)

BQirf <- irf(
 BQmodel,
 impulse = "GDPC1",
 response = c("GDPC1", "UNRATE"),
 ortho = TRUE,
 cumulative = FALSE,
 n.ahead = 12,
 runs = 1000  # 1000
)

# png("BQirf.png", width = 2000, height = 1000, res = 200, pointsize = 12)
par(mfrow = c(1, 2))
plot(BQirf, plot.type = "single", main = "")
par(mfrow = c(1, 1))
# dev.off()

library(rugarch)

SP500 <- fredr_series_observations(
 series_id = "SP500",
 observation_start = as.Date("2018-01-06"),
 observation_end = as.Date("2023-01-06"),
 frequency = "d",
 units = "lin" 
)

df3 <- SP500 |>
 mutate(Time = as.Date(date, "%Y/%m/%d")) |>
 drop_na() |>
 mutate(Return = (log(value) - log(lag(value))) * 100) |>
 drop_na()

g7 <- ggplot(df3, aes(x = Time, y = Return)) +
 geom_line(col="navy") +
 xlab("Time") + ylab("Daily return %") +
 theme_bw()
g7

GARCH1 <- ugarchspec(
 variance.model = list(model = "sGARCH",
                       garchOrder = c(1, 1)),
 mean.model = list(armaOrder = c(0, 0), 
                   include.mean = TRUE),
 distribution.model = "norm"
)

GARCH2 <- ugarchspec(
 variance.model = list(model = "sGARCH",
                       garchOrder = c(1, 1)),
 mean.model = list(armaOrder = c(0, 0), 
                   include.mean = TRUE),
 distribution.model = "std"
)

EGARCH1 <- ugarchspec(
 variance.model = list(model = "eGARCH",
                       garchOrder = c(1, 1)),
 mean.model = list(armaOrder = c(0, 0), 
                   include.mean = TRUE),
 distribution.model = "norm"
)

EGARCH2 <- ugarchspec(
 variance.model = list(model = "eGARCH",
                       garchOrder = c(1, 1)),
 mean.model = list(armaOrder = c(0, 0), 
                   include.mean = TRUE),
 distribution.model = "std"
)

GARCHmodel1 <- ugarchfit(
 spec = GARCH1,
 data = df3$Return
)

GARCHmodel2 <- ugarchfit(
 spec = GARCH2,
 data = df3$Return
)

EGARCHmodel1 <- ugarchfit(
 spec = EGARCH1,
 data = df3$Return
)

EGARCHmodel2 <- ugarchfit(
 spec = EGARCH2,
 data = df3$Return
)

infocriteria(GARCHmodel1)
infocriteria(EGARCHmodel1)
infocriteria(GARCHmodel2)
infocriteria(EGARCHmodel2)
