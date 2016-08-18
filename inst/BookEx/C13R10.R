library(FRAPO)
library(fGarch)
library(fMultivar)
library(sn)
library(fPortfolio)
library(PerformanceAnalytics)
## Preparing FX-data / returns
data(ESCBFX)
FXDaily <- timeSeries(ESCBFX, charvec = rownames(ESCBFX))
DDates <- time(FXDaily)
WedDays <- isWeekday(DDates, wday = 3)
FirstWed <- head(which(WedDays, arr.ind = TRUE), 1)
LastWed <- tail(which(WedDays, arr.ind = TRUE), 1)
AllWedDays <- timeSequence(from = DDates[FirstWed],
                           to = DDates[LastWed],
                           by = "week")
DumWed <- timeSeries(rep(1, length(AllWedDays)),
                     charvec = AllWedDays)
FXWeekly <- interpNA(cbind(DumWed, FXDaily),
                     method = "before")[AllWedDays, -1]
assetsNames <- Anames <- FxNames <- colnames(FXWeekly)
FxPrice <- 1 / FXWeekly
FxRet <- returns(FxPrice, percentage = TRUE,
                 type = "discrete", trim = TRUE)
FxRetSub <- window(FxRet,
                   start = start(FxRet), end = time(FxRet)[520])
## Setting parameters / initializing objects
J <- 1000
N <- ncol(FxPrice)
Eperiods <- time(FxRet)[-c(1:519)]
Speriods <- time(FxRet)[1:length(Eperiods)]
LengthBack <- length(Speriods)
WEP <- WMD <- WND <- matrix(NA, nrow = LengthBack, ncol = N)
x0 <- rep(0, N + 1)
pprior <- matrix(rep(1 / J, J), ncol = 1)

