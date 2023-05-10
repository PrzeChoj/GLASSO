if (!require("huge", quietly = TRUE)) {
  install.packages("huge")
}
library(huge)
data(stockdata)
x <- log(stockdata$data[2:1258, ] / stockdata$data[1:1257, ])
colnames(x) <- stockdata$info[, 1]
x
