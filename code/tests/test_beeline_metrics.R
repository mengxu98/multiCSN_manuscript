source("code/functions/functions.R")
truth <- data.frame(regulator = c("a", "b", "c"), target = c("b", "c", "a"))
short <- data.frame(regulator = "a", target = "b", weight = -5)
metric <- fast_standard_metric_summary(short, truth)
stopifnot(metric$EarlyPrecision == 1, metric$EPR == 2)

tied <- data.frame(regulator = c("a", "b", "a", "c"),
                   target = c("b", "c", "c", "b"), weight = c(5, 4, 4, 4))
metric <- fast_standard_metric_summary(tied, truth)
stopifnot(metric$EarlyPrecision == 0.5, metric$EPR == 1)
reordered <- fast_standard_metric_summary(tied[c(1, 4, 3, 2), ], truth)
stopifnot(identical(metric, reordered))

invalid <- rbind(short, data.frame(regulator = c("x", "a", "a", "b"),
                                  target = c("y", "a", "b", "c"),
                                  weight = c(100, 100, 2, Inf)))
stopifnot(identical(fast_standard_metric_summary(invalid, truth),
                    fast_standard_metric_summary(short, truth)))
zero <- short; zero$weight <- 0
stopifnot(fast_standard_metric_summary(zero, truth)$EarlyPrecision == 0,
          fast_standard_metric_summary(short[0, ], truth)$EPR == 0)

# A mathematical Youden tie is resolved at the highest score threshold.
pred <- data.frame(regulator = c("a", "a", "b", "b", "c", "c"),
                   target = c("b", "c", "c", "a", "a", "b"),
                   weight = 6:1)
metric <- fast_standard_metric_summary(pred, truth)
stopifnot(metric$Precision == 1, metric$Recall == 1/3)
empty <- data.frame(regulator = character(), target = character(), weight = numeric())
all_positive <- data.frame(regulator = c("A", "B"), target = c("B", "A"))
for (pred in list(empty, transform(all_positive, weight = 0))) {
  metric <- fast_standard_metric_summary(pred, all_positive)
  stopifnot(all(unlist(metric[c("Precision", "Recall", "F1", "ACC", "JI", "SI")]) == 0))
}
truth <- data.frame(regulator = c("a\rb", "a"), target = c("c", "b\rc"))
metric <- fast_standard_metric_summary(transform(truth[1, ], weight = 5), truth)
stopifnot(metric$JI == 0.5, metric$EPR == 6, metric$GoldEdges == 2)
cat("BEELINE metric regression tests passed\n")
