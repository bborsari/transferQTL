summary(m[m$is_eQTL == "y" & m$prediction == "y", "H3K27me3"])
summary(m[m$is_eQTL == "n" & m$prediction == "y", "sum"])
summary(m[m$is_eQTL == "y" & m$prediction == "y", "sum"])
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 & m$sum > 0), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 & m$sum > 4), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 & m$sum > 1), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
#m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$cv > 0.8 & m$sum > 4), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
#m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$sum > 4), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
#m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse(((m$cv < 0.8) | (m$tss_distance < 10000 & m$sum > 1 & m$cv > 5)), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
#m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse(((m$tss_distance < 10000 & m$sum > 1 & m$cv > 5) | (m$cv < 0.8)), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
#m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse(((m$tss_distance < 10000 & m$sum > 1 & m$cv > 5) | (m$cv < 0.8)), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$sum > 4), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | (m$cv > 5 & m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$tss_distance < 10000 & m$sum > 1 & m$cv > 5), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$sum > 4), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$tss_distance < 10000 & m$cv > 5 & ((m$sum == 0 &) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
m$prediction <- ifelse((m$tss_distance < 10000 & m$cv > 5 & ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5 | m$tss_distance < 10000), ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5 | m$tss_distance < 10000, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$tss_distance < 10000 & m$sum > 1), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$tss_distance < 10000 & ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5 | m$tss_distance < 10000, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1) | (m$sum > 3)), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
View(m[m$is_eQTL == "y" & m$prediction == "n", ])
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 | abs(m$slope) > 0.5), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1) | (m$sum > 3)), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & ((m$sum == 0) | (m$sum < 2 & m$H3K27me3 == 1))), "n", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 | abs(m$slope) > 0.5), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1) | (m$sum > 3)), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(((m$tss_distance < 10000 | m$cv > 5) & ((m$sum == 0))), "n", m$prediction) # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 | abs(m$slope) > 0.5), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1) | (m$sum > 3)), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((m$cv > 5) & (m$sum == 0), "n", m$prediction) # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse((m$cv < 0.8 | abs(m$slope) > 0.5), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse((((m$tss_distance < 10000 | m$cv > 5) & m$sum > 1) | (m$sum > 3)), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m$prediction <- ifelse(m$sum == 0, "n", m$prediction) # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
m[i, "prediction"] <- "y"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
if (m[i, "H3K27me3"] == 1) {
m[i, "prediction"] <- "n"
} else {
m[i, "prediction"] <- "y"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
if (m[i, "H3K27me3"] == 1) {
m[i, "prediction"] <- "n"
} else {
m[i, "prediction"] <- "y"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
# if (m[i, "H3K27me3"] == 1) {
#
#   m[i, "prediction"] <- "n"
#
# } else {
#
#   m[i, "prediction"] <- "y"
#
# }
m[i, "prediction"] <- "y"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- "y"
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
# if (m[i, "H3K27me3"] == 1) {
#
#   m[i, "prediction"] <- "n"
#
# } else {
#
#   m[i, "prediction"] <- "y"
#
# }
m[i, "prediction"] <- "y"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- "n"
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
# if (m[i, "H3K27me3"] == 1) {
#
#   m[i, "prediction"] <- "n"
#
# } else {
#
#   m[i, "prediction"] <- "y"
#
# }
m[i, "prediction"] <- "y"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- "y"
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- "n"
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
if (m[i, "tss_distance"] > 100000) {
m[i, "prediction"] <- "n"
} else {
m[i, "prediction"] <- "y"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.75) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ( m[i, "sum"] == 0 ) {
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.75) {
m[i, "prediction"] <- "y"
} else {
m[i, "prediction"] <- "n"
}
} else {
m[i, "prediction"] <- "y"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m$prediction <- ifelse(m$sum == 0, "n", "y")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m[m$is_eQTL != m$prediction, ]
table(m$is_eQTL)
m$prediction <- ifelse(m$cv < 0.8 | m$slope > 0.5, "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m$prediction <- ifelse(m$sum == 0, "n", "y")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m[m$is_eQTL != m$prediction, ]
m$prediction <- ifelse(m$sum == 0 & (m$cv < 0.8 | m$slope > 0.5), "y", m$prediction)
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m$prediction <- ifelse(m$sum == 0, "n", "y")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m[m$is_eQTL != m$prediction, ]
m$prediction <- ifelse(m$sum == 0 & (m$cv < 0.8 | m$slope > 0.5), "y", "n")
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else if (m[i, "sum"] == 0) {
m[i, "prediction"] <- "n"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else if (m[i, "sum"] == 0 & cv > 5) {
m[i, "prediction"] <- "n"
}
}
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else if ((m[i, "sum"] == 0) & (m[i, "cv"] > 5)) {
m[i, "prediction"] <- "n"
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if (m[i, "cv"] < 0.8 | m[i, "slope"] > 0.5) {
m[i, "prediction"] <- "y"
} else if ((m[i, "cv"] > 5)) {
if (m[i, "sum"] == 0) {
m[i, "prediction"] <- "n"
}
else {
m[i, "prediction"] <- "y"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
m <- m[m$cv < 0.8 | m$cv > 5, ]
m$prediction <- NA
for ( i in 1:nrow(m) ){
if ((m[i, "cv"] > 5)) {
if (m[i, "sum"] == 0) {
m[i, "prediction"] <- "n"
}
else {
m[i, "prediction"] <- "y"
}
}
}
# m$prediction <- ifelse(m$sum == 0, "n", "y") # + | (m$sum < 2 & m$H3K27me3 == 1))
confusionMatrix(reference = as.factor(m$is_eQTL), data = as.factor(m$prediction), positive = "y")
m <- m.safe
# good
m.pos <- m[m$cv < 0.8 | m$sum > 4, ] #| (m$tss_distance < 10000 & m$sum > 1 & m$cv > 5), ]
m.pos$prediction <- "y"
nrow(m.pos[m.pos$prediction == m.pos$is_eQTL, ]) / nrow(m.pos)
m.neg <- m[m$cv > 5 & ((m$sum == 0 & m$tss_distance < 10000) | (m$sum < 2 & m$H3K27me3 == 1)), ]
m.neg$prediction <- "n"
nrow(m.neg[m.neg$prediction == m.neg$is_eQTL, ]) / nrow(m.neg)
m.neg <-  m[m$cv > 2 | (m$tss_distance < 100000 & m$sum == 0), ]
m.neg$prediction <- "n"
nrow(m.neg[m.neg$prediction == m.neg$is_eQTL, ]) / nrow(m.neg)
m.neg <- m[m$cv > 5 & ((m$sum == 0 & m$tss_distance < 10000)), ] #  | (m$sum < 2 & m$H3K27me3 == 1)
m.neg$prediction <- "n"
nrow(m.neg[m.neg$prediction == m.neg$is_eQTL, ]) / nrow(m.neg)
savehistory("/nfs/users/rg/bborsari/eQTLs.model.nf/bin/simplified.history.R")
