
## SET WORKING DIRECTORY TO THE LOCATION OF **THIS** FILE

# read the data
tmp = read.csv("../bio/adult-indiv-carcass.csv", stringsAsFactors = F)

# drop out non-known females
tmp = tmp[tmp$sex == "F",]

# drop out non-known spawn status
tmp = tmp[tmp$prespawn != "Unk" & !is.na(tmp$prespawn),]

# coerce year to factor
tmp$year = as.factor(tmp$year)

# add a count variable: for use in summing records
tmp$count = 1

# aggregate data: total carcasses sampled by population and year
sampled = aggregate(count ~ population + year, tmp, sum)
colnames(sampled)[3] = "sampled"

# aggregate data: total carcasses sampled that were spawned out by population and year
success = aggregate(count ~ population + year, subset(tmp, prespawn == "Spawned"), sum)
colnames(success)[3] = "yes"

# merge these two data sets
tmp = merge(sampled, success, all = TRUE)
tmp$yes[is.na(tmp$yes)] = 0
tmp$no = tmp$sampled - tmp$yes

# fit the population-specific mixed models
# assume time-constant mean with logit-normal annual random effects
fit_list = list(
  CAT = glmmTMB::glmmTMB(cbind(yes, no) ~ 1 + (1|year), family = binomial, data = subset(tmp, population == "CAT")),
  LOS = glmmTMB::glmmTMB(cbind(yes, no) ~ 1 + (1|year), family = binomial, data = subset(tmp, population == "LOS")),
  MIN = glmmTMB::glmmTMB(cbind(yes, no) ~ 1 + (1|year), family = binomial, data = subset(tmp, population == "MIN")),
  UGR = glmmTMB::glmmTMB(cbind(yes, no) ~ 1 + (1|year), family = binomial, data = subset(tmp, population == "UGR"))
)

# obtain model-predicted values
preds = lapply(fit_list, function(fit) {
  pred_data = data.frame(year = unique(fit$frame$year))
  pop = stringr::str_extract(as.character(fit$call)[3], "([A-Z][A-Z][A-Z])")
  preds = predict(fit, newdata = pred_data, type = "response")
  out = data.frame(year = pred_data$year, preds)
  colnames(out)[2] = pop
  out
})

# format model predicted values
preds = merge(preds[[1]], preds[[2]], by = "year", all = TRUE) |>
  merge(preds[[3]], by = "year", all = TRUE) |>
  merge(preds[[4]], by = "year", all = TRUE)
preds$year = as.numeric(as.character(preds$year))

# replace any missing values with the mean of all non-missing values
preds = as.data.frame(apply(preds, 2, function(x) {x[is.na(x)] = mean(x, na.rm = TRUE); x}))

# cap any values less than 0.4: only happens for several UGR years
preds[preds < 0.4] = 0.4

# round
preds[,2:5] = round(preds[,2:5], 3)

# reformat
preds = reshape2::melt(preds, id.vars = "year", value.name = "prespawn_surv", variable.name = "population")

# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(preds)[1] = "brood_year"

# export the data file
write.csv(preds, "../bio/prespawn-surv.csv", row.names = FALSE)
