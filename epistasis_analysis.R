### Inputs ###

directory = "C:/Users/Karol Buda/Desktop/Program/R/Epistasis Model/Palmer_et_al"
file = "palmer_et_al.csv"
error = F
p_on = T

### Library ###

library(gtools)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

### Script ###

setwd(directory)

d = read.csv(file)

## Set up data differently if error = T

if(error) {
  d1 = d[5:dim(d)[1],1:3]
  data = data.frame(row.names = 1:(dim(d1)[1]*3))
  data = cbind(data, rep(d1[,1], each = 3))
  replicates = c()
  for(i in 1:dim(d1)[1]){
    replicates = c(replicates, rnorm(3, mean = as.numeric(d1[,2][i]), sd = as.numeric(d1[,3][i])))
  }
  replicates = log10(replicates/as.numeric(d1[which(d1[,1] == d[1,1]), 2]))
  data = cbind(data, replicates)
} else {
  data = d[5:dim(d)[1],1:2]
}

codes = data.frame()
order = data.frame()
simplex_chart = data.frame()

for(j in 1:dim(data)[1]){
  code = c()
  chr = character()
  for(i in 1:length(strsplit(d[1,1], "")[[1]])){
    if(strsplit(data[j,1], "")[[1]][i] == strsplit(d[1,1], "")[[1]][i]){
      code = c(code, -1)
      chr = paste0(chr, strsplit(data[j,1], "")[[1]][i])
    } else {
      code = c(code, 1)
      chr = paste0(chr, strsplit(data[j,1], "")[[1]][i])
    }
  }
  codes = rbind(codes, code)
  simplex_chart = rbind(simplex_chart, c(as.numeric(code), chr))
}

simplex_chart = unique(simplex_chart) ## Filter out repeats

codes = cbind(codes, data[,2])
colnames(codes) = c(paste("p", d[3,][!is.na(d[3,])], sep = ""), "effect")
colnames(simplex_chart) = c(paste("p", d[3,][!is.na(d[3,])], sep = ""), "genotype")

test_codes = gtools::permutations(2, (dim(codes)[2] - 1), c(-1, 1), repeats = TRUE)

for(j in 1:dim(test_codes)[1]) {
  passed = F
  for(i in 1:dim(codes)[1]){
    if(all(test_codes[j,] == codes[-dim(codes)[2]][i,])) {
      passed = T
    }
  }
  if(!passed){
    print(paste0("Data for ", paste(test_codes[j,], collapse = " "), " is missing"))
  }
}

## Loop through each model order and see if its significant

vars = paste(colnames(codes)[1:length(colnames(codes))-1], collapse="+")

## No interaction terms model
current_model = lm(paste("effect ~ (",vars ,")", sep=""), codes)
order = rbind(order, c("order1", summary(current_model)$r.squared, 
                   summary(current_model)$adj.r.squared, 
                   summary(current_model)$r.squared))

for(k in 2:(dim(codes)[2]-1)){
  ## Needed to paste with **k because lm doesn't like the power term introduced otherwise
  next_model = lm(paste("effect ~ (",vars ,")**", k, sep=""), codes)
  
  if(p_on) {
    if(!is.na(anova(current_model, next_model)[6][2,]) & anova(current_model, next_model)[6][2,] < 0.05) {
      current_model = next_model
      # Add order into order list
      order = rbind(order, c(paste("order", k, sep=""), summary(current_model)$r.squared, 
                         summary(current_model)$adj.r.squared, 
                         summary(current_model)$r.squared - as.numeric(order[k-1,2])))
    } else{
      print(paste0("Highest Model Order: ", k-1))
      print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
      break
    }
    print(paste0("Analyzing Next Model Order: ", k))
  } else {
    if(!is.na(anova(current_model, next_model)[6][2,])) {
      current_model = next_model
      # Add order into order list
      order = rbind(order, c(paste("order", k, sep=""), summary(current_model)$r.squared, 
                             summary(current_model)$adj.r.squared, 
                             summary(current_model)$r.squared - as.numeric(order[k-1,2])))
    } else{
      print(paste0("Highest Model Order: ", k-1))
      print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
      break
    }
    print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
    print(paste0("Setting next model as current model...: "))
    print(paste0("Analyzing Next Model Order: ", k))
  }
}

## Prepare output

pos_out = as.data.frame(summary(current_model)$coef[,1]*2)
pos_out[1,1] = summary(current_model)$coef[,1][1] ## Don't multiply intercept by 2

pos_out$names = rownames(pos_out)
rownames(pos_out) = c()
colnames(pos_out) = c("effect", "indices")

colnames(order) = c("model_order", "R2", "(Adjusted)", "Delta R2")

## Convert : and p characters
pos_out$indices = gsub(':', '|', pos_out$indices)
pos_out$indices = gsub('p', '', pos_out$indices)
pos_out$indices[1] = "INTERCEPT"

## Genotype predicted out

vars_sep = unlist(strsplit(vars, "+", fixed=TRUE))
simplex_chart = cbind(apply(simplex_chart[-dim(simplex_chart)[2]], 2, as.numeric), simplex_chart[dim(simplex_chart)[2]])

preds = c()
for(i in 1:dim(simplex_chart)[1]) {
  preds = c(preds, predict(current_model, simplex_chart[-dim(simplex_chart)[2]][i,]))
}

pred_df = cbind(simplex_chart[dim(simplex_chart)[2]], preds)
colnames(pred_df) = c("genotype", "predicted effect")

## Truncated First order model predicted out vs real data

simplex_chart = simplex_chart %>%
  arrange_at(names(simplex_chart)[-dim(simplex_chart)[2]])

preds = c()
for(i in 1:dim(simplex_chart)[1]) {
  preds = c(preds, predict(lm(paste("effect ~ (",vars ,")", sep=""), codes), simplex_chart[-dim(simplex_chart)[2]][i,]))
}

codes$effect = as.numeric(codes$effect)

effects_vector = codes %>%
  group_by_at(names(codes)[-dim(codes)[2]]) %>%
  summarise(avg = mean(effect)) %>%
  pull(avg)

pred_compare_df = cbind(simplex_chart[dim(simplex_chart)[2]], preds, effects_vector)

colnames(pred_compare_df) = c("genotype", "truncated effect", "observed effect")

## Make plot nicer and auto_adjust scale to make more sense

ggplot(pred_compare_df, aes(x = `truncated effect`, y = `observed effect`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = 'dashed') +
  geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = 'dashed') +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0 + log10(1.5), linetype = 'dashed') +
  geom_vline(xintercept = 0 - log10(1.5), linetype = 'dashed') +
  geom_polygon(data=data.frame(x=c(log10(1.5), 
                                   max(pred_compare_df$`truncated effect`) + log10(1.5), 
                                   max(pred_compare_df$`truncated effect`) + log10(1.5)), 
                               y = c(0, 0, 1*(max(pred_compare_df$`truncated effect`)+log10(1.5)) - log10(1.5))), 
               mapping=aes(x=x, y=y), fill = "blue", alpha = 0.2) +
  geom_polygon(data=data.frame(x=c(log10(1.5), 
                                   log10(1.5), 
                                   max(pred_compare_df$`truncated effect`),
                                   max(pred_compare_df$`truncated effect`)), 
                               y = c(log10(1.5)*1 + log10(1.5), 
                                     Inf, 
                                     Inf, 
                                     1*max(pred_compare_df$`truncated effect`) + log10(1.5))), 
               mapping=aes(x=x, y=y), fill = "red", alpha = 0.2) +
  geom_polygon(data=data.frame(x=c(log10(1.5), 
                                   log10(1.5), 
                                   max(pred_compare_df$`truncated effect`),
                                   max(pred_compare_df$`truncated effect`)), 
                               y = c(0, 
                                     min(pred_compare_df$`observed effect`), 
                                     min(pred_compare_df$`observed effect`), 
                                     0)), 
               mapping=aes(x=x, y=y), fill = "darkblue", alpha = 0.3) +
  xlab("Predicted effect using 1st order model") +
  ylab("Observed effect") +
  xlim(NA, max(pred_compare_df$`truncated effect`) + log10(1.5))


plot(`observed effect` ~ `truncated effect`, pred_compare_df)
abline(0, 1)
abline(0 + log10(1.5), 1, lty = 3)
abline(0 - log10(1.5), 1, lty = 3)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0 + log10(1.5), lty = 3, lwd = 0.1)
abline(v = 0 - log10(1.5), lty = 3, lwd = 0.1)

## Writing csvs

write.csv(pos_out[,c(2,1)],"pos_out.csv", row.names = FALSE)
write.csv(order,"model_order.csv", row.names = FALSE)
write.csv(pred_df, "gen_out.csv", row.names = FALSE)
write.csv(pred_compare_df, "pred_out.csv", row.names = FALSE)

## Provide feedback in log.txt form if there are singularities because of missing data AND the p-value of next model if exists

if(length(which(is.na(current_model$coef))) > 0) {
  write(c("-------------", paste("Variable", names(current_model$coef[which(is.na(current_model$coef))]), "= NA")), "log.txt")
}
if(!is.na(anova(current_model, next_model)[6][2,])) {
  write(c("-------------", paste0("Highest Model Order: ", k-1), paste0("Next model p-value: ", anova(current_model, next_model)[6][2,])), "log.txt", append = T)
}


## Plotting

pos_out$order = c(NA, paste0("Order ", str_count(pos_out$indices[-1], "[|]") + 1, " (", formatC(as.numeric(order[str_count(pos_out$indices[-1], "[|]") + 1, 2]), digits = 3), ")"))

pos_out$indices = factor(pos_out$indices, levels = pos_out$indices)

png("epistatic_effect.png", width = 1200, height = 752)

ggplot(pos_out[-1,], aes(x = indices, y = effect, fill = order)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=formatC(effect, digits=2), y = effect + 0.1*sign(effect)), size = 2.5) +
  facet_wrap(~ order, scales = "free_x") + 
  labs(x = "Mutation(s)", y = "Fold Effect on Activity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()

## Qualitative Violin plot for examining clustering/distributions

## Mean Violin Plot

newer_codes = codes %>%
  group_by_at(names(codes)[-length(names(codes))]) %>%
  summarise(avg = mean(effect))

replace_occurence = function(x, replaced, value){
  x[x == replaced] = value
  return(x)
}

codes_effects = c()
codes_identity = c()
codes_position = c()

subtract = function(x) {
  if(length(x) == 1) {
    return(NA)
  } else {
    return(x[2] - x[1])
  }
} ## Very primitive function with little checks, however after taking the mean there should only be 1 or 2 values
## If there is 1 value we can't display it because there is no interaction to compare so we return an NA and remove it

for(i in 1:(length(names(newer_codes)) - 1)) {
  
  current = newer_codes %>%
    group_by_at(names(newer_codes)[(1:(length(codes)-1))[(1:(length(codes)-1)) != i]]) %>%
    summarise(idk = subtract(avg)) %>%
    na.omit()
  
  codes_effects = c(codes_effects, current %>% pull(idk))
  
  for(j in 1:(dim(current)[1])) {
    str_code = paste(as.character(replace_occurence(current[j, 1:(length(names(current)) - 1)], -1, 0)), collapse = "")
    codes_identity = c(codes_identity, paste(c(substr(str_code, 0, i-1), substr(str_code, i, (length(names(newer_codes))-2))), collapse="x"))
  }
  
  codes_position = c(codes_position, rep(names(newer_codes)[i], dim(current)[1]))

}

###

new_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 

new_codes$positions = factor(new_codes$positions, levels = unique(new_codes$positions))

png("mean_fold_effect.png", width = 1200, height = 752)

ggplot(new_codes, aes(x = positions, y = effects, color = positions)) +
  geom_violin(color = "darkgray", trim = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_jitter(position = position_jitter(0)) +
  geom_text_repel(aes(label = identity), force = 0.5, nudge_x = 0.3, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
  labs(x = "Mutant Positions", y = "Mean Fold Effect Change on Activity") +
  theme_classic()

dev.off()

####################
## IN DEVELOPMENT ##
####################

# Looking at ALL the average effects, not just first order

## Remove all memory for next analysis
rm(list=ls())