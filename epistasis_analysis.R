### Inputs ###

directory = ""
file = ""
error = F

### Library ###

library(gtools)
library(ggplot2)
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

## Writing csvs

write.csv(pos_out[,c(2,1)],"pos_out.csv", row.names = FALSE)
write.csv(order,"model_order.csv", row.names = FALSE)
write.csv(pred_df, "gen_out.csv", row.names = FALSE)

## Provide feedback in log.txt form if there are singularities because of missing data AND the p-value of next model if exists

if(length(which(is.na(current_model$coef))) > 0) {
  write(c("-------------", paste("Variable", names(current_model$coef[which(is.na(current_model$coef))]), "= NA")), "log.txt")
}
if(!is.na(anova(current_model, next_model)[6][2,])) {
  write(c("-------------", paste0("Highest Model Order: ", k-1), paste0("Next model p-value: ", anova(current_model, next_model)[6][2,])), "log.txt", append = T)
}


## Plotting

pos_out$order = c(NA, paste0("Order ", str_count(pos_out$indices[-1], "[|]") + 1))

pos_out$indices = factor(pos_out$indices, levels = pos_out$indices)

png("fold_effect.png", width = 1200, height = 752)

ggplot(pos_out[-1,], aes(x = indices, y = effect, fill = order)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=formatC(effect, digits=2), y = effect + 0.1*sign(effect)), size = 2.5) +
  facet_wrap(~ order, scales = "free_x") + 
  labs(x = "Mutation(s)", y = "Effect on -deldelG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()

## Qualitative Violin plots for examining clustering/distributions

codes_effects = c()
codes_position = c()

for(i in 1:(length(names(codes)) - 1)) {
  codes_effects = c(codes_effects, as.numeric(subset(codes, codes[,i] == 1)$effect))
  codes_position = c(codes_position, rep(names(codes)[i], length(as.numeric(subset(codes, codes[,i] == 1)$effect))))
}

new_codes = data.frame(positions = codes_position, effects = codes_effects) 

png("fold_effect.png", width = 1200, height = 752)

ggplot(new_codes, aes(x = positions, y = effects, color = positions)) +
  geom_violin(color = "darkgray", trim = FALSE) +
  geom_jitter(position = position_jitter(0.2)) +
  labs(x = "Mutant Positions", y = "Fold Effect Change on Activity") +
  theme_classic()

dev.off()

##

codes$effect = as.numeric(codes$effect)

newer_codes = codes %>%
  group_by_at(names(codes)[-length(names(codes))]) %>%
  summarise(avg = mean(effect))

replace_occurence = function(x, replaced, value){
  x[x == replaced] = value
  return(x)
}

codes_iden = c()

for(i in 1:dim(newer_codes)[1]) {
  codes_iden = c(codes_iden, paste(as.character(replace_occurence(newer_codes[i, 1:(length(names(newer_codes)) - 1)], -1, 0)), collapse = ""))
}

newer_codes$identity = codes_iden

codes_effects = c()
codes_identity = c()
codes_position = c()

for(i in 1:(length(names(newer_codes)) - 1)) {
  codes_effects = c(codes_effects, as.numeric(subset(newer_codes, newer_codes[,i] == 1)$avg))
  codes_identity = c(codes_identity, as.character(subset(newer_codes, newer_codes[,i] == 1)$identity))
  codes_position = c(codes_position, rep(names(newer_codes)[i], length(as.numeric(subset(newer_codes, newer_codes[,i] == 1)$avg))))
}

newer_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 

png("mean_fold_effect.png", width = 1200, height = 752)

ggplot(newer_codes, aes(x = positions, y = effects, color = positions)) +
  geom_violin(color = "darkgray", trim = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_jitter(position = position_jitter(0)) +
  geom_text_repel(aes(label = identity), force = 0.5, nudge_x = 0.5, direction = "y", box.padding = 0, segment.size = 0.2, size = 3) +
  labs(x = "Mutant Positions", y = "Mean Fold Effect Change on Activity") +
  theme_classic()

dev.off()

