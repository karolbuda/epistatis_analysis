### Inputs ###

directory = "C:/Users/karol/OneDrive - The University Of British Columbia (1)/PhD/Epistasis_Lit/AtzA/atrazine"
file = "noor_et_al_atz.csv"
error = F
p_on = F
sign = T
log_base = 10
flip = F

### Library ###

library(gtools)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(igraph) ## Network library
library(gtools) ## Permutation library
library(e1071) ## Hamming distance lib

### Functions ###

data_loading = function(error) {
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
    replicates = log(replicates/as.numeric(d1[which(d1[,1] == d[1,1]), 2]), log_base)
    data = cbind(data, replicates)
  } else {
    data = d[5:dim(d)[1],1:2]
  }
  
  codes = data.frame()
  #simplex_chart = data.frame()
  assign("simplex_chart", data.frame(), envir = .GlobalEnv)
  
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
  codes$effect = as.numeric(codes$effect)
  codes <<- codes
  simplex_chart <<- simplex_chart
}

epistasis_analysis = function(p_on, sign) {
  
  order = data.frame()
  
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
        print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
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
  
  #preds = c()
  #for(i in 1:dim(simplex_chart)[1]) {
  #  preds = c(preds, predict(lm(paste("effect ~ (",vars ,")", sep=""), codes), simplex_chart[-dim(simplex_chart)[2]][i,]))
  #}
  
  #effects_vector = codes %>%
  #  group_by_at(names(codes)[-dim(codes)[2]]) %>%
  #  summarise(avg = mean(effect)) %>%
  #  pull(avg)
  
  #pred_compare_df = cbind(simplex_chart[dim(simplex_chart)[2]], preds, effects_vector)
  
  #colnames(pred_compare_df) = c("genotype", "truncated effect", "observed effect")
  
  ## Writing csvs
  
  write.csv(pos_out[,c(2,1)],"pos_out.csv", row.names = FALSE)
  write.csv(order,"model_order.csv", row.names = FALSE)
  write.csv(pred_df, "gen_out.csv", row.names = FALSE)
  #write.csv(pred_compare_df, "pred_out.csv", row.names = FALSE)
  
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
  
  ggplot(pos_out[-1,], aes(x = indices, y = effect, fill = order)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=formatC(effect, digits=2), y = effect + 0.1*sign(effect)*abs(max(effect))), size = 3) +
    facet_wrap(~ order, scales = "free_x") + 
    labs(x = "Mutation(s)", y = "Fold Effect on Activity") +
    theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  ggsave("epistatic_effect.png", units = "in", width = 12, height = 9)
  
  ###########################
  ## MEAN FOLD EFFECT PLOT ##
  ###########################
  
  ## Combining points into averages 
  
  newer_codes = codes %>%
    group_by_at(names(codes)[-length(names(codes))]) %>%
    summarise(avg = log(signif(mean(log_base^effect), digits = 3), log_base))
  # Code above ensures that I don't take mean of the logs, first I transform them back with 10^x
  # I also use the signif function because too many digits makes log10(1) != 0 for some reason
  
  replace_occurence = function(x, replaced, value){
    x[x == replaced] = value
    return(x)
  }
  # Used to remove values in data set which don't have either a WT-like or mutant-like genotype... i.e. we can't show
  # x0001 if we don't have both 00001 and 10001 in our data, so we remove the incorrectly calculated x0001
  
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
    
    if(dim(current)[1] > 0) {
      codes_effects = c(codes_effects, current %>% pull(idk))
      
      for(j in 1:(dim(current)[1])) {
        str_code = paste(as.character(replace_occurence(current[j, 1:(length(names(current)) - 1)], -1, 0)), collapse = "")
        codes_identity = c(codes_identity, paste(c(substr(str_code, 0, i-1), substr(str_code, i, (length(names(newer_codes))-2))), collapse="x"))
      }
      
      codes_position = c(codes_position, rep(names(newer_codes)[i], dim(current)[1]))
    }
    
  }
  
  new_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 
  
  new_codes$mut = factor(str_count(new_codes$identity, "1") + 1)
  
  new_codes$positions = factor(new_codes$positions, levels = unique(new_codes$positions))
  
  new_codes %>%
    ggplot(aes(x = positions, y = effects, color = positions)) +
      geom_boxplot(color = "darkgray") +
      geom_hline(yintercept=0, linetype="dashed", color = "black") +
      geom_point() +
      stat_summary(fun=mean, colour="darkred", geom="crossbar", width=0.2) + 
      geom_text_repel(aes(label = identity), force = 0.5, nudge_x = 0.3, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
      labs(x = "Mutant Positions", y = paste(c("Delta log", log_base, "(effect) of mutation"), collapse = "")) +
      theme(text = element_text(size=14))
  
  ggsave("mean_fold_effect.png", units = "in", width = 12, height = 9)

  if(sign == T) {
    ## Genotype predicted out
    
    effects_vector = newer_codes %>% pull(avg)
    
    vars_sep = unlist(strsplit(vars, "+", fixed=TRUE))
    simplex_chart = cbind(apply(simplex_chart[-dim(simplex_chart)[2]], 2, as.numeric), simplex_chart[dim(simplex_chart)[2]])
    
    ## Arrange simplex chart the same way for all
    
    simplex_chart = simplex_chart %>%
      arrange_at(names(simplex_chart)[-dim(simplex_chart)[2]])
    
    preds = c()
    for(i in 1:dim(simplex_chart)[1]) {
      preds = c(preds, predict(current_model, simplex_chart[-dim(simplex_chart)[2]][i,]))
    }
    
    pred_df = cbind(simplex_chart, effects_vector, preds)
    
    ## Truncated First order model predicted out vs real data
    
    preds = c()
    for(i in 1:dim(simplex_chart)[1]) {
      preds = c(preds, predict(lm(paste("effect ~ (",vars ,")", sep=""), codes), simplex_chart[-dim(simplex_chart)[2]][i,]))
    }
    
    pred_df = cbind(pred_df, preds)
    
    ### WT vs average effect
    
    wt_effects = new_codes$effect[which(new_codes$identity %in% new_codes$identity[-grep("1", new_codes$identity)])]
    
    ## Adds NA to the places where variable is missing
    
    missing = which(!(vars_sep %in% as.character(new_codes[which(new_codes$identity %in% new_codes$identity[-grep("1", new_codes$identity)]),]$positions)))
    
    if(length(missing) > 0) {
      for(i in 1:length(missing)) {
        wt_effects = append(wt_effects, 0, after=(missing[i]-1))
      }
    }
  
    wt_vs_avg = data.frame(indices = pos_out[2:(1+length(wt_effects)),]$indices, average_effect = pos_out[2:(1+length(wt_effects)),]$effect, wt_effect = wt_effects)
    
    wt_vs_avg_plot = data.frame(index = rep(wt_vs_avg$indices, 2), effects = c(wt_vs_avg$average_effect, wt_vs_avg$wt_effect), 
                                type = c(rep("average", dim(wt_vs_avg)[1]), rep("wt", dim(wt_vs_avg)[1])))
    
    wt_codes = codes[-length(codes)]
    wt_codes[wt_codes == -1] = 0
    wt_codes = unique(wt_codes)
    
    wt_codes = wt_codes %>%
      arrange_at(names(wt_codes))
    
    
    wt_model = lm(paste("effect ~ (",vars ,")", sep=""), codes)
    
    #int = codes %>%
      #filter_at(vars(-effect), all_vars(. == -1)) %>%
      #summarise(mean_effect = mean(effect)) %>%
      #pull(mean_effect)
    
    #wt_model$coefficients = c(as.numeric(int), wt_effects)
    wt_model$coefficients = c(0, wt_effects)
    names(wt_model$coefficients) = names(lm(paste("effect ~ (",vars ,")", sep=""), codes)$coefficients)
    
    preds = c()
    for(i in 1:dim(wt_codes)[1]) {
      preds = c(preds, predict(wt_model, wt_codes[i,]))
    }
    
    #for(i in 1:length(as.numeric(which(is.na(wt_model$coefficients))))) {
    #  preds[which(wt_codes[as.numeric(which(is.na(wt_model$coefficients)))[i] - 1] == 1)] = NA
    #}
    
    pred_df = cbind(pred_df, preds)
    
    colnames(pred_df) = c(colnames(simplex_chart), "observed effect", "epistatic effect", "average effect", "WT effect")
    
    if(length(missing) > 0) {
      for(i in 1:length(missing)) {
        pred_df$`WT effect`[which(pred_df[,missing[i]] == 1)] = NA
      }
    }

    mutations = c()
    for(i in 1:dim(pred_df)[1]){
      mutations = c(mutations, length(which(pred_df[i,1:dim(codes)[2]-1] == 1)))
    }
    
    pred_df$mutations = as.factor(mutations)
    
    #ggplot(pred_df, aes(`epistatic effect`, `observed effect`, color = mutations)) +
    #  geom_abline(slope = 1, intercept = 0) +
    #  geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = "dashed") +
    #  geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = "dashed") +
    #  geom_vline(xintercept = 0) +
    #  geom_hline(yintercept = 0) +
    #  geom_point() +
    #  geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
    #  theme(text = element_text(size=14))
    
    ### I stopped printing the models and outputting the average and epistatic figures as they're not too informative
    
    #print(summary(lm(`observed effect` ~ `epistatic effect`, data = pred_df)))
    
    #ggsave("observed_vs_epistatic.png", units = "in", width = 12, height = 9)
    
    #ggplot(pred_df, aes(`average effect`, `observed effect`, color = mutations)) +
    #  geom_abline(slope = 1, intercept = 0) +
    #  geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = "dashed") +
    #  geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = "dashed") +
    #  geom_vline(xintercept = 0) +
    #  geom_hline(yintercept = 0) +
    #  geom_point() +
    #  geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
    #  theme(text = element_text(size=14))
    
    #print(summary(lm(`observed effect` ~ `average effect`, data = pred_df)))
    
    #ggsave("observed_vs_average.png", units = "in", width = 12, height = 9)
    
    ggplot(pred_df, aes(`WT effect`, `observed effect`, color = mutations)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 1, intercept = 0 + log(1.5, log_base), linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0 - log(1.5, log_base), linetype = "dashed") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_point() +
      geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
      theme(text = element_text(size=14)) +
      xlab(paste(c("Predicted log", log_base, "(effect)"), collapse = "")) +
      ylab(paste(c("Observed log", log_base, "(effect)"), collapse = "")) +
      ggtitle("Observed vs Predicted Effect based on Single Mutants")
    
    #print(summary(lm(`observed effect` ~ `WT effect`, data = pred_df)))
    
    ggsave("observed_vs_wt.png", units = "in", width = 12, height = 9)
    
    ##################################################
    ## Quick ratio calculation for pymol csv output ##
    ###################################################
    
    pred_df <<- pred_df
    new_codes <<- new_codes
  }  
}

blanket_analysis = function() {
  all_stats = new_codes %>%
    group_by(positions) %>%
    summarise(mean = mean(effects),
              median = median(effects),
              top = max(effects), 
              bottom = min(effects),
              range = max(effects) - min(effects),
              skew = (3*(mean(effects) - median(effects)))/sd(effects),
              pos_pct = sum(effects > 0) / length(effects),
              neg_pct = sum(effects < 0) / length(effects))
  
  types = pred_df %>%
    filter((as.numeric(mutations) - 1) > 1) %>%
    group_by(genotype) %>%
    summarise(magnitude = log_base^(abs(`observed effect` - `WT effect`)),
              sign = `observed effect`*`WT effect`,
              negative = `observed effect` < `WT effect`) %>%
    mutate(magnitude = magnitude > 1.5,
           sign = sign < 0) %>%
    filter(magnitude == T | sign == T)
  
  # negative filters positive or negative, TRUE means negative epistasis
  # sign trumps magnitude so only count magnitude when sign is false and only count sign when magnitude is true
  # This ensures over counting and rest should always be positive or 0
  
  neg_mag = sum(types[types$negative == T & types$sign == F,]$magnitude)
  pos_mag = sum(types[types$negative == F & types$sign == F,]$magnitude)
  neg_sg = sum(types[types$negative == T & types$magnitude == T,]$sign)
  pos_sg = sum(types[types$negative == F & types$magnitude == T,]$sign) 
  rest = (dim(pred_df[(as.numeric(pred_df$mutations) - 1) > 1,])[1]) - neg_mag - pos_mag - pos_sg - neg_sg
  
  labs = c(pos_mag, neg_mag, pos_sg, neg_sg, rest)
  labs = replace(labs, labs == 0, NA)
  
  png(filename = "pie_chart.png", width = 790, height = 840, units = "px")
  
  pie(c(pos_mag, neg_mag, pos_sg, neg_sg, rest), labels = labs, col = 
        c("#ffe4a6", "#b8eaf5", "#a997cf", "#fcc4ae","#d9d9d9"))
  
  legend("topright", c("Positive Magnitude", "Negative Magnitude", "Positive Sign", "Negative Sign", "None"), cex = 0.8,
         fill = c("#ffe4a6", "#b8eaf5", "#a997cf", "#fcc4ae","#d9d9d9"))
  
  dev.off()
  
  write.csv(all_stats, "all_stats.csv", row.names = F)
  
  ## Pymol ratio export
  
  pymol_csv = simplex_chart[-length(simplex_chart)]
  
  ind = c()
  for(i in 1:dim(pymol_csv)[1]) {
    if(length(colnames(pymol_csv[i,][which(pymol_csv[i, ] == 1)])) > 0) {
      ind = c(ind, paste(colnames(pymol_csv[i,][which(pymol_csv[i, ] == 1)]), collapse = "|"))
    } else {
      ind = c("INTERCEPT")
    }
  }
  
  ind = gsub('p', '', ind)
  
  #pymol_csv = data.frame(indices = ind, effect = 10^(pred_df$`observed effect`)/10^(pred_df$`WT effect`)) %>% 
  #  filter(effect > 1.5 | effect < (1/1.5)) %>%
  #  mutate(effect = log10(effect))
  
  pymol_csv = data.frame(indices = ind, 
                         magnitude = log_base^(abs(pred_df$`observed effect` - pred_df$`WT effect`)),
                         sign = pred_df$`observed effect`*pred_df$`WT effect`,
                         negative = pred_df$`observed effect` < pred_df$`WT effect`) %>%
    mutate(effect = log(magnitude, log_base)*c(1,-1)[as.numeric(negative) + 1],
           magnitude = magnitude > 1.5,
           sign = as.numeric(sign < 0)) %>%
    filter(magnitude == T | sign == T) %>%
    select(-c(magnitude, negative))
  
  write.csv(pymol_csv, "pymol_csv.csv", row.names = F)
  
}

network_analysis = function(eff = "observed") {
  pred_df = pred_df[order(pred_df$mutations),]
  
  mut = length(strsplit(pred_df$genotype[1], "")[[1]]) ## Get length of WT as indicator of mutations
  
  ids = pred_df[,1:mut]
  ids[ids == -1] = 0
  ids = apply(ids, 1, function(x) paste(as.character(x), collapse = ""))
  
  
  if(eff == "WT") {
    codes = data.frame(id = ids, effect = pred_df$`WT effect`)
  } else if(eff == "average") {
    codes = data.frame(id = ids, effect = pred_df$`average effect`)
  } else if(eff == "epistatic") {
    codes = data.frame(id = ids, effect = pred_df$`epistatic effect`)
  } else if(eff == "observed") {
    codes = data.frame(id = ids, effect = pred_df$`observed effect`)
  } else if(eff == "ratio") {
    codes = data.frame(id = ids, effect = log_base^(pred_df$`observed effect`)/log_base^(pred_df$`WT effect`))
  }
  
  muts = c()
  for(i in 1:length(codes$id)) {
    muts = c(muts, sum(as.numeric(str_split(codes$id[i], "")[[1]])))
  }
  
  codes = codes[order(muts),] # order by mutation level
  codes$muts = muts[order(muts)]
  
  codes_perm = data.frame()
  
  for(i in 1:length(codes$id)) {
    codes_perm = rbind(codes_perm, as.numeric(str_split(codes$id[i], "")[[1]]))
  }
  
  codes_perm = as.matrix(codes_perm)
  
  from = c()
  to = c()
  
  for(i in 1:dim(codes_perm)[1]) {
    for(j in 1:dim(codes_perm)[1]) {
      if(hamming.distance(codes_perm[i, ], codes_perm[j, ]) == 1 & sum(codes_perm[i, ]) < sum(codes_perm[j, ])) {
        ## Hamming distance ensures only 1 change and sum ensures no back tracking
        from = c(from, codes$id[i])
        to = c(to, codes$id[j])
      }
    }
  }
  
  links = data.frame(from = from, to = to)
  
  ## change effects with data.frame that contains IDs and effect
  
  pascalTriangle <- function(h) {
    lapply(0:h, function(i) choose(i, 0:i))
  }
  
  # Define max vector
  
  weights = c()
  magnitude = c()
  
  pascal_terms = tail(pascalTriangle(mut),1)[[1]][-length(tail(pascalTriangle(mut),1)[[1]])]
  
  for(i in 1:dim(links)[1]) {
    
    dif = codes$effect[which(codes$id %in% links[i,2])] - codes$effect[which(codes$id %in% links[i,1])]
    magnitude = c(magnitude, dif)
    
    if(is.na(dif)) {
      weights = c(weights, 1)
    } else if(dif >= 0 & eff != "ratio") {
      weights = c(weights, 2)
    } else {
      weights = c(weights, 1)
    }
  }
  
  links = cbind(links, weights, magnitude)
  
  ## Currently spacing is based on pascal terms instead of existing terms... may need to change that for missing data
  
  ## x coordinates of nodes
  
  test = c(0)
  
  for(i in unique(codes$muts)[-1]) {
    test = c(test, rep(i, length(codes$muts[codes$muts == i])))
  }
  
  
  ## y coordinates of nodes
  
  pasc_terms = c()
  for(i in unique(codes$muts)) {
    pasc_terms = c(pasc_terms, length(codes$muts[codes$muts == i]))
  }
  
  pasc_terms = pasc_terms[-1]
  
  test_2 = c(max(pasc_terms)/2)
  
  for(x in 1:length(pasc_terms)) {
    test_2 = c(test_2, max(pasc_terms)-(max(pasc_terms)/(pasc_terms[x]+1))*(1:(pasc_terms[x])))
  }
  
  l <- matrix(c(test, test_2), nrow = (tail(cumsum(pasc_terms), 1) + 1), ncol = 2)
  
  if(flip) {
    links$weights[links$weights == 1] = 0
    links$weights[links$weights == 2] = 1
    links$weights[links$weights == 0] = 2
  }
  
  net_test = graph.data.frame(links, codes, directed = T)
  
  ## Calculate shortest paths
  
  #paths_2 = shortest_paths(net_test, from = codes$id[1], to = codes$id[length(codes$id)], 
  #                         weights = 1/(E(net_test)$magnitude + max(c(abs(min(E(net_test)$magnitude)), abs(max(E(net_test)$magnitude)))) + 1))
  
  if(eff != "ratio") {
  
    net_test_pos = graph.data.frame(links[links$weights > 1,], codes, directed = T)
    
    paths_pos = all_simple_paths(net_test_pos, from = codes$id[1], to = codes$id[which.max(codes$effect)])
    
    # Vector in order of edges that don't lead to dead ends
    
    if(length(paths_pos) > 0) {
      pos = rep(1, dim(links)[1])
      for(i in 1:length(paths_pos)) {
        paths_pos_cur = as_ids(paths_pos[[i]])
        for(j in 1:dim(links)[1]){
          for(i in 1:(length(paths_pos_cur) - 1)) {
            if(links[j,1] == paths_pos_cur[i] & links[j,2] == paths_pos_cur[i+1]) {
              pos = replace(pos, j, 2)
            }
          }
        }
      }
    }

  }
  ## Creates gradient of colors
  
  rbPal = colorRampPalette(c('#45d7ff','white', '#ff4e45'))
  
  ## This sets 0 in the center with max and min being equal to pos and neg of the largest value
  idk = max(c(abs(max(codes$effect, na.rm = T)), abs(min(codes$effect, na.rm = T))))
  col = rbPal(length(seq(-idk, idk, idk/50)))[as.numeric(cut(seq(-idk, idk, idk/50), breaks = length(seq(-idk, idk, idk/50))))]
  
  cols = c()
  #val = c()
  for(i in 1:length(codes$effect)) {
    if(!is.na(codes$effect[i])){
      cols = c(cols, col[which.min(abs(seq(-idk, idk, idk/50) - codes$effect[i]))])
      #val = c(val, codes$effect[which.min(abs(seq(-idk, idk, idk/50) - codes$effect[i]))])
    } else {
      cols = c(cols, 'black')
    }
  }
  
  codes$colors = cols
  
  ## Colors for ratio plot
  if(eff == "ratio") {
    
    ratio_effects = codes$effect
    
    no_imp = which(is.na(pred_df$`WT effect`))[which(which(is.na(pred_df$`WT effect`)) %in% which(is.na(pred_df$`observed effect`)))]
    
    ratio_effects[no_imp] = 1e-4 # Place holder value to definitely make these negative and positive in error checking
    ratio_effects[is.na(ratio_effects)] = 1e4
    
    cols = c()
    
    for(i in 1:length(ratio_effects)) {
      if(ratio_effects[i] > 1) {
        cols = c(cols, '#ff4e45')
      } else if(ratio_effects[i] < 1) {
        cols = c(cols, '#45d7ff')
      } else {
        cols = c(cols, 'white')
      }
    }
  }
  
  V(net_test)$frame.color <- "black"
  if(eff != "ratio") {
    V(net_test)$color <- codes$colors
  } else {
    V(net_test)$color <- cols
  }
  V(net_test)$size <- 15
  E(net_test)$lty <- c(2, 1)[links$weights]
  E(net_test)$width <- E(net_test)$weights
  E(net_test)$arrow.mode <- 0
  if(eff != "ratio") {
    if(length(paths_pos) > 0) {
      E(net_test)$color <- c("grey", "black")[pos]
    }
  } else {
    E(net_test)$color <- "grey"
  }
  
  #for (p in paths_2$vpath) { E(net_test, path=p)$color <- "red" }
  
  ## Finds the optimal path along each node
  
  ## Append to your optimal path algorithm
  ## If you encounter no good options, move back? This is a hard
  
  #logger = c()
  #counter = links$from[1]
  
  #while(T) {
  #  if(length(which(links$from == counter)) < 1) {
  #    break
  #  }
  #  if(max(links$magnitude[links$from == counter]) < 0) {
  #    break
  #  }
  #  tester = which(links$magnitude == max(links$magnitude[links$from == counter]) & links$from == counter)
  #  if(length(tester) > 0) {
  #    logger = c(logger, which(links$magnitude == max(links$magnitude[links$from == counter]) & links$from == counter))
  #    counter = links$to[tail(logger, 1)]
  #  } else{
  #    break
  #  }
  #}
  
  #if(eff != "ratio") {
    #E(net_test)$color[logger] = "orange"
  #}

  #layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))
  
  png(filename = paste("traj_", eff, ".png", sep=""), width = 790, height = 840, units = "px")
  
  if(length(names(simplex_chart)) - 1 > 5) {
    if(eff == "ratio") {
      
      ratio_effects = as.character(signif(codes$effect, 2))
      
      ratio_effects[no_imp] = '-'
      ratio_effects[is.na(ratio_effects)] = '+'
      
      plot(net_test, layout = l, vertex.label = c(ratio_effects, codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.9, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 10)
    } else {
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.9, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 10)
      
    }
  } else {
    if(eff == "ratio") {
      
      ratio_effects = as.character(signif(codes$effect, 2))
      
      ratio_effects[no_imp] = '-'
      ratio_effects[is.na(ratio_effects)] = '+'
      
      plot(net_test, layout = l, vertex.label = c(ratio_effects, codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])))
    } else {
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])))
      
    }
  }

  dev.off()
  
  #legend_image = as.raster(rbPal(100), ncol=1)
  #plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '')
  #text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5))
  #rasterImage(legend_image, 0, 0, 1,1)
  
  if(eff == "observed") {
    write.table(paste("There are", length(paths_pos), "paths possible"), "network.txt", row.names = F)
  }
}

## Analysis

outputs = data_loading(error)
epistasis_analysis(p_on, T)
blanket_analysis()
network_analysis()
network_analysis("WT")
network_analysis("ratio")
