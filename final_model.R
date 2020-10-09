### Inputs ###

directory = "C:/Users/karol/OneDrive/Program/R/Epistasis Model/weinrich_et_al/p_removed"
file = "weinrich_et_al.csv"
error = F
p_on = F
sign = T

### Library ###

library(gtools)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(igraph) ## Netowrk library
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
    replicates = log10(replicates/as.numeric(d1[which(d1[,1] == d[1,1]), 2]))
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
  
  ggplot(pos_out[-1,], aes(x = indices, y = effect, fill = order)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=formatC(effect, digits=2), y = effect + 0.1*sign(effect)*abs(max(effect))), size = 3) +
    facet_wrap(~ order, scales = "free_x") + 
    labs(x = "Mutation(s)", y = "Fold Effect on Activity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  ggsave("epistatic_effect.png", units = "in", width = 12, height = 9)
  
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
  
  new_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 
  
  new_codes$positions = factor(new_codes$positions, levels = unique(new_codes$positions))
  
  ggplot(new_codes, aes(x = positions, y = effects, color = positions)) +
    geom_boxplot(color = "darkgray") +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_jitter(position = position_jitter(0)) +
    stat_summary(fun=mean, colour="darkred", geom="crossbar", width=0.2) + 
    geom_text_repel(aes(label = identity), force = 0.5, nudge_x = 0.3, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5) +
    labs(x = "Mutant Positions", y = "Mean Fold Effect Change on Activity") +
    theme_classic()
  
  ggsave("mean_fold_effect.png", units = "in", width = 12, height = 9)

  if(sign == T) {
    ## Genotype predicted out
    
    effects_vector = codes %>%
      group_by_at(names(codes)[-dim(codes)[2]]) %>%
      summarise(avg = mean(effect)) %>%
      pull(avg)
    
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
    
    ## Getting WT Effect
    
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
    
    new_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 
    
    new_codes$positions = factor(new_codes$positions, levels = unique(new_codes$positions))
    
    ### WT vs average effect
    
    wt_effects = new_codes$effect[which(new_codes$identity %in% new_codes$identity[-grep("1", new_codes$identity)])]
    wt_vs_avg = data.frame(indices = pos_out[2:(1+length(wt_effects)),]$indices, average_effect = pos_out[2:(1+length(wt_effects)),]$effect, wt_effect = wt_effects)
    
    wt_vs_avg_plot = data.frame(index = rep(wt_vs_avg$indices, 2), effects = c(wt_vs_avg$average_effect, wt_vs_avg$wt_effect), 
                                type = c(rep("average", dim(wt_vs_avg)[1]), rep("wt", dim(wt_vs_avg)[1])))
    
    ggplot(wt_vs_avg_plot, aes(x = index, y = effects, fill = type)) +
      geom_bar(stat = "identity", position = position_dodge())
    
    wt_codes = codes[-length(codes)] 
    wt_codes[wt_codes == -1] = 0
    wt_codes = unique(wt_codes)
    
    wt_codes = wt_codes %>%
      arrange_at(names(wt_codes))
    
    wt_model = lm(paste("effect ~ (",vars ,")", sep=""), codes)
    
    int = codes %>%
      filter_at(vars(-effect), all_vars(. == -1)) %>%
      summarise(mean_effect = mean(effect)) %>%
      pull(mean_effect)
    
    wt_model$coefficients = c(as.numeric(int), wt_effects)
    names(wt_model$coefficients) = names(lm(paste("effect ~ (",vars ,")", sep=""), codes)$coefficients)
    
    preds = c()
    for(i in 1:dim(wt_codes)[1]) {
      preds = c(preds, predict(wt_model, wt_codes[i,]))
    }
    
    pred_df = cbind(pred_df, preds)
    
    colnames(pred_df) = c(colnames(simplex_chart), "observed effect", "epistatic effect", "average effect", "WT effect")
    
    mutations = c()
    for(i in 1:dim(pred_df)[1]){
      mutations = c(mutations, length(which(pred_df[i,1:dim(codes)[2]-1] == 1)))
    }
    
    pred_df$mutations = as.factor(mutations)
    
    ggplot(pred_df, aes(`epistatic effect`, `observed effect`, color = mutations)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = "dashed") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_point() +
      geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5)
    
    print(summary(lm(`observed effect` ~ `epistatic effect`, data = pred_df)))
    
    ggsave("observed_vs_epistatic.png", units = "in", width = 12, height = 9)
    
    ggplot(pred_df, aes(`average effect`, `observed effect`, color = mutations)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = "dashed") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_point() +
      geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5)
    
    print(summary(lm(`observed effect` ~ `average effect`, data = pred_df)))
    
    ggsave("observed_vs_average.png", units = "in", width = 12, height = 9)
    
    ggplot(pred_df, aes(`WT effect`, `observed effect`, color = mutations)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_abline(slope = 1, intercept = 0 + log10(1.5), linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0 - log10(1.5), linetype = "dashed") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_point() +
      geom_text_repel(aes(label = genotype), force = 1, direction = "y", box.padding = 0, segment.size = 0.2, size = 3.5)
    
    print(summary(lm(`observed effect` ~ `WT effect`, data = pred_df)))
    
    ggsave("observed_vs_wt.png", units = "in", width = 12, height = 9)
    
    pred_df <<- pred_df
  }  
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
  }
  
  pred_df
  
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
    
    if(dif > 0) {
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
  
  net_test = graph.data.frame(links, codes, directed = T)
  
  ## Calculate shortest paths
  
  paths_2 = shortest_paths(net_test, from = codes$id[1], to = codes$id[length(codes$id)], 
                           weights = 1/(E(net_test)$magnitude + max(c(abs(min(E(net_test)$magnitude)), abs(max(E(net_test)$magnitude)))) + 1))
  
  
  ## Creates gradient of colors
  
  rbPal = colorRampPalette(c('blue','white', 'red'))
  
  ## This sets 0 in the center with max and min being equal to pos and neg of the largest value
  idk = max(c(abs(max(codes$effect)), abs(min(codes$effect))))
  col = rbPal(length(seq(-idk, idk, idk/10)))[as.numeric(cut(seq(-idk, idk, idk/10), breaks = length(seq(-idk, idk, idk/10))))]
  
  cols = c()
  for(i in 1:length(codes$effect)) {
    cols = c(cols, col[which.min(abs(seq(-idk, idk, idk/10) - codes$effect[i]))])
  }
  
  codes$colors = cols
  
  V(net_test)$frame.color <- "black"
  V(net_test)$color <- codes$colors
  V(net_test)$size <- 15
  E(net_test)$lty <- c(2, 1)[weights]
  E(net_test)$width <- E(net_test)$weights
  E(net_test)$arrow.mode <- 0
  E(net_test)$color <- "grey"
  #for (p in paths_2$vpath) { E(net_test, path=p)$color <- "red" }
  
  ## Finds the optimal path along each node
  
  ## Append to your optimal path algorithm
  ## If you encounter no good options, move back? This is a hard
  
  logger = c()
  counter = links$from[1]
  
  while(T) {
    if(length(which(links$from == counter)) < 1) {
      break
    }
    if(max(links$magnitude[links$from == counter]) < 0) {
      break
    }
    tester = which(links$magnitude == max(links$magnitude[links$from == counter]) & links$from == counter)
    if(length(tester) > 0) {
      logger = c(logger, which(links$magnitude == max(links$magnitude[links$from == counter]) & links$from == counter))
      counter = links$to[tail(logger, 1)]
    } else{
      break
    }
  }
  
  E(net_test)$color[logger] <- "orange"
  
  plot(net_test, layout = l)
}

## Analysis

outputs = data_loading(error)
epistasis_analysis(p_on, T)
network_analysis()
network_analysis("average")
network_analysis("WT")
