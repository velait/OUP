size <- 10

degree_dist <- function(X) {
  
  dist <- rep(0, nrow(X)) 
  names(dist) <- 1:nrow(X) %>% as.character()
  
  rsums <- rowSums(X)
  
  
  dist[rsums %>% table %>% names] <- rsums %>% table %>% as.vector()
  
  
  dist
  
}

joint_degree_distribution <- function(X) {
  
  dist <- rep(0, nrow(X)) 
  names(dist) <- 1:nrow(X) %>% as.character()
  
  rsums <- rowSums(X)
  
  
  dist[rsums %>% table %>% names] <- rsums %>% table %>% as.vector()
  
  
  dist
  
}

adjacency <- function(x) {
  
  n <- length(x)
  
  f <- matrix(0, n, n)
  b <- matrix(0, n, n)
  
  # Forward
  for(i in 1:(n-1)) {
    
    larger <- x[(i+1):n] >= x[i]
    ind <- i + which(larger)[1]
    
    if(length(ind) != 0) {
      f[i, ind] <- 1
    }
    
  }
  
  # Backward
  for(i in n:2) {
    
    larger <- x[1:(i-1)] >= x[i]
    ind <- which(larger)[length(which(larger))]
    
    if(length(ind) != 0) {
      b[i, ind] <- 1
    }
  }
  
  f + t(f) + b + t(b)
  
  
}



x <- rnorm(10, 0, 1)
y <- rnorm(10, 0, 1)

A <- adjacency(x)
B <- adjacency(y)


joint_degree_distribution <- function(A, B) {
  
  freq <- cbind(rowSums(A), rowSums(B))
  
  df <- matrix(0, nrow(freq), nrow(freq))
  
  for(i in 1:nrow(freq)) {
    ind <- freq[i, ]
    df[ind[1], ind[2]] <- df[ind[1], ind[2]] + 1
  }
  
  
  return(df)
}

mutual_information <- function(A, B) {
  
  MI <- c(0)
  
  joint <- joint_degree_distribution(A, B)
  
  dA <- degree_dist(A)
  dB <- degree_dist(B)
  
  kA <- rowSums(A) %>% unique() %>% sort()
  kB <- rowSums(B) %>% unique() %>% sort()
  
  for(i in kA) {
    for(j in kB) {
      
      mi <- log(joint[i, j]^joint[i, j]/(dA[as.character(i)]*dB[as.character(j)]))
      
      MI <- c(MI, mi)
      
    }
  }
  
  
  return(sum(MI))
}



## Generate random and correlated series

df <- matrix(rnorm(1000, 0, 2), 10, 100)

df <- rbind(df, df[10, ] + rnorm(100, 0, 1))


cal_A <- lapply(1:nrow(df), FUN = function(i) {
  
  df[i, ] %>% adjacency()
  
})


MI_df <- lapply(1:length(cal_A), function(i) {
  
  lapply(1:length(cal_A), function(j) {
    
    mutual_information(cal_A[[i]], cal_A[[j]])
    
  }) %>% unlist()
}) %>% do.call(rbind, .)


## Try on real data

load("/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-David2014/David_phyloseq.Rdata")

# top 20 genera
seqA.genus <- seqA %>% aggregate_taxa("Genus")

seqA.genus.core <- seqA.genus %>% subset_taxa(Genus %in% top_taxa(seqA.genus, 50))

seqA.genus.top20 <- abundances(seqA.genus.core)



# Sliding window,
window_l <- 20
windows <- c(1, 10*(1:20))

sliding_MI_set <- lapply(10*1:10, function(s) {
    print(s)
    sliding_MI <- lapply(windows, function(w) {
      
      
      # Adjacency matrices
      cal_A <- lapply(1:nrow(seqA.genus.top20), FUN = function(i) {
        
        seqA.genus.top20[i, w:(w + s)] %>% adjacency()
        
      })
      
      
      # Mutual information
      MI_df <- lapply(1:length(cal_A), function(i) {
        
        lapply(1:length(cal_A), function(j) {
          
          mutual_information(cal_A[[i]], cal_A[[j]])
          
        }) %>% unlist()
      }) %>% do.call(rbind, .)
      
      
      # compute average MI
      MI_df <- MI_df - diag(diag(MI_df))
      
      MI_df_vec <- upper.triangle(MI_df) %>% as.vector()
      
      mean_MI <- MI_df_vec[MI_df_vec != 0] %>% mean
      
      return(mean_MI)
      
    })
    
    
    df <- cbind(MI = sliding_MI %>% unlist,
                time = windows, w = s) %>%
      as.data.frame()
    return(df)
    
})




sliding_MI_set %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = time, y = MI, group = w, color = as.factor(w)))



