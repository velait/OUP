# investigate prior

x_seq <- -0:100/100

data.frame(x = x_seq,
           rw = dnorm(x_seq, 1, .1),
           mixture = .75*dnorm(x_seq, 1, .25) + .25*dnorm(x_seq, 0, .25)) %>% 
  mutate(prior = rw*mixture) %>% 
  mutate(prior = prior/sum(prior),
         rw = rw/sum(rw),
         mixture = mixture/sum(mixture)) %>% 
  ggplot() +
  geom_line(aes(x = x, y = rw)) +
  geom_line(aes(x = x, y = mixture), linetype = "dashed") +
  geom_line(aes(x = x, y = prior), color = "red")


data.frame(x = rep(x_seq, 3),
           y = c(dnorm(x_seq, 0, .1), dnorm(x_seq, 0, .2), dnorm(x_seq, 0, .3)), 
           time = rep(1:3, each = length(x_seq)) %>% as.factor()) %>% 
  group_by(time) %>% mutate(y = y/sum(y)) %>% 
  ggplot(aes(x = x, y= y, color = time)) +
  geom_line()
  
           