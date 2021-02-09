# From Sueoka 1988 (https://www.pnas.org/content/pnas/85/8/2653.full.pdf), 
# after Sueoka 1962 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC220819/pdf/pnas00227-0104.pdf)
equilibrium GC content, p = v / (u + v), or
v / u = p / (1 - p)
u / v = (1 - p) / p
k = (1 - p) / p

# according to assaf et al (2017), p = 0.23
# so
k = (1 - 0.23) / 0.23
k = 0.77 / 0.23
k = 3.35

# upper and lower bounds of k according to assaf et al:
upper k = 0.79 / 0.21 = 3.76
lower k = 0.75 / 0.25 = 3
  
  
# old stuff
Lgc/Lat * k = 1
(1-Lat) / Lat = 1/k
k(1-Lat) = Lat
k(1-0.77) = 0.77
0.23k = 0.77
k = 3.35

