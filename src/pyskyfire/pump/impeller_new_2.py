# Just a simple calculator for now. Will be expanded later

n = 3000 #rpm
Q = 0.5
H = 20

beta_1 = 15 #check
beta_2 = 22.5

g = 9.81

# Specific speed
n_s = n*Q**(1/2)/(g*H)**(3/4)
Z = int(round(beta_2/3, 0))

#sigma = 1 - 

