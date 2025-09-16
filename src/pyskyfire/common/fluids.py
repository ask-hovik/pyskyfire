# Class for storing multi component propellants. 
# Feels like I could offload some functionality unto this. 
# Maybe use this class to unify coolprop and cea?

class Fluid:
    def __init__(self, type, propellants, fractions): 
        self.type = type
        self.propellants = propellants
        self.fractions = fractions
         


    