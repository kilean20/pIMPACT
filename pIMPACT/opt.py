import string
import random
import os
from differential_evolution import differential_evolution
from multi_minimization import multi_minimization
        
def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))