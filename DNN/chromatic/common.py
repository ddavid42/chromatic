import path
import sys
 
dev_directory = path.Path(__file__).abspath()
sys.path.append(dev_directory.parent.parent.parent)
 
# importing
from num2 import ChNbr

