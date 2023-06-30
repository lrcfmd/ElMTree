import os
from site import getsitepackages

wd = os.getcwd()
sitepack =getsitepackages()

from app import app

if __name__== "__main__":
    app.run(debug=False)

