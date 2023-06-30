import os
from site import getsitepackages

wd = os.getcwd()
sitepack =getsitepackages()

#raise Exception(f"wd: {wd}\nsp: {sitepack}")

from app import app

if __name__== "__main__":
    app.run(debug=True)

