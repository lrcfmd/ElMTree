from flask import Flask
from flask_bootstrap import Bootstrap4
from flask_restful import Resource, Api

app = Flask(__name__, template_folder="templates")

app.config['SECRET_KEY'] = "CHANGETHISTOSOMETHINGSECURE"
app.config['BOOTSTRAP_BOOTSWATCH_THEME'] = 'sandstone'

api = Api(app)
bootstrap = Bootstrap4(app)

from app import routes