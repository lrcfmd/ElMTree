from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField, RadioField
from wtforms.validators import DataRequired
# from wtforms.fields import RadioField

class SearchForm(FlaskForm):
    search_term = StringField("SearchTerm")
    submit = SubmitField("Advanced Search")
    metric_choice = RadioField("Elemental Scale", choices=["mod_petti", "atomic", "mendeleev", "petti", "oliynyk", "olinyk_sc", "jarvis", "jarvis_sc", "magpie", "magpie_sc", "cgcnn", "elemnet", "mat2vec", "matscholar", "megnet16", "random_200"], default="mod_petti")
    comp1 = StringField("Composition 1")
    comp2 = StringField("Composition 2")
    structural = BooleanField()
    experimental = BooleanField()
    must_include = StringField("Elements separated by commans")
    must_exclude = StringField("Elements separated by commas")
