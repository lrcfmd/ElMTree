import os
from collections import defaultdict
from time import time
from datetime import datetime

import numpy as np
import pickle as pk
import pandas as pd

from flask import render_template
from flask import request as flask_request

from app import app, api
from app.forms import SearchForm

from ElMD import ElMD

from logging.config import dictConfig

from ElMTree import ElMTree
from ElMTree import Entry 

from jinja2 import BaseLoader, TemplateNotFound, ChoiceLoader
from urllib import request, parse

from flask_restful import Resource, Api, reqparse

import plotly.express as px
from umap import UMAP

import os

print("loading tree")

tree = pk.load(open("indexed_LC.pk", "rb"))
print("Loading lookup")

lookup = pk.load(open("mtree_lookup.pk", "rb"))

n_comps = len(lookup)
n_records = sum([len(lookup[comp][dataset]) if isinstance(lookup[comp][dataset], list) and dataset != "compound_formula" else 1 for comp in lookup.keys() for dataset in lookup[comp].keys()])

class UrlLoader(BaseLoader):
    def __init__(self, url_prefix):
        self.url_prefix = url_prefix

    def get_source(self, environment, template):
        url = parse.urljoin(self.url_prefix, template)
        try:
            t = request.urlopen(url)
            if t.getcode() is None or t.getcode() ==200:
                return t.read().decode('utf-8'), None, None
        except IOError:
            pass

        raise TemplateNotFound(template)

# Add this to ensure can access static files dynamically
app.jinja_loader = ChoiceLoader([app.jinja_loader, UrlLoader('https://lmds.liv.ac.uk/static/')])

dictConfig({"version": 1,
            "disable_existing_loggers": False,
            "formatters": {"default": {
                        "format": '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
                }},
            "handlers": {
                "wsgi": {
                    "class": "logging.StreamHandler",
                    "stream": "ext://flask.logging.wsgi_errors_stream",
                    "formatter": "default",
                    }
                },

            "root": {"level": "DEBUG", "handlers": ["wsgi"]},
            })

@app.route("/search", methods=['GET', 'POST'])
def search():
    form = SearchForm()
    advanced_form = SearchForm()

    if form.validate_on_submit():
        app.logger.debug(form.search_term.data)
        query = str(form.search_term.data)

        comp = ElMD(query, metric="fast").pretty_formula

        if query == "":
            return render_template("search.html", form=form, composition="out func")

        ts = time()
        query = ElMD(query, metric="fast")
        search_result = tree.knn(query, 100)
         
        app.logger.debug(search_result)
        database_ids = [lookup[c[0].pretty_formula] for c in search_result if c[0].pretty_formula in lookup]

        for i in range(len(database_ids)):
            if "mpds" in database_ids[i]:
                database_ids[i]["mpds"] = list(set([x.split("-")[0] for x in database_ids[i]["mpds"]]))

        database_ids = [[f"{k}: {v}" for k, v in comp_dict.items()] for comp_dict in database_ids ]
        
        result = [(x[0][0].pretty_formula, np.round(x[0][1], 3), x[1]) for i, x in enumerate(zip(search_result, database_ids))]
        
        time_taken = time() - ts
	
        return render_template("search.html", search_result=result, form=form, advanced_form=advanced_form, time_taken=time_taken,n_comps=n_comps, n_records=n_records, composition=comp)

    if advanced_form.validate_on_submit():
        query = str(advanced_form.search_term.data)

        comp = ElMD(query, metric="fast").pretty_formula

        if query == "":
            return render_template("search.html", form=form, advanced_form=advanced_form, composition="out func")

        ts = time()
        query = ElMD(query, metric="fast")

        must_include = set(split(advanced_form.must_include.data, ","))
        must_exclude = set(split(advanced_form.must_exclude.data, ","))

        search_dict = {"experimental": advanced_form.experimental.data,
                       "structural": advanced_form.structural.data,
                       "must_include": must_include,
                       "must_exlude": must_exclude}

        search_result = tree.knn(query, 100, advanced_search=search_dict)
         
        database_ids = [lookup[c[0].pretty_formula] if c[0].pretty_formula in lookup else "" for c in search_result] 

        for i in range(len(database_ids)):
            if "mpds" in database_ids[i]:
                database_ids[i]["mpds"] = list(set([x.split("-")[0] for x in database_ids[i]["mpds"]]))


        database_ids = [str(x)[1:-1].replace("'", "").replace("[", "").replace(".cif", "").split("],") for x in database_ids]
        
        result = [(x[0][0].pretty_formula, np.round(x[0][1], 3), x[1]) for i, x in enumerate(zip(search_result, database_ids))]
        

        time_taken = time() - ts

        return render_template("search.html", search_result=result, form=form, advanced_form=advanced_form, time_taken=time_taken,n_comps=n_comps, n_records=n_records, composition=comp)

    return render_template("search.html", form=form, advanced_form=advanced_form, composition="out func")

@app.route("/elm2d", methods=['GET', 'POST'])
def elm2d():
    form = SearchForm()
    
    if form.validate_on_submit():
        query = str(form.search_term.data)
        truncated = False

        if query == "":
            return render_template("elm2d.html", form=form, composition="out func")

        if "," not in query:
            comp = ElMD(query, metric="fast").pretty_formula
            
            ts = time()
            q = ElMD(query, metric="fast")
            search_result = tree.knn(q, 100)
            
            comps_elmd = [q] + [x[0] for x in search_result]
            comps = [x.pretty_formula for x in comps_elmd]

        else:
            if len(comps) < 4:
                n = len(comps)
                comps = comps + comps
                truncated = True

            if len(comps) > 100:
                return "Error: Must enter fewer than 100 compostions"
            
            comps = [x for x in query.split(",")]
            comps_elmd = [ElMD(x) for x in comps]
        
        dm = np.zeros((len(comps), len(comps)))
        
        try:
            for i in range(len(comps)-1):
                for j in range(i + 1, len(comps)):
                    d = comps_elmd[i].elmd(comps_elmd[j])
                    dm[i][j] = d
                    dm[j][i] = d

        except Exception as e:
            app.logger.debug(e)
                
        mapper = UMAP(metric="precomputed")
        X = mapper.fit_transform(dm)
        app.logger.debug(X[:10])

        if truncated:
            df_data = {"x": X[:n, 0],
                       "y": X[:n, 1],
                       "comps": [c for i, c in enumerate(comps[:n])]}
        else:
            df_data = {"x": X[:, 0], 
                  "y": X[:, 1], 
                   "comps": [c for c in comps]}
        
        if "," in query:
            df = pd.DataFrame(df_data)
    
            fig = px.scatter(df, x="x", y="y", hover_data=["comps"], template="simple_white", labels={"x": "", "y": ""})
 
        if "," not in query:
            df_data["colors"] = ["query"] + ["other" for x in range(len(comps) -1)]
            df_data["distance"] = np.round(dm[:, 0], 3)

            df = pd.DataFrame(df_data)

            fig = px.scatter(df, x="x", y="y", hover_data=["comps", "distance"], color="colors", template="simple_white", labels={"x": "", "y": ""})
       
        fig.update_xaxes(showticklabels=False, visible=False)
        fig.update_yaxes(showticklabels=False, visible=False)
        fig.update_traces(marker_size=20)
        plot = fig.to_html()

        fig = px.imshow(dm, 
                        labels=dict(x="Composition 1", 
                                y="Composition 2", 
                                color="Distance"),
                        x=comps, # [x.pretty_formula for x in comps],
                        y=comps) # [x.pretty_formula for x in comps])
        fig.update_xaxes(side="top")
        dm_plot = fig.to_html()

        return render_template("elm2d.html", plot=plot, dm_plot=dm_plot, form=form, dm=dm)
    
    return render_template("elm2d.html", form=form)


class ElMTreeApiEndpoint(Resource):
    def get(self):
        return {"Hello": "world"}

    def put(self):
        query = flask_request.form['query']

        if query[0] == "[":
            query = query[1:]
        if query[-1] == "]":
            query = query[:-1]

        comps = [comp.strip('\"') for comp in query.split(",") if comp != ""]
        app.logger.debug(f"API Search")

        queries = [ElMD(comp, metric="fast").pretty_formula for comp in comps]
        results = []

        ts = time()

        for query in queries:
            query = ElMD(query, metric="fast")
            search_result = tree.knn(query, 10)
             
            database_ids = [lookup[c[0].pretty_formula] if c[0].pretty_formula in lookup else "" for c in search_result] 

            for i in range(len(database_ids)):
                if "mpds" in database_ids[i]:
                    database_ids[i]["mpds"] = list(set([x.split("-")[0] for x in database_ids[i]["mpds"]]))

            results.append([(x[0][0].pretty_formula, np.round(x[0][1], 3), x[1]) for i, x in enumerate(zip(search_result, database_ids))])

            time_taken = time() - ts

        return {"ElMTree Searches": [str(x) for x in results]}
 
api.add_resource(ElMTreeApiEndpoint, "/elmtree_api")

@app.route("/elmd", methods=['GET', 'POST'])
def elmd():
    form = SearchForm()
    
    if form.validate_on_submit():
        app.logger.debug(form.comp1.data)
        app.logger.debug(form.comp2.data)
        metric = form.metric_choice.data

        elmd_comp1 = ElMD(str(form.comp1.data), metric=metric)
        elmd_comp2 = ElMD(str(form.comp2.data), metric=metric)
        
        distance, plan = elmd_comp1.elmd(elmd_comp2, return_assignments=True)
        
        plan = plan.reshape((len(elmd_comp1.composition.keys()), len(elmd_comp2.composition.keys())))

        fig = px.imshow(plan, 
                        labels=dict(x="Composition 2", 
                                y="Composition 1", 
                                color="Mass Transported"),
                        y=[el.split("0")[0] for el in elmd_comp1.pretty_formula.split(" ")],
                        x=[el.split("0")[0] for el in elmd_comp2.pretty_formula.split(" ")])
        fig.update_xaxes(side="top")
        plot = fig.to_html()

        return render_template("elmd.html", distance=distance, form=form, elmd_comp1=elmd_comp1.pretty_formula, elmd_comp2=elmd_comp2.pretty_formula, plot=plot)

    return render_template("elmd.html", form=form)
