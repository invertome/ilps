#!/usr/bin/env python3
# Purpose: Flask web app to explore pipeline results interactively
# Usage: python results_explorer.py --port 5000
# Notes: Serves HTML pages with embedded visualizations (trees, alignments, 3D structures)

from flask import Flask, render_template, request
import os
import argparse

app = Flask(__name__, template_folder='output/templates', static_folder='output')

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Results Explorer Web App')
parser.add_argument('--port', type=int, default=5000, help='Port to run the Flask app on')
args = parser.parse_args()

# Route for the main page
@app.route('/')
def index():
    types = ['prepro', 'pro', 'mature']
    return render_template('index.html', types=types)

# Route for type-specific pages
@app.route('/<type>')
def type_page(type):
    if type not in ['prepro', 'pro', 'mature']:
        return "Invalid type", 404
    # List of available files
    tables = [f for f in os.listdir(f'output/{type}') if f.endswith('.csv')]
    plots = [f for f in os.listdir(f'output/{type}') if f.endswith('.png')]
    trees = [f for f in os.listdir('analysis') if f.startswith(f'{type}_') and f.endswith('.tre')]
    alignments = [f for f in os.listdir(f'analysis/aligned_{type}') if f.endswith('.fasta')]
    structures = [f for f in os.listdir('output/figures') if f.startswith(f'{type}_') and f.endswith('.html')]
    return render_template('type.html', type=type, tables=tables, plots=plots, trees=trees, alignments=alignments, structures=structures)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=args.port)
