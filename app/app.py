from flask import Flask, render_template, request, redirect
import yaml
import os

app = Flask(__name__)

CONFIG_FILE = 'config.yaml'

@app.route('/', methods=['GET', 'POST'])
def edit_config():
    if request.method == 'POST':
        config_data = {}
        for key in request.form:
            value = request.form[key]
            if key in ['samples', 'sc_samples']:
                config_data[key] = [item.strip() for item in value.split(',')]
            else:
                config_data[key] = value

        with open(CONFIG_FILE, 'w') as file:
            yaml.dump(config_data, file, default_flow_style=False, sort_keys=False)
        return redirect('/')
    else:
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as file:
                config_data = yaml.safe_load(file)
        else:
            config_data = {}
        return render_template('edit_config.html', config=config_data)

if __name__ == '__main__':
    app.run(debug=True)
