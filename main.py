import pandas as pd
import requests, os
from flask import request, Response, send_file
from flask import Flask, render_template
from flask import redirect, url_for, jsonify , flash, session
from time import sleep
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from json import dumps
from io import StringIO
from datetime import datetime
from werkzeug.utils import secure_filename
import os
import platform
import subprocess
import math
import base64
import uuid
import hashlib
import helpers

app = Flask(__name__)

app.secret_key = "MaUnSe"

@app.route('/')
def index():
    list_string = ""
    file_list = os.listdir('data/')
    for file in file_list:
        list_string += f"<option value='{file}' ondblclick='openDataset(this)'>{file}</option>"
    return render_template('index.html', list_string=list_string)


@app.route("/upload", methods=["POST"])
def upload_file():
    if "file" not in request.files:
        flash("No file part in the request.")
        return redirect(url_for("index"))
    file = request.files["file"]
    if file.filename == "":
        flash("No file selected.")
        return redirect(url_for("index"))
    filename = secure_filename(file.filename)
    save_path = "data/" + filename
    file.save(save_path)
    flash(f"File '{filename}' uploaded successfully.")
    return redirect(url_for("index"))

@app.route("/add_molecule", methods=['POST'])
def add_molecule():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    filename = str(data.get("filename"))
    smiles = str(data.get("smiles")).strip()
    if not smiles:
        return jsonify({"error": "SMILES empty"}), 400
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400
    new_row = {col: None for col in df.columns}
    if "smiles" in df.columns:
        new_row["smiles"] = smiles
    new_row["index"] = -1
    new_row["molecule"] = mol = Chem.MolFromSmiles(smiles, sanitize=True)
    new_row["id"] = helpers.smiles_to_id(smiles)
    new_row["name"] = f"added the {datetime.now().strftime('%d.%m.%Y')}"
    
    descs = helpers.compute_pesticide_descriptors_from_mol(mol)
    for desc_name, value in descs.items():
        if desc_name in new_row:
            new_row[desc_name] = float(value) if value is not None else pd.NA
    
    df.loc[-1] = new_row
    df = df.fillna('')
    df = df.sort_index().reset_index(drop=True)
    df["index"] = df["index"] + 1
    globals()[session_uid] = df 
    redirect_url = url_for("table", session_uid=session_uid, filename=filename)
    return jsonify({"status": "ok", "redirect_url": redirect_url})

@app.route("/remove_molecule", methods=['POST'])
def remove_molecule():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    filename = str(data.get("filename"))
    id = str(data.get("id")).strip()
    if not id:
        return jsonify({"error": "id empty"}), 400
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400
    df = df[df["id"] != id]
    df = df.sort_index().reset_index(drop=True)
    df["index"] = df.index
    globals()[session_uid] = df 
    redirect_url = url_for("table", session_uid=session_uid, filename=filename)
    return jsonify({"status": "ok", "redirect_url": redirect_url})


@app.route("/add_column", methods=['POST'])
def add_column():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    filename = str(data.get("filename"))
    column_name = str(data.get("column_name")).strip()
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400
    df[column_name] = ""
    globals()[session_uid] = df 
    redirect_url = url_for("table", session_uid=session_uid, filename=filename)
    return jsonify({"status": "ok", "redirect_url": redirect_url})

@app.route("/remove_column", methods=['POST'])
def remove_column():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    filename = str(data.get("filename"))
    column_name = str(data.get("column_name")).strip()
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400
    df = df.drop(columns=[column_name])
    globals()[session_uid] = df 
    redirect_url = url_for("table", session_uid=session_uid, filename=filename)
    return jsonify({"status": "ok", "redirect_url": redirect_url})

@app.route("/update_cell", methods=["POST"])
def update_cell():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    row = str(data.get("row"))
    col = str(data.get("column"))
    value = str(data.get("value"))
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400
    if col not in df.columns:
        return jsonify({"status": "error", "message": "row or column invalid"}), 400
    try:
        if pd.api.types.is_numeric_dtype(df[col].dtype):
            if value == "" or value is None:
                df.loc[df["id"] == row, col] = pd.NA
            else:
                df.loc[df["id"] == row, col] = float(value)
        else:
            df.loc[df["id"] == row, col] = value
    except Exception as e:
        app.logger.info("expection when trying to put value in dataframe" , e)
        df.loc[df["id"] == row, col] = value
    globals()[session_uid] = df
    return jsonify({
        "status": "ok",
        "value": value,
    })
    
@app.route("/save_data", methods=["POST"])
def save_data():
    data = request.get_json()
    session_uid = str(data.get("session_uid"))
    filename = data.get("filename")
    if not filename:
        return jsonify({"status": "error", "message": "Filename missing"}), 400

    filename = secure_filename(filename)
    df = globals()[session_uid]
    if df is None:
        return jsonify({"error": "unknow session_uid"}), 400

    csv_path = "data/" + filename
    try:
        #df2 = df.drop(columns=["index"])
        df2 = df.drop(columns=["index","molecule"])
        df2.to_csv(csv_path, index=False, sep="\t")
    except Exception as exc:
        app.logger.info("Error saving CSV:", exc)
        return jsonify({
            "status": "error",
            "message": "Failed to save file"
        }), 500

    return jsonify({
        "status": "ok",
        "filename": filename
    })


@app.route("/export_data/", methods = ['POST'])
def export_data():
    session_uid = str(request.form.get('session_uid'))  
    filename = str(request.form.get('filename'))
    df = globals()[session_uid]
    df2 = df.drop(columns=["index","molecule"])
    buffer = StringIO()
    df2.to_csv(buffer, index=False)
    csv_data = buffer.getvalue()

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename={filename}"
        },
    )

@app.route("/set-significant", methods=["POST"])
def set_significant():
    data = request.get_json()
    value = bool(data.get("significant"))
    session["use_significant"] = value
    return jsonify({"status": "ok", "value": value})

@app.route("/open-data-folder", methods=["GET"])
def open_data_folder():
    folder_path = os.path.abspath("data")
    system = platform.system()
    try:
        if system == "Darwin":  # macOS
            subprocess.Popen(["open", folder_path])
        elif system == "Windows":
            subprocess.Popen(["explorer", folder_path])
        elif system == "Linux":
            subprocess.Popen(["xdg-open", folder_path])
        else:
            return jsonify({"error": "Unsupported OS"}), 500
        return jsonify({"status": "ok"})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@app.route('/load/<filename>')
def load(filename):
    csv_path = "data/" + filename
    session_uid = uuid.uuid4().hex
    if csv_path.endswith('.csv'):
        app.logger.info('Loading csv')
        df = pd.read_csv(csv_path, sep="\t", nrows= 99999)
        #PandasTools.AddmoleculeColumnToFrame(df, 'SMILES', 'molecule')
        df['molecule'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x, sanitize=True) )
        # if molecule is None then molecule become empty molecule
        df['molecule'] = df['molecule'].apply(lambda molecule: molecule if molecule else Chem.MolFromSmiles('') )
    elif csv_path.endswith('.sdf'):
        app.logger.info('Loading sdf')
        mol_supplier = Chem.MultithreadedSDMolSupplier(csv_path, numWriterThreads=5)
        molecules = []
        for mol in mol_supplier:
            if mol is not None:
                props = mol.GetPropsAsDict()
                props["Title"] = mol.GetProp("_Name")
                props["molecule"] = mol
                molecules.append(props)
        df = pd.DataFrame(molecules)

    #if("Name" not in df.columns):
    #    df['Name'] = df['molecule'].apply(lambda x: Chem.InchiToInchiKey(Chem.MolToInchi(x, options='-KET -15T')))

    df.insert(0, 'molecule', df.pop('molecule'))
    df.reset_index(inplace=True)
    globals()[session_uid] = df 
    return redirect(url_for("table", session_uid=session_uid, filename=filename))



@app.route('/table')
def table():
    use_significant = session.get("use_significant", True)
    session_uid = request.args.get("session_uid")
    filename = request.args.get("filename")
    sort_key = request.args.get("sort")
    sort_order = request.args.get("order", "asc")

    df = globals().get(session_uid)
    if df is None:
        return redirect(url_for("load",  filename=filename))
    
    if sort_key is not None and sort_key in df.columns:
        ascending = sort_order != "desc"
        df = df.sort_values(by=sort_key, ascending=ascending)
    table = "<table class='maintable'>"
    #for index, row in df.iterrows():
    for row_number, (index, row) in enumerate(df.iterrows()):
        if (row_number == 0 ):
            table += "<tr>"
            for column in df.columns:
                if column == 'index': 
                    table += "<th  class='stickyHeaderIndex' onclick='toggleContextMenu(event)' title='click to open menu options'>•••</th>"        
                elif column == 'molecule':
                    table += "<th  class='stickyHeaderMolecule'>molecule</th>"     
                elif column == 'id':
                    sortableHeader = f"sortableHeader sorted {sort_order}"  if column == sort_key else "sortableHeader"        
                    table += "<th  class='stickyHeaderName " + sortableHeader +"' data-sort-key='"+ str(column) + "'>id</th>"  
                else:
                    sortableHeader = f"sortableHeader sorted {sort_order}"  if column == sort_key else "sortableHeader"        
                    column_name = str(column).replace('_',' ').replace('*',' ').replace('.','<br>')    
                    table += "<th name='"+ str(column) + "' class='" + sortableHeader +"' data-sort-key='"+ str(column) + "'  title='Click to sort by this column'>" + column_name + "</th>"        
            table += "</tr>"    
        table += f"<tr id='{str(index )}' row_number='{str(row_number+1 )}' class='observable' ></tr>"
    table += "</table>"
    
    check_significant = "checked" if (use_significant) else ""
    return render_template('main.html', table=table, session_uid = session_uid, filename= filename , significant = check_significant)
    


@app.route("/readrow/" , methods = ['POST'])
def readrow():
    session_uid = str(request.form.get('session_uid'))  
    index = int(str(request.form.get('index')))
    row_number = (str(request.form.get('row_number')))
    try:
        df = globals()[session_uid]
        row = df.iloc[index] #.query('index == ' + str(index) )
        return renderRow(row, df.columns, row_number)
    except:
        return Response("", status=204)
    


def renderRow(row, columns, row_number):
    table = "<tr>"
    use_significant = session.get("use_significant", True)
    for column in columns:
        if ( column == "molecule"):
            try:
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 100)
                #drawer.drawOptions().useBWAtomPalette()
                drawer.SetLineWidth(0.9)
                drawer.drawOptions().minFontSize = 9
                color=(245/255.0,245/255.0, 245/255.0, 245/255.0)
                drawer.drawOptions().setBackgroundColour(color)
                rdMolDraw2D.PrepareAndDrawMolecule(drawer, row['molecule'])                
                table += f"""<td class='stickyMolecule'><div style=' height: 100px; width: 200px;'><img src='data:image/png;base64, {base64.b64encode(drawer.GetDrawingText()).decode('utf8')}' /></div></td>"""
            except Exception as ex :
                app.logger.info("An exception occurred " + str(ex))
                table += f"<td class='stickyMolecule'><div style=' display: table-cell;  border: 0px solid black; height: 104px; width: 204px;'><img src='data:image/png;base64,xxx' /></div></td>"
        elif ( column == "smiles"):
            table += "<td class='smiles' >" + str(row[column]) + "</td>"
        elif ( column == "index"):
            table += "<td class='stickyIndex' >" + row_number + "</td>"
        elif ( column == "id"):
            table += "<td class='stickyName'>" + str(row[column]) + "</td>"
        
        else:
            
            if isinstance(row[column],float):
                if math.isnan(row[column]):
                    tdvalue = ""
                else:
                    try:
                        if use_significant:
                            tdvalue = "<span title='" + str(row[column]) + "'>" + str(helpers.round_to_one_significant_decimal(row[column])) + "</span>"
                        else:
                            tdvalue = "<span >" + str(row[column]) + "</span>"

                    except Exception as ex :
                        app.logger.info("An exception occurred " + str(ex))
                        tdvalue = "<span style='color:red' >" + str(row[column]) + "</span>"
            else:
                tdvalue = str(row[column] )
            table += "<td title='double click to edit the value of this cell' ondblclick='editCell(this)' data-row='" + row["id"] + "' data-col='" + column + "' ><span>" + tdvalue + "</span></td>"
    table += "</tr>"
    return table


