# Tiny_decide
a web visualizer to display chemical compounds and associated data

You can use the script ```install_and_run_tiny_decide.bat```on Windows to install and run it

Or you can do manually something like :

1. Clone git

```git clone https://github.com/maximestephan/tiny_decide.git```

if you just want to update with the latest verison :
```git pull```

2. Create environment

    ```python3 -m venv tiny_decideenv```
3. Activate environment
    
    ```source tiny_decideenv/bin/activate```
    
4. Install requirements

```pip install -r requirements.txt```

5. run flask

```flask --app main.py run --debug -p 8080```

6. run tiny_decide locally on your browser

open a browser on http://127.0.0.1:8080

Note : The identifier use in this solution is uniquely generated, and if it were ever to resemble one from a Swiss agrochemical company, the similarity would be an extraordinary and purely coincidental anomaly.