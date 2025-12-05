# Tiny_decide
a web visualizer to display chemical compounds and associated data

1. Clone git

```git clone https://github.com/maximestephan/tiny_decide.git```

2. Create environment

on Windows:
    ```py -m venv tiny_decideenv```

on Linux/macOS:
    ```python3 -m venv tiny_decideenv```

3. Activate environment

on Windows:
    ```.\tiny_decideenv\Scripts\activate```
    
on Linux/macOS: 
    ```./tiny_decideenv/bin/activate```
    
( if you get a permission deny, execute this command 
    ```chmod u+x ./tiny_decideenv/bin/activate```
)

4. Install requirements

```pip install -r requirements.txt```

5. run flask

```flask --app main.py run --debug -p 8080```

6. run tiny_decide locally on your browser

open a browser on http://127.0.0.1:8080




Note: the molecules are the one present on Wikipedia and put it together thanks to https://wikipedia.cheminfo.org/