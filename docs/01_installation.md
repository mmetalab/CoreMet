# Installation

This file will guide you in installing CoreMet on your computer.

If you are **experienced** with Python, you can just follow the steps in a short-form [down here](#short-form).

## Step 1 - Install Anaconda

## Step 2 - Download this repository

## Step 3 - Install the requirements to run CoreMet

## Step 4 - Run the server

## Step 5 - Access CoreMet in your browser

---

## Short form (for experienced users)

- Clone this repo in a folder of your choice on your pc:  
`cd full/path/to/folder/`  
`git clone https://github.com/mmetalab/mpi-vgae-web.git`
- Go inside the newly created folder:  
`cd ./mpi-vgae-web`
- Create a new conda environment to keep things tidy:  
`conda create -n cormet Python==3.10`
- Activate the newly created environment:  
`conda activate cormet`
- Install the required packages:  
`pip install -r requirements.txt`
- Run the server:  
`python run.py`
- The script will create a local server on `http://127.0.0.1:8080/`. Open that address with any browser.

After the first installation, you can run CoreMet by:

```bash
conda activate cormet
cd full/path/to/mpi-vgae-web/
python run.py
```
