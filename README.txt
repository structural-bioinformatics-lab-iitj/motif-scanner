# Motif Scanner — Setup & Launch Guide

A web-based protein structure analysis tool for scanning X***X motifs,
Ω-loops, α-turns, β-turns, and γ-turns.

---

## What You Need

Make sure all of these files are in the same folder before you start:

```
your-folder/
├── app.py
├── demo.pdb
├── launch.bat          (Windows)
├── launch.sh           (Linux)
├── motif_backend.py
├── omega_backend.py
└── turn_backend.py
```

---

## Windows Instructions

### Step 1 — Install Python (one time only)

1. Go to https://www.python.org/downloads/ and download Python
2. Run the installer
3. **Important:** tick the box that says "Add Python to PATH" before clicking Install
4. That's it — you never need to do this again

> Python is the only thing you need to install manually.
> The launcher handles everything else automatically.

### Step 2 — Every time you want to use the app

1. Open your folder
2. Double-click **`launch.bat`**
3. It will automatically install all required packages on first run
4. Your browser will open at the scanner

> **Note:** Do not close the black window while using the app.
> When you are done, close the black window to stop the app.

---

## Linux Instructions

### First time only

Open a terminal in your folder and run this once to give the launcher
permission to execute:

```bash
chmod +x launch.sh
```

### Every time you want to use the app

Open a terminal in your folder and type:

```bash
bash launch.sh
```

Press Enter. The script will:
- Automatically install any missing packages
- Start the app
- Open your browser at http://127.0.0.1:5000

> **Note:** Do not close the terminal while using the app.
> When you are done, press Ctrl + C in the terminal to stop the app.

### Opening a terminal in your folder (Linux tip)

- Open your file manager and navigate to your folder
- Right-click on an empty space inside the folder
- Select **"Open Terminal Here"** if available
- Or open a terminal anywhere and type:
  ```bash
  cd ~/Downloads
  ```
  replacing Downloads with wherever your folder is

---

## Using the App

Once the browser opens:

1. **Upload a file** by dragging and dropping onto the upload area,
   or click to browse. Supported formats:
   - Structure files: `.pdb` `.cif` `.mmcif`
   - Sequence files: `.txt` `.fasta` `.fa` `.faa`

2. Click **Run** to start the analysis

3. Results appear in tabs:
   - **Sequence motifs** — X***X pattern hits
   - **Ω-loops** — structural loop motifs
   - **α/β/γ turns** — turn classifications

4. Click any result row to highlight it in the 3D viewer

5. Click **↗** on the viewer to open the fullscreen 3D explorer

6. Use **↓ CSV** buttons to download results as spreadsheets

### Try the demo

Click **"▶ Try the demo"** below the upload area to load a prebuilt
example structure (7OM2, 684 residues) with a guided tour that explains
each section of the results.

### Definitions

Click **"? definitions"** in the top right of the page to see the exact
criteria and literature sources used for each motif type.

---

## Troubleshooting

**"No module named flask" error**
```bash
python3 -m pip install flask biopython freesasa --break-system-packages
```

**"dssp: command not found" warning (Linux)**
```bash
sudo apt install dssp
```

**Browser does not open automatically**
Manually go to http://127.0.0.1:5000 in your browser while the terminal
is running.

**Permission denied (Linux)**
```bash
sudo bash launch.sh
```

**App worked before but not opening now**
Make sure you are running the launcher from inside the correct folder
where all the files are saved together.

---

## Stopping the App

- **Windows:** Close the black terminal window
- **Linux:** Press `Ctrl + C` in the terminal

---

## Contact

If you encounter issues not covered here, contact the project supervisor.
