#!/bin/bash

# ── Motif Scanner Launcher ─────────────────────────────────────────────────
# Double-click this file (or run it in terminal) to start the scanner.
# It will install everything it needs automatically on first run.
# ──────────────────────────────────────────────────────────────────────────

# Go to the folder this script lives in (so app.py can be found)
cd "$(dirname "$0")"

echo ""
echo "  ┌─────────────────────────────────┐"
echo "  │       Motif Scanner             │"
echo "  │  Starting up, please wait...    │"
echo "  └─────────────────────────────────┘"
echo ""

# ── 1. Check Python is installed ───────────────────────────────────────────
if ! command -v python3 &>/dev/null; then
    echo "  ✗  Python3 not found. Installing..."
    sudo apt-get update -qq
    sudo apt-get install -y python3 python3-pip
fi

# ── 2. Install Python packages if missing ──────────────────────────────────
echo "  Checking Python packages..."

python3 -c "import flask" 2>/dev/null || {
    echo "  Installing Flask..."
    python3 -m pip install flask --quiet --break-system-packages 2>/dev/null \
    || python3 -m pip install flask --quiet --user
}

python3 -c "import Bio" 2>/dev/null || {
    echo "  Installing BioPython..."
    python3 -m pip install biopython --quiet --break-system-packages 2>/dev/null \
    || python3 -m pip install biopython --quiet --user
}

python3 -c "import freesasa" 2>/dev/null || {
    echo "  Installing FreeSASA..."
    python3 -m pip install freesasa --quiet --break-system-packages 2>/dev/null \
    || python3 -m pip install freesasa --quiet --user
}

# ── 3. Install DSSP if missing ─────────────────────────────────────────────
if ! command -v dssp &>/dev/null && ! command -v mkdssp &>/dev/null; then
    echo "  Installing DSSP (needs your password)..."
    sudo apt-get install -y dssp 2>/dev/null || true
fi

# ── 4. Check app.py exists ─────────────────────────────────────────────────
if [ ! -f "app.py" ]; then
    echo ""
    echo "  ✗  ERROR: app.py not found."
    echo "     Make sure all files are in the same folder as this launcher."
    echo ""
    read -p "  Press Enter to close..."
    exit 1
fi

# ── 4b. Check .env exists (needed for the "email results" feature) ────────
if [ ! -f ".env" ]; then
    echo ""
    echo "  ⚠  No .env file found -- Gmail 'email results' feature won't work"
    echo "     on this machine until it's set up (this is per-machine, it"
    echo "     doesn't come with the code)."
    echo ""
    read -p "  Set it up now? [Y/n] " setup_now
    if [[ "$setup_now" != "n" && "$setup_now" != "N" ]]; then
        ./setup_env.sh
    else
        echo "  Skipping -- run ./setup_env.sh later to enable email results."
    fi
fi

# ── 5. Free port 5000 if something is already using it ────────────────────
PORT=5000

if command -v systemctl &>/dev/null && systemctl is-active --quiet motifscanner 2>/dev/null; then
    echo "  Stopping motifscanner service to free port $PORT..."
    sudo systemctl stop motifscanner
    sleep 1
fi

PIDS=$(lsof -ti tcp:$PORT 2>/dev/null || fuser $PORT/tcp 2>/dev/null)
if [ -n "$PIDS" ]; then
    echo "  Port $PORT still in use by PID(s): $PIDS -- stopping..."
    kill -9 $PIDS 2>/dev/null
    sleep 1
fi

# ── 6. Start the app ───────────────────────────────────────────────────────
echo "  All good! Starting Motif Scanner..."
echo ""
echo "  ► Opening at http://127.0.0.1:5000"
echo "  ► Close this window to stop the app."
echo ""

# Open browser after a short delay (gives Flask time to start)
(sleep 2 && xdg-open "http://127.0.0.1:5000") &

# Run the app
python3 app.py
