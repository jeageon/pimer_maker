#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
VENV_DIR="$PROJECT_DIR/.venv"
cd "$PROJECT_DIR"

if [ ! -d "$VENV_DIR" ]; then
  if [ -x "$(command -v python3)" ]; then
    python3 -m venv "$VENV_DIR"
  else
    echo "python3 not found. Please install Python 3."
    exit 1
  fi
fi

source "$VENV_DIR/bin/activate"

PYTHON_BIN="$VENV_DIR/bin/python3"
if [ ! -x "$PYTHON_BIN" ]; then
  PYTHON_BIN="$VENV_DIR/bin/python"
fi

if [ ! -x "$PYTHON_BIN" ]; then
  echo "Python executable not found in virtual environment."
  exit 1
fi

$PYTHON_BIN -m pip install -r requirements.txt
STREAMLIT_BROWSER_GATHER_USAGE_STATS=false "$PYTHON_BIN" -m streamlit run src/webui.py
