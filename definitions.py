import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
PY_SRC_DIR = os.path.join(ROOT_DIR, "arepo_helper")
DATA_DIR = os.path.join(PY_SRC_DIR, "data")
C_SRC_DIR = os.path.join(PY_SRC_DIR, "libs")
AREPO_SRC_DIR = os.path.join(DATA_DIR, "source")
