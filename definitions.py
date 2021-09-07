import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
PY_SRC_DIR = os.path.join(ROOT_DIR, "arepo_helper")
DATA_DIR = os.path.join(PY_SRC_DIR, "data")
LIBS_DIR = os.path.join(PY_SRC_DIR, "libs")
C_SRC_DIR = os.path.join(LIBS_DIR, "src")
C_INC_DIR = os.path.join(LIBS_DIR, "headers")
AREPO_SRC_DIR = os.path.join(DATA_DIR, "source")
