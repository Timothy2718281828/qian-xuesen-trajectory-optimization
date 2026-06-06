# ============================================================
#  Makefile - Final Project: Qian Xuesen Extended Trajectory
# ============================================================
#  Usage:
#    make all    - Generate ALL outputs (plots + MP4 + report PDF)
#    make pdf    - Compile report.tex only (faster)
#    make plots  - Generate all figures from Python scripts
#    make video  - Generate trajectory animation MP4
#    make clean  - Remove ALL generated files
# ============================================================

PYTHON := python
LATEX  := xelatex
FFMPEG := ffmpeg

SRC_DIR   := src
DATA_DIR  := data
BUILD_DIR := build

# All source scripts
PLOT_SCRIPTS := $(SRC_DIR)/visualize.py

.PHONY: all pdf plots video clean

# ============================================================
#  Default target
# ============================================================
all: plots video pdf
	@echo "=== make all complete ==="

# ============================================================
#  Compile LaTeX report
# ============================================================
pdf: report.tex
	@echo "Compiling report.tex with XeLaTeX..."
	mkdir -p $(BUILD_DIR)
	$(LATEX) -output-directory=$(BUILD_DIR) -interaction=nonstopmode report.tex
	$(LATEX) -output-directory=$(BUILD_DIR) -interaction=nonstopmode report.tex
	cp $(BUILD_DIR)/report.pdf .
	@echo "report.pdf generated."

# ============================================================
#  Generate all figures
# ============================================================
plots:
	@echo "Generating figures..."
	$(PYTHON) $(SRC_DIR)/visualize.py
	@echo "Figures generated."

# ============================================================
#  Generate trajectory animation
# ============================================================
video:
	@echo "Generating trajectory animation..."
	$(PYTHON) $(SRC_DIR)/visualize.py --video
	@echo "Animation generated."

# ============================================================
#  Clean
# ============================================================
clean:
	@echo "Cleaning..."
	rm -rf $(BUILD_DIR)
	rm -f report.pdf
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.synctex.gz *.fls *.fdb_latexmk *.xdv
	rm -f *.mp4 *.gif
	rm -rf __pycache__ $(SRC_DIR)/__pycache__
	@echo "Clean complete."
