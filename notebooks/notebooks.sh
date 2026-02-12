#!/bin/bash
# Unified notebook management tool for UCLCHEM
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

cd "$(dirname "$0")"

check_deps() {
    command -v jupytext &>/dev/null || { echo -e "${RED}Error: jupytext not found${NC}"; exit 1; }
}

show_help() {
    echo -e "${BLUE}Usage:${NC} ./notebooks.sh [--sync|--check|--fix|--generate]"
    echo ""
    echo "  --sync      Full two-way synchronization"
    echo "  --check     Check no .ipynb files are staged"
    echo "  --fix       Convert all .ipynb → .py (creates pairing if needed)"
    echo "  --generate  Generate .ipynb from .py"
}

cmd_sync() {
    echo -e "${BLUE}Full sync...${NC}"
    jupytext --sync *.ipynb 2>/dev/null || true
    echo -e "${GREEN}✓ Done${NC}"
}

cmd_check() {
    echo -e "${BLUE}Checking for staged .ipynb...${NC}"
    staged=$(git diff --cached --name-only --diff-filter=ACM | grep 'notebooks/.*\.ipynb$' | grep -v '.ipynb_checkpoints' || true)
    if [ -n "$staged" ]; then
        echo -e "${RED}✗ .ipynb files should not be committed${NC}"
        echo "$staged"
        echo -e "${YELLOW}Fix with: ./notebooks.sh --fix${NC}"
        exit 1
    fi
    echo -e "${GREEN}✓ No .ipynb files staged${NC}"
}

cmd_fix() {
    echo -e "${BLUE}Converting .ipynb → .py...${NC}"
    for nb in *.ipynb; do
        [ -f "$nb" ] || continue
        py="${nb%.ipynb}.py"
        
        # Check if already paired
        if [ -f "$py" ] && grep -q "jupyter:" "$py" 2>/dev/null; then
            # Already paired, just sync
            jupytext --sync "$nb" && echo -e "${GREEN}✓ $nb (synced)${NC}"
        else
            # New notebook, create pairing
            jupytext --set-formats ipynb,py:percent "$nb" && echo -e "${GREEN}✓ $nb (paired)${NC}"
        fi
    done
    echo -e "${GREEN}✓ All converted${NC}"
}

cmd_generate() {
    echo -e "${BLUE}Generating .ipynb from .py...${NC}"
    for py in *.py; do
        grep -q "jupyter:" "$py" 2>/dev/null || continue
        nb="${py%.py}.ipynb"
        jupytext --to notebook "$py" && echo -e "${GREEN}✓ $nb${NC}"
    done
    echo -e "${GREEN}✓ Done${NC}"
}

check_deps

case "${1:---help}" in
    --sync) cmd_sync ;;
    --check) cmd_check ;;
    --fix) cmd_fix ;;
    --generate) cmd_generate ;;
    --help|-h) show_help ;;
    *) echo -e "${RED}Unknown: $1${NC}"; show_help; exit 1 ;;
esac
