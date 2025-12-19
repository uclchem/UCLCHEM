#!/bin/bash
# Development environment setup script for UCLCHEM contributors

set -e

echo "Setting up UCLCHEM development environment..."

# Install package with dev dependencies
echo "Installing package with dev dependencies..."
pip install -e ".[dev]"

# Install pre-commit hooks
echo "Installing pre-commit hooks..."
pre-commit install

# Run pre-commit on all files to ensure everything is formatted
echo "Running pre-commit checks on all files..."
pre-commit run --all-files || true

echo ""
echo "✓ Development environment setup complete!"
echo ""
echo "Pre-commit hooks are now installed and will run automatically before each commit."
echo "You can manually run all checks with: pre-commit run --all-files"
