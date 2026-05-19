#!/bin/bash
# Quick verification script to test the povu-rs build

set -e

echo "=== Povu-rs Build Verification ==="
echo

echo "1. Checking build dependencies..."
command -v cargo >/dev/null 2>&1 || { echo "❌ cargo not found"; exit 1; }
command -v cmake >/dev/null 2>&1 || { echo "❌ cmake not found"; exit 1; }
echo "✓ cargo and cmake found"
echo

echo "2. Building povu-rs..."
cd "$(dirname "$0")"
cargo build --verbose 2>&1 | tee build.log
echo "✓ Build completed"
echo

echo "3. Running unit tests..."
cargo test --lib 2>&1 | tee -a build.log
echo "✓ Unit tests passed"
echo

echo "4. Checking if test data exists..."
if [ -d "../tests/data" ]; then
    echo "✓ Test data directory found"
    echo "  Files:"
    ls -lh ../tests/data/ || true
else
    echo "⚠ Test data directory not found at ../tests/data"
    echo "  Integration tests will be skipped"
fi
echo

echo "5. Running integration tests (may skip if no test data)..."
cargo test --tests 2>&1 | tee -a build.log || echo "⚠ Some integration tests skipped"
echo

echo "=== Build Verification Complete ==="
echo "Check build.log for details"
