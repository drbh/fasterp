#!/bin/bash
set -e

echo "🦀 Building WASM..."
cargo build --target wasm32-unknown-unknown --no-default-features --features wasm --lib --release

echo "🔧 Generating JS bindings..."
wasm-bindgen target/wasm32-unknown-unknown/release/fasterp.wasm \
  --out-dir pkg \
  --target web \
  --typescript

echo "📦 Copying to book playground..."
mkdir -p book/src/playground/pkg
cp -r pkg/* book/src/playground/pkg/

echo "📚 Building mdBook..."
cd book
mdbook build

echo "✅ Done! Playground available at:"
echo "   - Local: book/book/playground/index.html"
echo "   - Serve: mdbook serve --open (from book/ directory)"
